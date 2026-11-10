# =============================================================================
# 全量双向 MR-PRESSO 分析（128 组）
# -----------------------------------------------------------------------------
# 设计：
#   2 个本地生殖表型（MI / AS，FinnGen gz） × 32 个癌症 GWAS（OpenGWAS 在线）
#   × 2 个方向（正向：生殖表型 -> 癌症；反向：癌症 -> 生殖表型）= 128 组
#
# 可复现性：
#   * 每组分析前 set.seed(1234)，MR-PRESSO 随机模拟完全可复现
#   * 工具变量 / QC 参数与手稿主分析完全一致：
#       p1 = 1e-5；LD clump r2 = 0.001, kb = 10000；
#       harmonise action = 2；剔除 palindromic SNP；剔除 MAF < 0.01
#   * MR-PRESSO：NbDistribution = 1000，SignifThreshold = 0.05
#
# 兼容性：
#   * 全部路径使用正斜杠 "/"（Windows 下 R 可直接识别）
#   * 任何单组失败（SNP 不足 / API 报错 / 包内错误）均自动跳过并记录，
#     不会中断整体循环
#
# 依赖：TwoSampleMR、MRPRESSO、data.table、ieugwasr、dplyr（tidyverse 系）
# =============================================================================

# -----------------------------------------------------------------------------
# 0. 基础环境
# -----------------------------------------------------------------------------
options(stringsAsFactors = FALSE)

library(data.table)
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(MRPRESSO)

# ---- 0.1 OpenGWAS 令牌（新版 ieugwasr 通过环境变量 OPENGWAS_JWT 自动携带）----
# OpenGWAS 令牌：从环境变量 OPENGWAS_JWT 读取，严禁在源码中硬编码 (见 .gitignore)
my_token <- Sys.getenv("OPENGWAS_JWT")
if (nchar(my_token) == 0) {
  stop("环境变量 OPENGWAS_JWT 未设置，请先 Sys.setenv(OPENGWAS_JWT=\"<your token>\") 或写入 .Renviron")
}
# ieugwasr 新版会自动读取 OPENGWAS_JWT 环境变量；此处仅确保其已设置
Sys.setenv(OPENGWAS_JWT = my_token)
# ---- 0.2 固定分析参数（与手稿一致，请勿改动）----
p1_thr       <- 1e-5    # 工具变量显著性阈值
clump_r2     <- 0.001   # LD 聚类 r2
clump_kb     <- 10000   # LD 聚类窗口 kb
maf_thr      <- 0.01    # 次要等位基因频率下限
use_proxies  <- FALSE   # 不使用 LD 代理
nb_dist      <- 1000    # MR-PRESSO 模拟次数
signif_thr   <- 0.05    # MR-PRESSO 离群判定阈值
set_seed     <- 1234    # 全局随机种子

# ---- 0.3 路径与输出目录 ----
script_dir <- "D:/研究生科研小组/小论文/meta+MR/MR/MR-PRESSO"
out_dir    <- file.path(script_dir, "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
final_csv  <- "D:/研究生科研小组/小论文/MR_PRESSO_full_results.csv"   # 任务指定导出路径
log_file   <- file.path(out_dir, "run_log.txt")                      # 详细运行日志
log_msg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

# -----------------------------------------------------------------------------
# 1. 32 个癌症数据集清单（ID 与补充表 S5、S8–S11 100% 对应）
# -----------------------------------------------------------------------------
cancer_map <- data.frame(
  cancer_subgroup = c(
    rep("Colon cancer", 3),
    rep("Colorectal cancer", 4),
    rep("Liver cancer", 3),
    rep("Oesophagus cancer", 1),
    rep("Gastric cancer", 2),
    rep("Bladder cancer", 3),
    rep("Prostate cancer", 2),
    rep("Testicular cancer", 1),
    rep("Leukaemia", 4),
    rep("Non-Hodgkin lymphoma", 1),
    rep("Hodgkin lymphoma", 1),
    rep("Melanoma (skin cancer)", 5),
    rep("Thyroid cancer (endocrine)", 2)
  ),
  cancer_id = c(
    "ukb-b-20145", "ukb-d-C3_COLON", "finn-b-C3_COLON",
    "ebi-a-GCST90013862", "ebi-a-GCST90013866", "finn-b-C3_COLORECTAL", "ieu-b-4965",
    "ebi-a-GCST90018638", "ebi-a-GCST90018858", "ieu-b-4953",
    "finn-b-C3_OESOPHAGUS",
    "ebi-a-GCST90018629", "ebi-a-GCST90018849",
    "ukb-b-8193", "ieu-b-4874", "finn-b-C3_BLADDER",
    "finn-b-C3_PROSTATE", "ieu-b-85",
    "finn-b-C3_TESTIS",
    "ieu-b-4914", "finn-b-C3_AML", "ieu-b-4956", "ieu-b-4958",
    "finn-b-CD2_NONHODGKIN_NAS",
    "finn-b-CD2_HODGKIN_LYMPHOMA",
    "finn-b-C3_MELANOMA_SKIN", "ukb-d-C3_MELANOMA_SKIN", "ieu-b-4969",
    "ukb-a-58", "ukb-b-12915",
    "ebi-a-GCST90018929", "ieu-a-1082"
  )
)
stopifnot(nrow(cancer_map) == 32)

# -----------------------------------------------------------------------------
# 2. 两个本地 FinnGen 生殖表型的读取与标准化（只读一次，缓存复用）
# -----------------------------------------------------------------------------
# FinnGen 原始列：#chrom, pos, ref, alt, rsids, beta, sebeta, pval,
#                 频率列 R12 为 af_alt，R5 为 maf（均为 alt/效应等位频率）
trait_map <- data.frame(
  trait_name = c("Male infertility", "Abnormal spermatozoa"),
  trait_id   = c("finngen_R12_N14_MALEINFERT", "finngen_R5_R18_ABNORMAL_SPERMATOZ"),
  file_path  = c(
    "D:/研究生科研小组/小论文/GWAS数据库/finn数据/finngen_R12_N14_MALEINFERT.gz",
    "D:/研究生科研小组/小论文/GWAS数据库/finn数据/finngen_R5_R18_ABNORMAL_SPERMATOZ.gz"
  )
)

# 缓存环境：避免在函数内重复读取大文件（R12 约 768MB，R5 约 551MB）
cache_env <- new.env(parent = emptyenv())

read_finngen_std <- function(trait_id) {
  if (!is.null(cache_env[[trait_id]])) return(cache_env[[trait_id]])

  path <- trait_map$file_path[trait_map$trait_id == trait_id]
  log_msg("首次读取 FinnGen 文件：", path)

  # 先读表头，自动识别频率列（af_alt 或 maf）
  header <- names(fread(path, nrows = 0, showProgress = FALSE))
  eaf_col <- if ("af_alt" %in% header) "af_alt" else if ("maf" %in% header) "maf" else
    stop("未找到频率列 af_alt/maf：", path)

  fg <- fread(path,
              select = c("rsids", "ref", "alt", "pval", "beta", "sebeta", eaf_col),
              showProgress = FALSE)
  setnames(fg, eaf_col, "eaf")

  # rsids 清洗：多 rsID 取第一个；剔除缺失 / 占位符；去重
  first_rs <- trimws(vapply(strsplit(fg$rsids, "[,;]"),
                            function(z) if (length(z) >= 1) z[1] else NA_character_,
                            character(1)))
  fg$SNP <- first_rs
  fg <- fg[!is.na(fg$SNP) & nzchar(fg$SNP) & grepl("^rs[0-9]+$", fg$SNP), ]
  fg <- fg[!duplicated(fg$SNP), ]

  # 标准化为 TwoSampleMR 兼容字段
  std <- data.frame(
    SNP           = fg$SNP,
    beta          = fg$beta,
    se            = fg$sebeta,
    pval          = fg$pval,
    effect_allele = fg$alt,    # alt 为效应等位
    other_allele  = fg$ref,    # ref 为其他等位
    eaf           = fg$eaf,
    stringsAsFactors = FALSE
  )
  cache_env[[trait_id]] <- std
  log_msg(trait_id, " 标准化完成，可用 SNP 数：", nrow(std))
  std
}

# 取“作为暴露”的本地 FinnGen 工具变量（p1 过滤 + 在线 LD 聚类）
get_local_exposure <- function(trait_id, trait_name) {
  std <- read_finngen_std(trait_id)
  exp <- format_data(std, type = "exposure",
                     snp_col = "SNP", beta_col = "beta", se_col = "se",
                     effect_allele_col = "effect_allele",
                     other_allele_col = "other_allele",
                     pval_col = "pval", eaf_col = "eaf")
  exp <- exp[exp$pval.exposure < p1_thr, ]
  if (nrow(exp) == 0) return(exp)
  exp <- clump_data(exp, clump_r2 = clump_r2, clump_kb = clump_kb)
  exp$id.exposure <- trait_id
  exp$exposure    <- trait_name
  exp
}

# 取“作为结局”的本地 FinnGen 数据（按暴露 SNP 子集提取）
get_local_outcome <- function(trait_id, trait_name, snps) {
  std <- read_finngen_std(trait_id)
  sub <- std[std$SNP %in% snps, ]
  out <- format_data(sub, type = "outcome",
                     snp_col = "SNP", beta_col = "beta", se_col = "se",
                     effect_allele_col = "effect_allele",
                     other_allele_col = "other_allele",
                     pval_col = "pval", eaf_col = "eaf")
  out$id.outcome <- trait_id
  out$outcome    <- trait_name
  out
}

# -----------------------------------------------------------------------------
# 3. 在线癌症数据：暴露（工具变量）与结局（按 SNP 提取）
# -----------------------------------------------------------------------------
get_cancer_exposure <- function(cid) {
  key <- paste0("exp_", cid)
  if (!is.null(cache_env[[key]])) return(cache_env[[key]])
  exp <- extract_instruments(outcomes = cid, p1 = p1_thr, clump = TRUE,
                             r2 = clump_r2, kb = clump_kb)
  cache_env[[key]] <- exp
  exp
}

get_cancer_outcome <- function(cid, snps) {
  extract_outcome_data(snps = snps, outcomes = cid, proxies = use_proxies)
}

# -----------------------------------------------------------------------------
# 4. 统一 QC：harmonise(action=2) + 去 palindromic + 去 MAF<0.01
# -----------------------------------------------------------------------------
harmonise_qc <- function(exp_dat, out_dat) {
  har <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat, action = 2)
  har <- har[har$mr_keep == TRUE, ]
  if (nrow(har) == 0) return(har)

  if ("palindromic" %in% names(har)) har <- har[!har$palindromic, ]
  if (nrow(har) == 0) return(har)

  if (!"eaf.exposure" %in% names(har)) har$eaf.exposure <- NA_real_
  if (!"eaf.outcome" %in% names(har))  har$eaf.outcome  <- NA_real_
  maf_e <- pmin(har$eaf.exposure, 1 - har$eaf.exposure)
  maf_o <- pmin(har$eaf.outcome,  1 - har$eaf.outcome)
  maf_min <- apply(cbind(maf_e, maf_o), 1, function(z) {
    z <- z[!is.na(z)]
    if (length(z) == 0) NA_real_ else min(z)
  })
  har[!is.na(maf_min) & maf_min >= maf_thr, ]
}

# -----------------------------------------------------------------------------
# 5. IVW 计算工具：返回 beta 与格式化 OR(95%CI) 字符串
# -----------------------------------------------------------------------------
ivw_summary <- function(har) {
  ivw <- mr(har, method_list = "mr_ivw")
  if (nrow(ivw) == 0) return(list(beta = NA_real_, ci = NA_character_))
  orr <- generate_odds_ratios(ivw)
  list(beta = orr$b[1],
       ci = sprintf("%.2f (%.2f, %.2f)", orr$or[1], orr$or_lci95[1], orr$or_uci95[1]))
}

# -----------------------------------------------------------------------------
# 6. 单组 MR-PRESSO（含全部异常处理）
# -----------------------------------------------------------------------------
run_one_pair <- function(exposure_name, exposure_id, outcome_name, outcome_id,
                         exp_dat, out_dat) {
  row <- data.frame(
    Exposure = exposure_name, Exposure_ID = exposure_id,
    Outcome = outcome_name, Outcome_ID = outcome_id,
    MR_PRESSO_global_P = NA_character_,
    Outlier_SNPs_identified = NA_character_,
    IVW_beta_before = NA_real_, IVW_OR_95CI_before = NA_character_,
    IVW_beta_after  = NA_real_, IVW_OR_95CI_after  = NA_character_,
    stringsAsFactors = FALSE
  )

  # ---- 6.1 数据调和与 QC ----
  har <- tryCatch(harmonise_qc(exp_dat, out_dat),
                  error = function(e) e)
  if (inherits(har, "error")) {
    row$MR_PRESSO_global_P <- "Insufficient SNPs"
    row$Outlier_SNPs_identified <- paste0("Insufficient SNPs (harmonise error: ",
                                          conditionMessage(har), ")")
    log_msg("  调和失败：", conditionMessage(har))
    return(row)
  }
  if (nrow(har) < 4) {
    # MR-PRESSO 至少需要 4 个独立工具变量
    row$MR_PRESSO_global_P <- "Insufficient SNPs"
    row$Outlier_SNPs_identified <- paste0("Insufficient SNPs (n = ", nrow(har), ")")
    log_msg("  SNP 不足（n = ", nrow(har), "），跳过 MR-PRESSO")
    return(row)
  }

  # ---- 6.2 剔除离群前的 IVW ----
  before <- tryCatch(ivw_summary(har), error = function(e) list(beta = NA_real_, ci = NA_character_))
  row$IVW_beta_before <- before$beta
  row$IVW_OR_95CI_before <- before$ci

  # ---- 6.3 MR-PRESSO 主分析（固定种子保证可复现）----
  # 本机 MRPRESSO 1.0 的函数签名：
  #   mr_presso(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data,
  #             OUTLIERtest, DISTORTIONtest, SignifThreshold,
  #             NbDistribution, seed)
  # 把行名设为 SNP rsID，便于把 Outlier Test 的行号映射回具体 SNP
  rownames(har) <- har$SNP
  set.seed(set_seed)
  presso <- tryCatch(
    mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
              SdOutcome = "se.outcome", SdExposure = "se.exposure",
              data = har,
              OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
              SignifThreshold = signif_thr,
              NbDistribution = nb_dist, seed = set_seed),
    error = function(e) e
  )
  if (inherits(presso, "error")) {
    row$MR_PRESSO_global_P <- "Insufficient SNPs"
    row$Outlier_SNPs_identified <- paste0("MR-PRESSO error: ", conditionMessage(presso))
    log_msg("  MR-PRESSO 报错：", conditionMessage(presso))
    return(row)
  }

  # 全局检验 P 值（模拟 P=0 时包内会返回 "<0.001" 形式的字符值，原样保留）
  g_p <- presso$`MR-PRESSO results`$`Global Test`$Pvalue
  row$MR_PRESSO_global_P <- as.character(g_p)

  # 离群 SNP：Outlier Test 的行名即 rsID；其 P 值可能是数值，也可能是
  # "<0.005" 这类字符（字符型一律视为显著）
  ot <- presso$`MR-PRESSO results`$`Outlier Test`
  outliers <- character(0)
  if (!is.null(ot) && is.data.frame(ot)) {
    p_num <- suppressWarnings(as.numeric(ot$Pvalue))
    hit_num <- which(!is.na(p_num) & p_num <= signif_thr)
    hit_chr <- which(is.na(p_num) & grepl("^<", as.character(ot$Pvalue)))
    outliers <- rownames(ot)[union(hit_num, hit_chr)]
  }
  if (length(outliers) == 0) {
    row$Outlier_SNPs_identified <- "No outlier detected"
  } else {
    row$Outlier_SNPs_identified <- paste(outliers, collapse = "; ")
  }

  # ---- 6.4 剔除离群后的 IVW ----
  if (length(outliers) == 0) {
    row$IVW_beta_after <- before$beta
    row$IVW_OR_95CI_after <- before$ci
  } else {
    har2 <- har[!har$SNP %in% outliers, ]
    if (nrow(har2) < 2) {
      row$IVW_OR_95CI_after <- "Insufficient SNPs after outlier removal"
      log_msg("  剔除离群后剩余 SNP 过少（n = ", nrow(har2), "）")
    } else {
      after <- tryCatch(ivw_summary(har2),
                        error = function(e) list(beta = NA_real_, ci = "IVW failed after removal"))
      row$IVW_beta_after <- after$beta
      row$IVW_OR_95CI_after <- after$ci
    }
  }
  log_msg("  全局 P = ", row$MR_PRESSO_global_P, " | 离群：",
          row$Outlier_SNPs_identified)
  row
}

# -----------------------------------------------------------------------------
# 7. 主循环：2 表型 × 32 癌症 × 2 方向 = 128 组
# -----------------------------------------------------------------------------
results <- list()
counter <- 0
t0 <- Sys.time()

for (ti in seq_len(nrow(trait_map))) {
  tname <- trait_map$trait_name[ti]
  tid   <- trait_map$trait_id[ti]

  # 本地 FinnGen 暴露工具变量（正向用，提前读一次）
  local_exp <- tryCatch(get_local_exposure(tid, tname), error = function(e) e)

  for (ci in seq_len(nrow(cancer_map))) {
    cname <- cancer_map$cancer_subgroup[ci]
    cid   <- cancer_map$cancer_id[ci]

    # ---------------- 正向：生殖表型 -> 癌症 ----------------
    counter <- counter + 1
    log_msg(sprintf("[%d/128] FORWARD %s (%s) -> %s (%s)",
                    counter, tname, tid, cname, cid))
    fwd <- tryCatch({
      if (inherits(local_exp, "error")) stop(conditionMessage(local_exp))
      if (nrow(local_exp) == 0) stop("本地暴露工具变量为 0")
      c_out <- get_cancer_outcome(cid, local_exp$SNP)
      if (is.null(c_out) || nrow(c_out) == 0) stop("癌症结局提取为空")
      run_one_pair(tname, tid, cname, cid, local_exp, c_out)
    }, error = function(e) {
      log_msg("  正向整体失败：", conditionMessage(e))
      data.frame(Exposure = tname, Exposure_ID = tid, Outcome = cname, Outcome_ID = cid,
                 MR_PRESSO_global_P = "Insufficient SNPs",
                 Outlier_SNPs_identified = paste0("Pipeline error: ", conditionMessage(e)),
                 IVW_beta_before = NA_real_, IVW_OR_95CI_before = NA_character_,
                 IVW_beta_after = NA_real_, IVW_OR_95CI_after = NA_character_,
                 stringsAsFactors = FALSE)
    })
    results[[length(results) + 1]] <- fwd

    # ---------------- 反向：癌症 -> 生殖表型 ----------------
    counter <- counter + 1
    log_msg(sprintf("[%d/128] REVERSE %s (%s) -> %s (%s)",
                    counter, cname, cid, tname, tid))
    rev <- tryCatch({
      c_exp <- get_cancer_exposure(cid)
      if (is.null(c_exp) || nrow(c_exp) == 0) stop("癌症暴露工具变量为 0")
      l_out <- get_local_outcome(tid, tname, c_exp$SNP)
      if (nrow(l_out) == 0) stop("本地结局提取为空")
      run_one_pair(cname, cid, tname, tid, c_exp, l_out)
    }, error = function(e) {
      log_msg("  反向整体失败：", conditionMessage(e))
      data.frame(Exposure = cname, Exposure_ID = cid, Outcome = tname, Outcome_ID = tid,
                 MR_PRESSO_global_P = "Insufficient SNPs",
                 Outlier_SNPs_identified = paste0("Pipeline error: ", conditionMessage(e)),
                 IVW_beta_before = NA_real_, IVW_OR_95CI_before = NA_character_,
                 IVW_beta_after = NA_real_, IVW_OR_95CI_after = NA_character_,
                 stringsAsFactors = FALSE)
    })
    results[[length(results) + 1]] <- rev

    # 每完成一组癌症的双向分析，落盘一次进度（防止中断丢结果）
    progress <- do.call(rbind, results)
    fwrite(progress, file.path(out_dir, "MR_PRESSO_progress.csv"))
  }

  # 释放当前表型的本地 FinnGen 全量缓存（每个文件上千万行，避免两个表型同时驻留）
  if (exists(tid, envir = cache_env)) {
    rm(list = tid, envir = cache_env)
    gc(verbose = FALSE)
    log_msg("已释放 ", tid, " 的本地缓存")
  }
}

# -----------------------------------------------------------------------------
# 8. 汇总写出（列名严格按要求）
# -----------------------------------------------------------------------------
final <- do.call(rbind, results)
final <- final[, c("Exposure", "Exposure_ID", "Outcome", "Outcome_ID",
                   "MR_PRESSO_global_P", "Outlier_SNPs_identified",
                   "IVW_beta_before", "IVW_OR_95CI_before",
                   "IVW_beta_after", "IVW_OR_95CI_after")]
stopifnot(nrow(final) == 128)

fwrite(final, final_csv)                       # 任务指定路径
fwrite(final, file.path(script_dir, "MR_PRESSO_full_results.csv"))  # 同时在 MR-PRESSO 目录留档
log_msg("全部完成，共 ", nrow(final), " 行；耗时 ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " 分钟")
log_msg("结果文件：", final_csv)
