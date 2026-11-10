# =============================================================================
# 全量双向 MR-PRESSO —— 断点续跑脚本
# -----------------------------------------------------------------------------
# 背景：
#   主脚本 run_bidirectional_MR_PRESSO.R 在运行到第 100 组时进程被外部终止
#   （R 层面无报错，日志无错误信息），output/MR_PRESSO_progress.csv 已保存
#   前 98 组结果。本脚本：
#     1) 读取已有进度；
#     2) 按 128 组的标准顺序生成完整任务计划；
#     3) 以 “Exposure_ID | Outcome_ID” 为唯一键，自动跳过已完成的组，
#        只运行缺失的组（本次预期为第 99–128 组，共 30 组）；
#     4) 合并新旧结果，严格按 128 组顺序写出最终 CSV（两个路径各一份）。
#
# 可复现性 / QC 参数与主脚本完全一致：
#   set.seed(1234)；p1 = 1e-5；clump r2 = 0.001, kb = 10000；
#   harmonise action = 2；剔除 palindromic；剔除 MAF < 0.01；
#   MR-PRESSO NbDistribution = 1000，SignifThreshold = 0.05。
#
# 本脚本可安全重复执行：已完成的组不会重跑，写入仅限输出文件。
# =============================================================================

# -----------------------------------------------------------------------------
# 0. 基础环境（与主脚本一致）
# -----------------------------------------------------------------------------
options(stringsAsFactors = FALSE)

library(data.table)
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(MRPRESSO)

# ---- 0.1 OpenGWAS 令牌（从环境变量注入，ieugwasr 自动携带）----
# OpenGWAS 令牌：从环境变量 OPENGWAS_JWT 读取，严禁在源码中硬编码 (见 .gitignore)
my_token <- Sys.getenv("OPENGWAS_JWT")
if (nchar(my_token) == 0) {
  stop("环境变量 OPENGWAS_JWT 未设置，请先 Sys.setenv(OPENGWAS_JWT=\"<your token>\") 或写入 .Renviron")
}
# ieugwasr 新版会自动读取 OPENGWAS_JWT 环境变量；此处仅确保其已设置
Sys.setenv(OPENGWAS_JWT = my_token)
# ---- 0.2 固定分析参数（与主脚本一致，请勿改动）----
p1_thr       <- 1e-5
clump_r2     <- 0.001
clump_kb     <- 10000
maf_thr      <- 0.01
use_proxies  <- FALSE
nb_dist      <- 1000
signif_thr   <- 0.05
set_seed     <- 1234

# ---- 0.3 路径与输出目录 ----
script_dir <- "D:/研究生科研小组/小论文/meta+MR/MR/MR-PRESSO"
out_dir    <- file.path(script_dir, "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
final_csv  <- "D:/研究生科研小组/小论文/MR_PRESSO_full_results.csv"
log_file   <- file.path(out_dir, "run_log.txt")
progress_file <- file.path(out_dir, "MR_PRESSO_progress.csv")

log_msg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

# -----------------------------------------------------------------------------
# 1. 32 个癌症数据集清单（与主脚本完全一致）
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
# 2. 两个本地 FinnGen 生殖表型（与主脚本一致）
# -----------------------------------------------------------------------------
trait_map <- data.frame(
  trait_name = c("Male infertility", "Abnormal spermatozoa"),
  trait_id   = c("finngen_R12_N14_MALEINFERT", "finngen_R5_R18_ABNORMAL_SPERMATOZ"),
  file_path  = c(
    "D:/研究生科研小组/小论文/GWAS数据库/finn数据/finngen_R12_N14_MALEINFERT.gz",
    "D:/研究生科研小组/小论文/GWAS数据库/finn数据/finngen_R5_R18_ABNORMAL_SPERMATOZ.gz"
  )
)

cache_env <- new.env(parent = emptyenv())

read_finngen_std <- function(trait_id) {
  if (!is.null(cache_env[[trait_id]])) return(cache_env[[trait_id]])

  path <- trait_map$file_path[trait_map$trait_id == trait_id]
  log_msg("首次读取 FinnGen 文件：", path)

  header <- names(fread(path, nrows = 0, showProgress = FALSE))
  eaf_col <- if ("af_alt" %in% header) "af_alt" else if ("maf" %in% header) "maf" else
    stop("未找到频率列 af_alt/maf：", path)

  fg <- fread(path,
              select = c("rsids", "ref", "alt", "pval", "beta", "sebeta", eaf_col),
              showProgress = FALSE)
  setnames(fg, eaf_col, "eaf")

  first_rs <- trimws(vapply(strsplit(fg$rsids, "[,;]"),
                            function(z) if (length(z) >= 1) z[1] else NA_character_,
                            character(1)))
  fg$SNP <- first_rs
  fg <- fg[!is.na(fg$SNP) & nzchar(fg$SNP) & grepl("^rs[0-9]+$", fg$SNP), ]
  fg <- fg[!duplicated(fg$SNP), ]

  std <- data.frame(
    SNP           = fg$SNP,
    beta          = fg$beta,
    se            = fg$sebeta,
    pval          = fg$pval,
    effect_allele = fg$alt,
    other_allele  = fg$ref,
    eaf           = fg$eaf,
    stringsAsFactors = FALSE
  )
  cache_env[[trait_id]] <- std
  log_msg(trait_id, " 标准化完成，可用 SNP 数：", nrow(std))
  std
}

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
# 3. 在线癌症数据（与主脚本一致）
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
# 4. 统一 QC（与主脚本一致）
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
# 5. IVW 计算工具（与主脚本一致）
# -----------------------------------------------------------------------------
ivw_summary <- function(har) {
  ivw <- mr(har, method_list = "mr_ivw")
  if (nrow(ivw) == 0) return(list(beta = NA_real_, ci = NA_character_))
  orr <- generate_odds_ratios(ivw)
  list(beta = orr$b[1],
       ci = sprintf("%.2f (%.2f, %.2f)", orr$or[1], orr$or_lci95[1], orr$or_uci95[1]))
}

# -----------------------------------------------------------------------------
# 6. 单组 MR-PRESSO（与主脚本一致）
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
    row$MR_PRESSO_global_P <- "Insufficient SNPs"
    row$Outlier_SNPs_identified <- paste0("Insufficient SNPs (n = ", nrow(har), ")")
    log_msg("  SNP 不足（n = ", nrow(har), "），跳过 MR-PRESSO")
    return(row)
  }

  before <- tryCatch(ivw_summary(har), error = function(e) list(beta = NA_real_, ci = NA_character_))
  row$IVW_beta_before <- before$beta
  row$IVW_OR_95CI_before <- before$ci

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

  g_p <- presso$`MR-PRESSO results`$`Global Test`$Pvalue
  row$MR_PRESSO_global_P <- as.character(g_p)

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
# 7. 构造 128 组的完整任务计划（顺序与主脚本严格一致）
# -----------------------------------------------------------------------------
plan <- list()
k <- 0
for (ti in seq_len(nrow(trait_map))) {
  tname <- trait_map$trait_name[ti]
  tid   <- trait_map$trait_id[ti]
  for (ci in seq_len(nrow(cancer_map))) {
    cname <- cancer_map$cancer_subgroup[ci]
    cid   <- cancer_map$cancer_id[ci]

    # 正向：生殖表型 -> 癌症
    k <- k + 1
    plan[[k]] <- data.frame(pair = k, direction = "FORWARD",
                            exp_name = tname, exp_id = tid,
                            out_name = cname, out_id = cid,
                            stringsAsFactors = FALSE)
    # 反向：癌症 -> 生殖表型
    k <- k + 1
    plan[[k]] <- data.frame(pair = k, direction = "REVERSE",
                            exp_name = cname, exp_id = cid,
                            out_name = tname, out_id = tid,
                            stringsAsFactors = FALSE)
  }
}
plan <- do.call(rbind, plan)
stopifnot(nrow(plan) == 128)

# 唯一键：暴露 ID + 结局 ID（方向已隐含在二者顺序中）
plan$key <- paste(plan$exp_id, plan$out_id, sep = "|")

# -----------------------------------------------------------------------------
# 8. 读取已有进度，映射到 128 个位置
# -----------------------------------------------------------------------------
# all_rows[[i]] 存放第 i 组的结果行；先全部初始化为 NULL
all_rows <- vector("list", 128)

if (file.exists(progress_file)) {
  done <- fread(progress_file, showProgress = FALSE)
  log_msg("读取已有进度：", nrow(done), " 行")
  done$key <- paste(done$Exposure_ID, done$Outcome_ID, sep = "|")
  for (i in seq_len(nrow(plan))) {
    hit <- done[done$key == plan$key[i], ]
    if (nrow(hit) >= 1) {
      hit <- hit[1, ]
      hit$key <- NULL
      all_rows[[i]] <- as.data.frame(hit, stringsAsFactors = FALSE)
    }
  }
}

done_keys <- plan$key[!vapply(all_rows, is.null, logical(1))]
log_msg("已完成 ", length(done_keys), " 组，待运行 ", 128 - length(done_keys), " 组")

# 组装当前完整表（用于落盘），统一列顺序
assemble_final <- function() {
  rows <- all_rows[!vapply(all_rows, is.null, logical(1))]
  df <- do.call(rbind, rows)
  df[, c("Exposure", "Exposure_ID", "Outcome", "Outcome_ID",
         "MR_PRESSO_global_P", "Outlier_SNPs_identified",
         "IVW_beta_before", "IVW_OR_95CI_before",
         "IVW_beta_after", "IVW_OR_95CI_after")]
}

# -----------------------------------------------------------------------------
# 9. 仅运行缺失组（按计划顺序，逻辑与主循环一致）
# -----------------------------------------------------------------------------
t0 <- Sys.time()

for (ti in seq_len(nrow(trait_map))) {
  tname <- trait_map$trait_name[ti]
  tid   <- trait_map$trait_id[ti]

  # 本地暴露延迟计算：只有该表型存在缺失的正向组时才需要
  local_exp <- NULL
  local_exp_err <- NULL

  for (ci in seq_len(nrow(cancer_map))) {
    cname <- cancer_map$cancer_subgroup[ci]
    cid   <- cancer_map$cancer_id[ci]

    # ---- 正向 pair 索引 ----
    fwd_idx <- (ti - 1) * 64 + (ci - 1) * 2 + 1
    rev_idx <- fwd_idx + 1

    # ---------------- 正向：生殖表型 -> 癌症 ----------------
    if (is.null(all_rows[[fwd_idx]])) {
      log_msg(sprintf("[%d/128] FORWARD %s (%s) -> %s (%s)",
                      fwd_idx, tname, tid, cname, cid))
      fwd <- tryCatch({
        if (is.null(local_exp)) {
          local_exp <- tryCatch(get_local_exposure(tid, tname),
                                error = function(e) e)
        }
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
      all_rows[[fwd_idx]] <- fwd
      fwrite(assemble_final(), progress_file)
    }

    # ---------------- 反向：癌症 -> 生殖表型 ----------------
    if (is.null(all_rows[[rev_idx]])) {
      log_msg(sprintf("[%d/128] REVERSE %s (%s) -> %s (%s)",
                      rev_idx, cname, cid, tname, tid))
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
      all_rows[[rev_idx]] <- rev
      fwrite(assemble_final(), progress_file)
    }
  }

  # 释放当前表型的本地缓存（与主脚本一致）
  if (exists(tid, envir = cache_env)) {
    rm(list = tid, envir = cache_env)
    gc(verbose = FALSE)
    log_msg("已释放 ", tid, " 的本地缓存")
  }
}

# -----------------------------------------------------------------------------
# 10. 汇总写出（列名 / 行数严格校验）
# -----------------------------------------------------------------------------
final <- assemble_final()
stopifnot(nrow(final) == 128)

fwrite(final, final_csv)
fwrite(final, file.path(script_dir, "MR_PRESSO_full_results.csv"))
log_msg("全部完成，共 ", nrow(final), " 行；本次续跑耗时 ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " 分钟")
log_msg("结果文件：", final_csv)
