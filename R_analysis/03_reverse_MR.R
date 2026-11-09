# =============================================================================
# 批量反向 Mendelian Randomization（反向 MR）
# 方向：各类癌症（Exposure） -> 男性不育 MI（Outcome, finngen_R12_N14_MALEINFERT）
# 数据来源：全程 ieugwasr + TwoSampleMR 在线 OpenGWAS API
#           结局 MI 若在线不可用，自动回退到本地 FinnGen R12 文件（仅结局回退）
#
# 严格匹配手稿的核心参数（请勿改动）：
#   工具变量阈值 p1 = 1e-5
#   聚类 clump_r2 = 0.001, clump_kb = 10000
#   harmonise_data(action = 2)
#   去除 palindromic / ambiguous SNP
#   去除 MAF < 0.01 的 SNP
#   方法：IVW、MR-Egger、Weighted median、Simple mode、Weighted mode
#   异质性 mr_heterogeneity()；多效性 mr_pleiotropy_test()
#
# 本脚本可直接整体复制到 RStudio 运行，无需手动 auth()，无需修改任何参数。
# =============================================================================


# -----------------------------------------------------------------------------
# 0. 依赖包检查与加载（缺失则自动安装）
# -----------------------------------------------------------------------------
options(stringsAsFactors = FALSE)

# 0.1 CRAN 依赖
.cran_pkgs <- c("data.table", "remotes")
for (pkg in .cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("安装 CRAN 包：", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# 0.2 GitHub 上的 TwoSampleMR 与 ieugwasr（OpenGWAS 版本）
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  message("从 GitHub 安装 TwoSampleMR ...")
  remotes::install_github("MRCIEU/TwoSampleMR", upgrade = "never")
}
if (!requireNamespace("ieugwasr", quietly = TRUE)) {
  message("从 GitHub 安装 ieugwasr ...")
  remotes::install_github("MRCIEU/ieugwasr", upgrade = "never")
}

library(TwoSampleMR)
library(ieugwasr)
library(data.table)


# -----------------------------------------------------------------------------
# 1. 固定参数（与手稿完全一致）
# -----------------------------------------------------------------------------
p1_thr   <- 1e-5      # 工具变量显著性阈值（1e-5，不是 5e-8）
clump_r2 <- 0.001     # LD 聚类 r2
clump_kb <- 10000     # LD 聚类窗口 kb
maf_thr  <- 0.01      # 次要等位基因频率下限
use_proxies <- FALSE  # 结局匹配不使用 LD 代理 SNP（严格 QC，可复现）
                        # 若手稿确需代理，将此处改为 TRUE 即可

# 固定结局：男性不育
outcome_id   <- "finngen_R12_N14_MALEINFERT"
# 本地 FinnGen 备份（仅当在线 API 不可用时自动回退）
local_fg_path <- "D:/研究生科研小组/小论文/GWAS数据库/finn数据/finngen_R12_N14_MALEINFERT.gz"

# OpenGWAS API token — NEVER hard-code a token in version-controlled code.
# Read the token from the environment variable OPENGWAS_JWT (fallback OPENGWAS_TOKEN).
# Set it once with Sys.setenv(OPENGWAS_JWT = "your-token"), or define it permanently
# in the Windows environment variables.
my_token <- Sys.getenv("OPENGWAS_JWT", unset = Sys.getenv("OPENGWAS_TOKEN"))
if (!nzchar(my_token)) {
  stop("OpenGWAS token not found. Set the OPENGWAS_JWT environment variable before running.")
}


# -----------------------------------------------------------------------------
# 2. OpenGWAS non-interactive authentication (no browser login)
# -----------------------------------------------------------------------------
# Current ieugwasr (1.2.x) removed ieugwasr::auth(); online requests read the token
# automatically from the OPENGWAS_JWT environment variable.
Sys.setenv(OPENGWAS_JWT = my_token)
message("OpenGWAS token loaded from the OPENGWAS_JWT environment variable.")


# -----------------------------------------------------------------------------
# 3. 癌症暴露清单（cancer_subgroup 与 exposure_id 一一对应，共 32 个数据集）
# -----------------------------------------------------------------------------
cancer_map <- data.frame(
  cancer_subgroup = c(
    # Colon cancer (3)
    rep("Colon cancer", 3),
    # Colorectal cancer (4)
    rep("Colorectal cancer", 4),
    # Liver cancer (3)
    rep("Liver cancer", 3),
    # Oesophagus cancer (1)
    rep("Oesophagus cancer", 1),
    # Gastric cancer (2)
    rep("Gastric cancer", 2),
    # Bladder cancer (3)
    rep("Bladder cancer", 3),
    # Prostate cancer (2)
    rep("Prostate cancer", 2),
    # Testicular cancer (1)
    rep("Testicular cancer", 1),
    # Leukaemia (4)
    rep("Leukaemia", 4),
    # Non-Hodgkin lymphoma (1)
    rep("Non-Hodgkin lymphoma", 1),
    # Hodgkin lymphoma (1)
    rep("Hodgkin lymphoma", 1),
    # Melanoma (skin cancer) (5)
    rep("Melanoma (skin cancer)", 5),
    # Thyroid cancer (endocrine) (2)
    rep("Thyroid cancer (endocrine)", 2)
  ),
  exposure_id = c(
    # Colon cancer
    "ukb-b-20145", "ukb-d-C3_COLON", "finn-b-C3_COLON",
    # Colorectal cancer
    "ebi-a-GCST90013862", "ebi-a-GCST90013866", "finn-b-C3_COLORECTAL", "ieu-b-4965",
    # Liver cancer
    "ebi-a-GCST90018638", "ebi-a-GCST90018858", "ieu-b-4953",
    # Oesophagus cancer
    "finn-b-C3_OESOPHAGUS",
    # Gastric cancer
    "ebi-a-GCST90018629", "ebi-a-GCST90018849",
    # Bladder cancer
    "ukb-b-8193", "ieu-b-4874", "finn-b-C3_BLADDER",
    # Prostate cancer
    "finn-b-C3_PROSTATE", "ieu-b-85",
    # Testicular cancer
    "finn-b-C3_TESTIS",
    # Leukaemia
    "ieu-b-4914", "finn-b-C3_AML", "ieu-b-4956", "ieu-b-4958",
    # Non-Hodgkin lymphoma
    "finn-b-CD2_NONHODGKIN_NAS",
    # Hodgkin lymphoma
    "finn-b-CD2_HODGKIN_LYMPHOMA",
    # Melanoma (skin cancer)
    "finn-b-C3_MELANOMA_SKIN", "ukb-d-C3_MELANOMA_SKIN", "ieu-b-4969",
    "ukb-a-58", "ukb-b-12915",
    # Thyroid cancer (endocrine)
    "ebi-a-GCST90018929", "ieu-a-1082"
  )
)
stopifnot(nrow(cancer_map) == 32)
message("待分析癌症暴露数据集总数：", nrow(cancer_map))


# -----------------------------------------------------------------------------
# 4. 输出目录：自动创建 output 文件夹
#    最终结果存放于 MI 目录下的 output 子文件夹
# -----------------------------------------------------------------------------
base_dir <- "D:/研究生科研小组/小论文/meta+MR/MR/MI"   # 如需当前工作目录，改为 "."
out_dir  <- file.path(base_dir, "output")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
message("结果输出目录：", out_dir)

# 文件名安全清洗（替换非法字符）
sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

# 安全取首个元素（防止长度为 0 时下标越界）
safe_first <- function(x) {
  if (length(x) >= 1) x[1] else NA
}


# -----------------------------------------------------------------------------
# 5. 结局 MI 数据获取：在线优先，失败自动回退本地 FinnGen
# -----------------------------------------------------------------------------

# 5.1 探测在线结局是否可用
message("探测在线 OpenGWAS 结局 ", outcome_id, " ...")
outcome_online <- tryCatch(
  {
    info <- ieugwasr::gwasinfo(outcome_id)
    !is.null(info) && length(info) > 0
  },
  error = function(e) FALSE
)

if (outcome_online) {
  message("在线结局可用，将通过 OpenGWAS API 提取结局数据。")
} else {
  message("在线结局不可用（网络/权限/ID 问题），将自动回退到本地文件：")
  message(local_fg_path)
}

# 5.2 本地 FinnGen 全量读取（只读一次并缓存）
# 用独立环境持有缓存，避免函数内普通赋值只生成局部变量、导致每次调用都
# 重新 fread 的作用域问题。全量 2000 多万行只读一次，后续 ID 直接复用。
.fg_env <- new.env(parent = emptyenv())
.fg_env$dat <- NULL
load_local_fg <- function() {
  if (!is.null(.fg_env$dat)) {
    message("复用已缓存的本地 FinnGen 结局数据（无需重新读取文件）。")
    return(.fg_env$dat)
  }
  message("首次本地加载 FinnGen（文件较大，约 768MB 压缩包，请耐心等待）...")

  # 仅读取需要的列，降低内存占用
  fg <- data.table::fread(
    input       = local_fg_path,
    select      = c("rsids", "ref", "alt", "pval", "beta", "sebeta", "af_alt"),
    data.table  = FALSE,
    showProgress = FALSE
  )

  # rsids 清洗：多个 rsID 以逗号/分号分隔时取第一个；剔除空值与非 rsID 占位符
  # 注意：空 rsids 经 strsplit 会得到 character(0)，此时必须返回字符型 NA
  # （NA_character_），裸 NA 是 logical 类型，会让 vapply 报类型冲突错误。
  first_rs <- trimws(vapply(
    strsplit(fg$rsids, "[,;]"),
    function(z) if (length(z) >= 1) z[1] else NA_character_,
    character(1)
  ))
  fg$SNP <- first_rs
  fg <- fg[!is.na(fg$SNP) & nzchar(fg$SNP) & grepl("^rs[0-9]+$", fg$SNP), ]
  fg <- fg[!duplicated(fg$SNP), ]

  # 格式化为标准 outcome 数据框
  # 新版 TwoSampleMR 0.7.13 的 format_data 已无 maf_col 参数，仅保留 eaf_col；
  # af_alt 为 effect allele(alt) 的等位基因频率，直接作为 eaf 读入，
  # 后续 MAF 统一由 eaf 经 pmin(eaf, 1-eaf) 推导。
  out_local <- TwoSampleMR::format_data(
    fg,
    type              = "outcome",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "sebeta",
    effect_allele_col = "alt",   # alt 为效应等位
    other_allele_col  = "ref",   # ref 为其他等位
    pval_col          = "pval",
    eaf_col           = "af_alt" # alt(effect allele) 的频率
  )
  out_local$id.outcome <- outcome_id
  out_local$outcome    <- "Male infertility (FinnGen R12)"

  .fg_env$dat <- out_local
  message("本地 FinnGen 加载完成，可用 SNP 数：", nrow(out_local))
  return(out_local)
}

# 5.3 按暴露 SNP 集合获取结局数据（在线 / 本地统一入口）
get_outcome_dat <- function(snps) {
  if (outcome_online) {
    od <- tryCatch(
      TwoSampleMR::extract_outcome_data(
        snps     = snps,
        outcomes = outcome_id,
        proxies  = use_proxies
      ),
      error = function(e) NULL
    )
    if (!is.null(od) && nrow(od) > 0) {
      return(od)
    }
    message("在线结局本次未返回有效数据，自动切换本地 FinnGen ...")
  }

  full <- load_local_fg()
  od <- full[full$SNP %in% snps, ]
  if (nrow(od) == 0) {
    stop("本地 FinnGen 结局中也未匹配到任何暴露 SNP")
  }
  return(od)
}


# -----------------------------------------------------------------------------
# 6. 五种 MR 方法（同时运行）
# -----------------------------------------------------------------------------
# 新版 mr() 的 method_list 必须是“方法名字符串”（mr_method_list()$obj），
# 不能再传函数对象。以下 5 个字符串与手稿要求的 5 种方法一一对应。
method_list <- c(
  "mr_ivw",               # Inverse variance weighted
  "mr_egger_regression",  # MR-Egger
  "mr_weighted_median",   # Weighted median
  "mr_simple_mode",       # Simple mode
  "mr_weighted_mode"      # Weighted mode
)


# -----------------------------------------------------------------------------
# 7. 汇总容器
# -----------------------------------------------------------------------------
summary_rows       <- list()  # 全部方法结果
hetero_pleio_rows  <- list()  # 异质性 + 多效性
ivw_rows           <- list()  # 仅 IVW


# -----------------------------------------------------------------------------
# 8. 主循环：每个癌症暴露独立在线提取自己的工具 SNP（严禁复用其他 SNP）
# -----------------------------------------------------------------------------
for (i in seq_len(nrow(cancer_map))) {

  cid  <- cancer_map$exposure_id[i]
  csub <- cancer_map$cancer_subgroup[i]
  safe <- sanitize_name(cid)

  message("\n============================================================")
  message(sprintf("[%d/%d] %s | exposure ID = %s", i, nrow(cancer_map), csub, cid))

  tryCatch(
    {
      # 8.1 独立提取“该癌症自己”的工具变量
      # 注意：extract_instruments 的参数名为 r2 / kb（即 clump_r2 / clump_kb）
      exp_dat <- TwoSampleMR::extract_instruments(
        outcomes = cid,
        p1       = p1_thr,      # 1e-5
        clump    = TRUE,
        r2       = clump_r2,    # 0.001
        kb       = clump_kb     # 10000
      )

      if (is.null(exp_dat) || nrow(exp_dat) == 0) {
        message(sprintf("WARNING: 暴露 %s 在线提取到 0 个工具 SNP，跳过该 ID。", cid))
        next
      }
      message("提取到工具 SNP 数（clump 后）：", nrow(exp_dat))

      # 8.2 保存该癌症的工具 SNP 原始表
      data.table::fwrite(
        exp_dat,
        file.path(out_dir, paste0("exposure_instr_", safe, ".csv"))
      )

      # 8.3 获取固定结局 MI 的对应数据
      out_dat <- get_outcome_dat(exp_dat$SNP)

      # 8.4 协调数据（action = 2：正链推断，无法确定则剔除）
      har_dat <- TwoSampleMR::harmonise_data(
        exposure_dat = exp_dat,
        outcome_dat  = out_dat,
        action       = 2
      )

      # 仅保留 harmonise 通过的 SNP
      har_dat <- har_dat[har_dat$mr_keep == TRUE, ]

      # 8.5 去除 palindromic / ambiguous SNP
      if ("palindromic" %in% names(har_dat)) {
        n_before <- nrow(har_dat)
        har_dat <- har_dat[!har_dat$palindromic, ]
        message("去除 palindromic SNP：", n_before - nrow(har_dat), " 个")
      }

      # 8.6 剔除 MAF < 0.01（暴露与结局两侧 MAF 的最小值需 >= 0.01）
      # 新版只有 eaf 列、没有 maf 列，MAF 由 eaf 自行推导：maf = pmin(eaf, 1-eaf)
      # 先确保两侧 eaf 列存在（缺失则补 NA），避免下标异常
      if (!"eaf.exposure" %in% names(har_dat)) har_dat$eaf.exposure <- NA_real_
      if (!"eaf.outcome" %in% names(har_dat))  har_dat$eaf.outcome  <- NA_real_
      maf_exposure <- pmin(har_dat$eaf.exposure, 1 - har_dat$eaf.exposure)
      maf_outcome  <- pmin(har_dat$eaf.outcome,  1 - har_dat$eaf.outcome)
      maf_min <- apply(
        cbind(maf_exposure, maf_outcome),
        1,
        function(z) {
          z <- z[!is.na(z)]
          if (length(z) == 0) NA_real_ else min(z)
        }
      )
      n_before <- nrow(har_dat)
      har_dat <- har_dat[!is.na(maf_min) & maf_min >= maf_thr, ]
      message("剔除 MAF < 0.01 的 SNP：", n_before - nrow(har_dat), " 个")

      if (nrow(har_dat) == 0) {
        message(sprintf("WARNING: 暴露 %s 经 palindromic/MAF 过滤后剩余 0 个 SNP，跳过该 ID。", cid))
        next
      }
      message("最终进入 MR 的 SNP 数：", nrow(har_dat))

      # 8.7 保存 harmonise 后完整数据集（含 MAF）
      # 新版 harmonise 结果本身只有 eaf，这里显式补上两侧 MAF 列，
      # 方便后续直接核对 MAF<0.01 的剔除情况。
      har_dat$maf.exposure <- pmin(har_dat$eaf.exposure, 1 - har_dat$eaf.exposure)
      har_dat$maf.outcome  <- pmin(har_dat$eaf.outcome,  1 - har_dat$eaf.outcome)
      data.table::fwrite(
        har_dat,
        file.path(out_dir, paste0("harmonised_", safe, ".csv"))
      )

      # 8.8 运行五种 MR 方法，并计算 OR 及 95% CI
      mr_res <- TwoSampleMR::mr(har_dat, method_list = method_list)
      mr_res_or <- TwoSampleMR::generate_odds_ratios(mr_res)

      data.table::fwrite(
        mr_res_or,
        file.path(out_dir, paste0("mr_raw_", safe, ".csv"))
      )

      # 8.9 异质性检验（Q、Q_df、Q_pval）
      het <- TwoSampleMR::mr_heterogeneity(har_dat)

      # 8.10 多效性检验（MR-Egger intercept 及 P 值）
      pleio <- TwoSampleMR::mr_pleiotropy_test(har_dat)

      # 8.11 汇总：全部方法结果（用于 Supplementary Table S8）
      sdf <- data.frame(
        cancer_subgroup = csub,
        exposure_id     = cid,
        outcome_id      = outcome_id,
        nsnp            = mr_res_or$nsnp,
        method          = mr_res_or$method,
        beta            = mr_res_or$b,
        se              = mr_res_or$se,
        pval            = mr_res_or$pval,
        or              = mr_res_or$or,
        or_lci95        = mr_res_or$or_lci95,
        or_uci95        = mr_res_or$or_uci95
      )
      summary_rows[[length(summary_rows) + 1]] <- sdf

      # 8.12 汇总：异质性 + 多效性宽表（每个 ID 一行）
      het_ivw   <- het[het$method == "Inverse variance weighted", ]
      het_egger <- het[het$method == "MR Egger", ]

      hpdf <- data.frame(
        cancer_subgroup = csub,
        exposure_id     = cid,
        outcome_id      = outcome_id,
        nsnp            = nrow(har_dat),
        Q               = safe_first(het_ivw$Q),
        Q_df            = safe_first(het_ivw$Q_df),
        Q_pval          = safe_first(het_ivw$Q_pval),
        Q_egger         = safe_first(het_egger$Q),
        Q_df_egger      = safe_first(het_egger$Q_df),
        Q_pval_egger    = safe_first(het_egger$Q_pval),
        egger_intercept = safe_first(pleio$egger_intercept),
        intercept_se    = safe_first(pleio$se),
        intercept_pval  = safe_first(pleio$pval)
      )
      hetero_pleio_rows[[length(hetero_pleio_rows) + 1]] <- hpdf

      # 8.13 汇总：仅 IVW 结果（用于批量计算 MDOR，Supplementary Table S11）
      ivw_df <- mr_res_or[grepl("inverse variance", mr_res_or$method, ignore.case = TRUE), ]
      if (nrow(ivw_df) >= 1) {
        ivw_rows[[length(ivw_rows) + 1]] <- data.frame(
          cancer_subgroup = csub,
          exposure_id     = cid,
          nsnp            = safe_first(ivw_df$nsnp),
          ivw_beta        = safe_first(ivw_df$b),
          ivw_se          = safe_first(ivw_df$se),
          ivw_pval        = safe_first(ivw_df$pval)
        )
      } else {
        message("注意：", cid, " 未产生 IVW 结果，ivw_only 中不包含该 ID。")
      }
    },
    error = function(e) {
      # 网络错误 / API 报错 / 提取失败：打印信息并跳过，不中断整个循环
      message(sprintf("WARNING: ID %s 处理失败，已跳过。错误信息：%s",
                      cid, conditionMessage(e)))
    }
  )
}


# -----------------------------------------------------------------------------
# 9. 循环结束：写出三个汇总 CSV
# -----------------------------------------------------------------------------
save_bind <- function(lst, filename) {
  path <- file.path(out_dir, filename)
  if (length(lst) == 0) {
    message("没有成功结果，输出空表：", filename)
    data.table::fwrite(data.frame(), path)
  } else {
    data.table::fwrite(do.call(rbind, lst), path)
  }
  message("已写出：", path)
}

# 9.1 全部方法汇总
save_bind(summary_rows,      "reverse_MR_summary_all.csv")
# 9.2 异质性 + 多效性汇总
save_bind(hetero_pleio_rows, "reverse_MR_hetero_pleio.csv")
# 9.3 仅 IVW 汇总
save_bind(ivw_rows,          "reverse_MR_ivw_only.csv")


# -----------------------------------------------------------------------------
# 10. 控制台结束提示
# -----------------------------------------------------------------------------
message("\n============================================================")
message("全部反向 MR 分析完成！")
message("所有结果已保存至：", out_dir, "（即 MI 目录下的 ./output 文件夹）")
message("请使用 reverse_MR_ivw_only.csv 批量计算 MDOR（更新 Supplementary Table S11）。")
message("请核对筛选参数：p1 = 1e-5、clump_r2 = 0.001、clump_kb = 10000、")
message("已剔除 palindromic/ambiguous SNP 及 MAF < 0.01 的 SNP。")
message("============================================================")
