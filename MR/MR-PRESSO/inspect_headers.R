# 核对 FinnGen 两个文件的表头，并测试 OpenGWAS 连通性
library(data.table)

mi <- "D:/研究生科研小组/小论文/GWAS数据库/finn数据/finngen_R12_N14_MALEINFERT.gz"
as <- "D:/研究生科研小组/小论文/GWAS数据库/finn数据/finngen_R5_R18_ABNORMAL_SPERMATOZ.gz"

h1 <- names(fread(mi, nrows = 3, showProgress = FALSE))
h2 <- names(fread(as, nrows = 3, showProgress = FALSE))
cat("MI header:\n"); print(h1)
cat("AS header:\n"); print(h2)

# 注入令牌（沿用既有脚本中的 token）
# OpenGWAS 令牌：从环境变量 OPENGWAS_JWT 读取，严禁在源码中硬编码 (见 .gitignore)
my_token <- Sys.getenv("OPENGWAS_JWT")
if (nchar(my_token) == 0) {
  stop("环境变量 OPENGWAS_JWT 未设置，请先 Sys.setenv(OPENGWAS_JWT=\"<your token>\") 或写入 .Renviron")
}
# ieugwasr 新版会自动读取 OPENGWAS_JWT 环境变量；此处仅确保其已设置
Sys.setenv(OPENGWAS_JWT = my_token)
library(ieugwasr)
library(TwoSampleMR)
# 测试 API：取一个数据集的元信息
info <- tryCatch(gwasinfo("finn-b-C3_PROSTATE"), error = function(e) e)
if (inherits(info, "error")) {
  cat("gwasinfo ERROR:", conditionMessage(info), "\n")
} else {
  cat("gwasinfo OK, sample size:", info$sample_size, "\n")
}
# 小样本提取测试
x <- tryCatch(
  extract_outcome_data(snps = c("rs902773", "rs103294"), outcomes = "finn-b-C3_PROSTATE", proxies = FALSE),
  error = function(e) e
)
if (inherits(x, "error")) {
  cat("extract ERROR:", conditionMessage(x), "\n")
} else {
  cat("extract OK, rows:", nrow(x), "\n")
}
