# 冒烟测试：仅 1 个表型 × 1 个癌症 × 双向，验证函数签名与结果结构
options(stringsAsFactors = FALSE)
library(data.table); library(TwoSampleMR); library(ieugwasr); library(dplyr); library(MRPRESSO)
# OpenGWAS 令牌：从环境变量 OPENGWAS_JWT 读取，严禁在源码中硬编码 (见 .gitignore)
my_token <- Sys.getenv("OPENGWAS_JWT")
if (nchar(my_token) == 0) {
  stop("环境变量 OPENGWAS_JWT 未设置，请先 Sys.setenv(OPENGWAS_JWT=\"<your token>\") 或写入 .Renviron")
}
# ieugwasr 新版会自动读取 OPENGWAS_JWT 环境变量；此处仅确保其已设置
Sys.setenv(OPENGWAS_JWT = my_token)
# 直接复用主脚本中的函数与参数（source 到主循环之前会整体执行，故这里内联最小流程）
source_env <- new.env()
# 手动读取主脚本前 6 节（函数定义部分）：用 parse + 截断方式
lines <- readLines("D:/研究生科研小组/小论文/meta+MR/MR/MR-PRESSO/run_bidirectional_MR_PRESSO.R",
                   warn = FALSE, encoding = "UTF-8")
cut <- grep("^# 7\\. 主循环", lines)[1]
eval(parse(text = lines[1:(cut - 1)]), envir = source_env)

# 用小模拟次数快速验证
assign("nb_dist", 100, envir = source_env)

# 正向：AS -> PCa
le <- source_env$get_local_exposure("finngen_R5_R18_ABNORMAL_SPERMATOZ", "Abnormal spermatozoa")
cat("local exposure SNPs:", nrow(le), "\n")
co <- source_env$get_cancer_outcome("finn-b-C3_PROSTATE", le$SNP)
fwd <- source_env$run_one_pair("Abnormal spermatozoa", "finngen_R5_R18_ABNORMAL_SPERMATOZ",
                               "Prostate cancer", "finn-b-C3_PROSTATE", le, co)
cat("\n--- FORWARD result ---\n"); print(fwd)

# 反向：PCa -> AS
ce <- source_env$get_cancer_exposure("finn-b-C3_PROSTATE")
cat("cancer exposure SNPs:", nrow(ce), "\n")
lo <- source_env$get_local_outcome("finngen_R5_R18_ABNORMAL_SPERMATOZ", "Abnormal spermatozoa", ce$SNP)
rev <- source_env$run_one_pair("Prostate cancer", "finn-b-C3_PROSTATE",
                               "Abnormal spermatozoa", "finngen_R5_R18_ABNORMAL_SPERMATOZ",
                               ce, lo)
cat("\n--- REVERSE result ---\n"); print(rev)
cat("\nSMOKE TEST DONE\n")
