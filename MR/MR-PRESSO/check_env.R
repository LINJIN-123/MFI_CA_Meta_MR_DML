# 环境核对：R 版本、关键包、OpenGWAS 认证
cat("R version:", R.version.string, "\n")
pkgs <- c("TwoSampleMR", "MRPRESSO", "data.table", "ieugwasr", "dplyr")
for (p in pkgs) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) NA)
  cat(sprintf("pkg %-12s: %s\n", p, ifelse(is.na(v), "NOT INSTALLED", v)))
}
# OpenGWAS 令牌
cat("OPENGWAS_JWT env:", Sys.getenv("OPENGWAS_JWT", unset = "(empty)") != "", "\n")
tok <- tryCatch(ieugwasr::get_opengwas_jwt(), error = function(e) "")
if (!requireNamespace("ieugwasr", quietly = TRUE)) {
  cat("ieugwasr not installed, token check skipped\n")
} else {
  cat("token from ieugwasr:", nchar(tok) > 0, "| length:", nchar(tok), "\n")
}
# Renviron 位置
cat("R_USER:", Sys.getenv("R_USER"), "\n")
renv <- file.path(Sys.getenv("R_USER"), ".Renviron")
cat("Renviron exists:", file.exists(renv), renv, "\n")
if (file.exists(renv)) cat(paste(readLines(renv, warn = FALSE), collapse = "\n"), "\n")
