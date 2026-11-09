# =============================================================================
# 04_sensitivity_FDR_MDOR.R
# Post-processing of the 2SMR IVW results:
#   (1) Benjamini-Hochberg FDR correction across ALL primary IVW tests;
#   (2) Minimum detectable odds ratio (MDOR) for 80% power at alpha = 0.05.
# Reads the combined IVW result table; writes independent outputs, never edits source.
# =============================================================================

library(data.table)

# ---- Paths (edit if needed) ----
src_file <- "D:/研究生科研小组/小论文/MFI-cancer-meta-MR-DML/output_results/MR/MR_ivw_all.csv"
out_dir  <- "D:/研究生科研小组/小论文/MFI-cancer-meta-MR-DML/output_results/MR"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ivw <- data.table::fread(src_file, encoding = "UTF-8")
# Expected columns: direction, exposure_id, outcome_id, nsnp, ivw_beta, ivw_se, ivw_pval

# -----------------------------------------------------------------------------
# 1. Benjamini-Hochberg FDR
# IMPORTANT: the correction set must cover ALL primary IVW tests (do not first
# subset to significant rows). Monotonicity is enforced from the LARGEST p upward.
# -----------------------------------------------------------------------------
bh_fdr <- function(p) {
  n <- length(p)
  o <- order(p, decreasing = TRUE)        # start from the largest p
  q <- numeric(n)
  prev <- 1
  for (rank_desc in seq_len(n)) {
    idx <- o[rank_desc]
    asc_rank <- n - rank_desc + 1         # ascending rank of this p-value
    val <- p[idx] * n / asc_rank
    val <- min(val, prev)                 # enforce monotonicity, cap at 1
    q[idx] <- val
    prev <- val
  }
  pmin(q, 1)
}

ivw[, bh_fdr := bh_fdr(ivw_pval)]

# -----------------------------------------------------------------------------
# 2. MDOR = exp(2.8016 * ivw_se)
# 2.8016 = Z(1-alpha/2) + Z(power) = 1.96 + 0.8416
# -----------------------------------------------------------------------------
ivw[, MDOR := exp(2.8016 * ivw_se)]

# -----------------------------------------------------------------------------
# 3. Labels following the manuscript rule (strict < 0.05)
# -----------------------------------------------------------------------------
ivw[, Suggestion := ifelse(
  ivw_pval < 0.05,
  "Statistically significant \u2014 post-hoc power is for reference only and not recommended for formal reporting",
  "Report MDOR (non-significant, under-powered)")]

ivw[, Tip := ifelse(ivw_pval < 0.05, "/",
                    "limited statistical power / under-powered")]

# Special handling: effect sizes of implausible magnitude (|beta| in the tens)
# are not converted to OR/MDOR; only the under-powered attribute is retained.
ivw[abs(ivw_beta) > 10, `:=`(MDOR = NA_real_,
                             Suggestion = "Report MDOR (non-significant, under-powered)",
                             Tip = "limited statistical power / under-powered; implausible effect size")]

data.table::fwrite(ivw, file.path(out_dir, "MR_ivw_FDR_MDOR.csv"), bom = TRUE)
cat("Rows:", nrow(ivw),
    "| nominal P<0.05:", sum(ivw$ivw_pval < 0.05),
    "| FDR<0.05:", sum(ivw$bh_fdr < 0.05), "\n")
