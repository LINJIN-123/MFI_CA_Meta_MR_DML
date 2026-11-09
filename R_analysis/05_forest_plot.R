# =============================================================================
# 05_forest_plot.R
# Forest and funnel plots for the meta-analysis and the two-sample MR.
# Small publication-ready figures (300 dpi PNG). Raw figure files are git-ignored
# under output_results; keep the scripts and any small result tables in the repo.
# =============================================================================

suppressPackageStartupMessages({
  library(meta)
  library(TwoSampleMR)
  library(data.table)
})

out_dir <- "D:/研究生科研小组/小论文/MFI-cancer-meta-MR-DML/output_results"
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Meta-analysis forest and funnel (re-run the meta object or load it) ----
# Recreate the meta object from 01_meta_analysis.R outputs.
studies <- data.table::fread(file.path(out_dir, "meta/meta_individual_studies.csv"))
meta_hr <- metagen(TE = Log_aHR, seTE = SE, studlab = Study,
                   data = studies, sm = "HR", method.tau = "REML")

png(file.path(fig_dir, "meta_forest.png"), width = 2400, height = 1500, res = 300)
forest(meta_hr, backtransf = TRUE, prediction = TRUE,
       col.square = "navy", col.diamond = "maroon")
dev.off()

png(file.path(fig_dir, "meta_funnel.png"), width = 1800, height = 1800, res = 300)
funnel(meta_hr)
dev.off()

# ---- 2. MR forest plot (one exposure, all cancer outcomes) ----
mr_res <- data.table::fread(file.path(out_dir, "MR/forward_MR_all_methods.csv"))
ivw <- mr_res[method == "Inverse variance weighted"]
png(file.path(fig_dir, "MR_forest_IVW.png"), width = 2400, height = 3000, res = 300)
# TwoSampleMR forest expects a harmonised-data + result object; a simple
# manual forest is drawn from the IVW estimates for portability.
beta <- ivw$b; se <- ivw$se
ord <- order(beta)
plot(beta[ord], seq_along(beta), xlim = range(c(beta - 1.96 * se, beta + 1.96 * se)),
     yaxt = "n", ylab = "", xlab = "IVW logOR (95% CI)", pch = 15)
arrows(beta[ord] - 1.96 * se[ord], seq_along(beta),
       beta[ord] + 1.96 * se[ord], seq_along(beta), length = 0.02, angle = 90, code = 3)
abline(v = 0, lty = 2)
axis(2, at = seq_along(beta), labels = ivw$id.outcome[ord], las = 1)
dev.off()

cat("Figures saved to:", fig_dir, "\n")
