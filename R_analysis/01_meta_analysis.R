# =============================================================================
# 01_meta_analysis.R
# Meta-analysis of adjusted hazard ratios (aHR) from observational cohort studies.
# Reconstructed from the original analysis; produces pooled estimates, heterogeneity,
# leave-one-out sensitivity, publication bias / trim-and-fill, and subgroup RVE.
# Packages: readxl, meta, metafor, robumeta, clubSandwich
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Load packages and read the extracted data table
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(meta)
  library(metafor)
})

# Edit this path to the extracted aHR table (Excel or CSV)
input_xlsx <- "D:/研究生科研小组/小论文/meta+MR/成功版本/3 meta分析/数据提取.xlsx"
mydata <- read_excel(input_xlsx)

# Coerce the effect estimate and its 95% CI bounds to numeric
mydata$aHR_num <- as.numeric(mydata$`aHR`)
mydata$LCI     <- as.numeric(mydata$`Lower limit of 95% confidence interval`)
mydata$UCI     <- as.numeric(mydata$`95% confidence interval upper limit`)

# Log-transform: effect on log scale and its standard error
mydata$log_HR  <- log(mydata$aHR_num)
mydata$log_LCI <- log(mydata$LCI)
mydata$log_UCI <- log(mydata$UCI)
mydata$SE      <- (mydata$log_UCI - mydata$log_LCI) / (2 * 1.96)

# Study label, e.g. "Author Year"
mydata$study_label <- paste(gsub("\n", " ", mydata$Author), mydata$Year)

# Output directory
out_dir <- "D:/研究生科研小组/小论文/MFI-cancer-meta-MR-DML/output_results/meta"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Generic inverse-variance meta-analysis (fixed + random, REML tau^2)
# -----------------------------------------------------------------------------
meta_hr <- metagen(TE = log_HR,
                   seTE = SE,
                   studlab = study_label,
                   data = mydata,
                   sm = "HR",
                   fixed = TRUE,
                   random = TRUE,
                   method.tau = "REML")
summary(meta_hr)

# Per-study effects and weights
studies_res <- data.frame(
  Study = meta_hr$studlab,
  aHR = exp(meta_hr$TE),
  LCI = exp(meta_hr$TE - 1.96 * meta_hr$seTE),
  UCI = exp(meta_hr$TE + 1.96 * meta_hr$seTE),
  Log_aHR = meta_hr$TE,
  SE = meta_hr$seTE,
  Weight_Fixed = meta_hr$w.fixed / sum(meta_hr$w.fixed) * 100,
  Weight_Random = meta_hr$w.random / sum(meta_hr$w.random) * 100)

# Pooled summary and heterogeneity
pooled_res <- data.frame(
  Model = c("Fixed effect model", "Random effects model"),
  aHR = c(exp(meta_hr$TE.fixed), exp(meta_hr$TE.random)),
  LCI = c(exp(meta_hr$lower.fixed), exp(meta_hr$lower.random)),
  UCI = c(exp(meta_hr$upper.fixed), exp(meta_hr$upper.random)),
  P_value = c(meta_hr$pval.fixed, meta_hr$pval.random),
  I2 = meta_hr$I2 * 100,
  Tau2 = meta_hr$tau2,
  P_heterogeneity = meta_hr$pval.Q)

write.csv(studies_res, file.path(out_dir, "meta_individual_studies.csv"), row.names = FALSE)
write.csv(pooled_res, file.path(out_dir, "meta_pooled_summary.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 2. Leave-one-out sensitivity analysis
# -----------------------------------------------------------------------------
loo <- metainf(meta_hr, pooled = "random")
write.csv(data.frame(Study = loo$studlab,
                    aHR = exp(loo$TE.random),
                    LCI = exp(loo$lower.random),
                    UCI = exp(loo$upper.random)),
          file.path(out_dir, "meta_leave_one_out.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Publication bias: funnel plot, Egger test and trim-and-fill
# -----------------------------------------------------------------------------
egger <- metabias(meta_hr, method.bias = "linreg")
tf <- trimfill(meta_hr)
write.csv(data.frame(Egger_p = egger$p.value,
                    Trimfill_aHR = exp(tf$TE.random),
                    Trimfill_LCI = exp(tf$lower.random),
                    Trimfill_UCI = exp(tf$upper.random)),
          file.path(out_dir, "meta_bias_trimfill.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Subgroup analysis with cluster-robust variance estimation (RVE)
# Multiple non-independent outcomes can share one cohort; robumeta handles clustering.
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(robumeta); library(clubSandwich) })
sub <- read_excel("D:/研究生科研小组/小论文/meta+MR/成功版本/3 meta分析/亚组分析表.xlsx", sheet = "Sheet1")
colnames(sub)[colnames(sub) == "Subgroups"] <- "subgroup"
colnames(sub)[colnames(sub) == "Study"] <- "study"
colnames(sub)[colnames(sub) == "Diseases"] <- "disease"
colnames(sub)[colnames(sub) == "aHR"] <- "HR"
colnames(sub)[colnames(sub) == "Lower limit of 95% confidence interval"] <- "lo"
colnames(sub)[colnames(sub) == "95% confidence interval upper limit"] <- "hi"

# Numeric study id (cluster variable)
sub$studyid <- NA
sub$studyid[grepl("Eisenberg", sub$study)] <- 1
sub$studyid[grepl("Zarén", sub$study)] <- 2
sub$studyid[grepl("Al-Jebari", sub$study)] <- 3
sub$studyid[grepl("Hanson", sub$study)] <- 4
sub$studyid[grepl("Walsh", sub$study)] <- 5
sub <- sub[!is.na(sub$studyid), ]
# Remove nested/early-onset sub-types
sub <- sub[!(sub$studyid == 5 & grepl("low-grade|high-grade", sub$disease, ignore.case = TRUE)), ]
sub <- sub[!(sub$studyid == 3 & grepl("early-onset", sub$disease, ignore.case = TRUE)), ]

sub$yi <- log(sub$HR)
sub$se <- (log(sub$hi) - log(sub$lo)) / (2 * 1.96)
sub$vi <- sub$se^2

robu_fit <- robu(yi ~ 0 + subgroup, data = sub,
                 studynum = studyid, var.eff.size = vi, small = TRUE)
print(robu_fit)

reg <- robu_fit$reg.results
subgroup_res <- data.frame(Subgroup = reg$terms,
                           logHR = reg$beta,
                           se_robust = reg$se,
                           p_robust = reg$pval,
                           HR = exp(reg$beta),
                           CI_lower = exp(reg$beta - 1.96 * reg$se),
                           CI_upper = exp(reg$beta + 1.96 * reg$se))
write.csv(subgroup_res, file.path(out_dir, "meta_subgroup_RVE.csv"), row.names = FALSE)
cat("Meta-analysis complete. Outputs saved to:", out_dir, "\n")
