# =============================================================================
# 02_forward_MR.R
# Forward two-sample MR: male factor infertility traits -> cancers.
# Same pipeline and locked parameters as the reverse analysis; only the
# exposure/outcome roles are swapped. Uses the OpenGWAS online API.
#
# Locked parameters (do not change):
#   p1 = 1e-5; clump r2 = 0.001, kb = 10000; harmonise action = 2;
#   remove palindromic/ambiguous SNPs and MAF < 0.01;
#   methods: IVW, MR-Egger, weighted median, simple mode, weighted mode.
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
})

# ---- Token: read from environment, never hard-code ----
my_token <- Sys.getenv("OPENGWAS_JWT", unset = Sys.getenv("OPENGWAS_TOKEN"))
if (!nzchar(my_token)) stop("Set the OPENGWAS_JWT environment variable before running.")
Sys.setenv(OPENGWAS_JWT = my_token)

p1_thr <- 1e-5; clump_r2 <- 0.001; clump_kb <- 10000; maf_thr <- 0.01

# ---- Exposure list (male factor infertility traits) ----
exposure_ids <- c("finngen_R12_N14_MALEINFERT",
                  "finn-b-R18_ABNORMAL_SPERMATOZ")

# ---- Cancer outcome list (same 32 datasets used in the reverse analysis) ----
cancer_map <- data.table::fread(
  "D:/研究生科研小组/小论文/MFI-cancer-meta-MR-DML/output_results/MR/cancer_map.csv")
# If the map is unavailable, define outcome IDs directly:
# cancer_map <- data.frame(exposure_id = c("ieu-b-...", "finn-b-..."),
#                         cancer_subgroup = c("...", "..."))

out_dir <- "D:/研究生科研小组/小论文/MFI-cancer-meta-MR-DML/output_results/MR"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

all_res <- list()
for (exp_id in exposure_ids) {
  # 1. Exposure instruments
  exp_dat <- extract_instruments(outcomes = exp_id, p1 = p1_thr,
                                 clump_r2 = clump_r2, clump_kb = clump_kb)
  for (i in seq_len(nrow(cancer_map))) {
    out_id <- cancer_map$exposure_id[i]
    message("Forward: ", exp_id, " -> ", out_id)
    out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = out_id)
    har <- harmonise_data(exp_dat, out_dat, action = 2)
    har <- har[har$mr_keep, ]
    har <- har[!har$palindromic, ]
    har <- har[!is.na(har$eaf.outcome) & har$eaf.outcome >= maf_thr &
               har$eaf.outcome <= 1 - maf_thr, ]
    if (nrow(har) < 3) next

    mr <- mr(har, method_list = c("mr_ivw", "mr_egger_regression",
                                 "mr_weighted_median", "mr_simple_mode",
                                 "mr_weighted_mode"))
    het <- mr_heterogeneity(har)
    ple <- mr_pleiotropy_test(har)
    mr$id.exposure <- exp_id
    mr$id.outcome <- out_id
    all_res[[length(all_res) + 1]] <- mr
  }
}

res <- data.table::rbindlist(all_res, fill = TRUE)
data.table::fwrite(res, file.path(out_dir, "forward_MR_all_methods.csv"), bom = TRUE)
cat("Forward MR complete. Rows:", nrow(res), "\n")
