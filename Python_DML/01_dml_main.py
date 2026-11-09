# -*- coding: utf-8 -*-
"""
01_dml_main.py
Summary-level double machine learning (DML) with cross-fitting.
Primary learner: LassoCV; sensitivity learner: RandomForestRegressor.
Converted from the original Jupyter notebook (DML.ipynb); parameters are unchanged.

Workflow:
  load FinnGen .gz -> select IVs at p < 1e-5 -> harmonise by chrom+pos / allele
  -> 5-fold cross-fitting nuisance models -> orthogonal theta estimate
  -> SE, z, p, OR, 95% CI; random-forest sensitivity.
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import KFold
from sklearn.linear_model import LassoCV
from sklearn.ensemble import RandomForestRegressor

# ----------------------------------------------------------------------
# 0. Paths and parameters
# ----------------------------------------------------------------------
# The final analysis uses FinnGen R12 (the exploratory notebook referenced R5).
EXPOSURE_GZ = "D:/研究生科研小组/小论文/GWAS数据库/finn数据/finngen_R12_..._exposure.gz"
OUTCOME_GZ = "D:/研究生科研小组/小论文/GWAS数据库/finn数据/finngen_R12_..._outcome.gz"
OUT_DIR = "D:/研究生科研小组/小论文/MFI-cancer-meta-MR-DML/output_results/DML"
os.makedirs(OUT_DIR, exist_ok=True)

P_THRESHOLD = 1e-5
N_FOLDS = 5
SEED = 1234


# ----------------------------------------------------------------------
# 1. Load summary statistics and select instrumental variables
# ----------------------------------------------------------------------
def load_sumstats(path):
    """Read a FinnGen .gz summary-statistic file."""
    df = pd.read_csv(path, compression="gzip", sep="\t")
    return df


def select_iv(df, p_col="pval", threshold=P_THRESHOLD):
    """Keep SNPs meeting the p < 1e-5 instrument threshold."""
    return df[df[p_col] < threshold].copy()


def harmonise(exp, out):
    """Merge on chromosome/position and align alleles; build beta/SE columns."""
    merged = pd.merge(
        exp, out,
        left_on=["chrom", "pos"], right_on=["chrom", "pos"],
        suffixes=("_exp", "_out"), how="inner")
    # Flip the exposure effect when the effect alleles do not match.
    flip = merged["alt_exp"] != merged["alt_out"]
    merged["beta_exp_aligned"] = np.where(flip, -merged["beta_exp"], merged["beta_exp"])
    keep = ["beta_exp_aligned", "sebeta_exp", "beta_out", "sebeta_out"]
    merged = merged.dropna(subset=keep)
    return merged.reset_index(drop=True)


# ----------------------------------------------------------------------
# 2. Cross-fitting DML
# ----------------------------------------------------------------------
def cross_fit_dml(data, nuisance_model, n_folds=N_FOLDS, seed=SEED):
    """
    Cross-fitted nuisance models. Each fold is predicted by models fit on the
    other folds, so residuals are out-of-sample and avoid overfitting bias.
    """
    x = data["beta_exp_aligned"].values.reshape(-1, 1)
    # Covariates available at summary level (here: the aligned exposure SE as a
    # weight feature). The design is intentionally small and transparent.
    w = data[["sebeta_exp", "sebeta_out"]].values
    y = data["beta_out"].values

    kf = KFold(n_splits=n_folds, shuffle=True, random_state=seed)
    res_x = np.zeros_like(x.ravel())
    res_y = np.zeros_like(y)

    for train_idx, test_idx in kf.split(w):
        mx = nuisance_model()
        my = nuisance_model()
        mx.fit(w[train_idx], x[train_idx].ravel())
        my.fit(w[train_idx], y[train_idx])
        res_x[test_idx] = x[test_idx].ravel() - mx.predict(w[test_idx])
        res_y[test_idx] = y[test_idx] - my.predict(w[test_idx])

    # Orthogonal (IVW-style) point estimate
    theta = float(np.sum(res_x * res_y) / np.sum(res_x ** 2))
    n = len(y)
    residual = res_y - theta * res_x
    se = float(np.sqrt(np.sum(residual ** 2) / (n - 1)) / np.sqrt(np.sum(res_x ** 2)))
    z = theta / se
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return {"theta": theta, "se": se, "z": z, "p": p,
            "OR": np.exp(theta),
            "OR_LCI": np.exp(theta - 1.96 * se),
            "OR_UCI": np.exp(theta + 1.96 * se)}


def lasso_model():
    # Matches the notebook: 3-fold internal CV, 5000 iterations
    return LassoCV(cv=3, max_iter=5000, random_state=SEED)


def rf_model():
    # Sensitivity learner: 200 trees, depth 5
    return RandomForestRegressor(n_estimators=200, max_depth=5, random_state=SEED)


# ----------------------------------------------------------------------
# 3. Run
# ----------------------------------------------------------------------
def main():
    exp = select_iv(load_sumstats(EXPOSURE_GZ))
    out = load_sumstats(OUTCOME_GZ)
    data = harmonise(exp, out)
    print("Harmonised instruments:", len(data))

    lasso = cross_fit_dml(data, lasso_model)
    rf = cross_fit_dml(data, rf_model)

    table = pd.DataFrame([
        {"learner": "LassoCV", **lasso},
        {"learner": "RandomForest", **rf}])
    table.to_csv(os.path.join(OUT_DIR, "DML_estimates.csv"), index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
