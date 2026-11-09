# -*- coding: utf-8 -*-
"""
02_figure3_plot.py
Figure 3 visualisation panel (three subplots) for the DML analysis:
  (1) cross-fitted residual scatter (residual Y vs residual X);
  (2) original beta scatter (outcome vs exposure);
  (3) OLS vs DML estimate comparison.
A separate sensitivity forest compares the Lasso and random-forest learners.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

OUT_DIR = "D:/研究生科研小组/小论文/MFI-cancer-meta-MR-DML/output_results/DML"
FIG_DIR = os.path.join(OUT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

sns.set_style("whitegrid")


def plot_figure3(data, residuals, dml_theta, ols_theta):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # Panel 1: residual scatter with the DML slope
    ax = axes[0]
    rx, ry = residuals
    ax.scatter(rx, ry, s=18, alpha=0.7, color="#2c7fb8")
    xs = np.linspace(rx.min(), rx.max(), 100)
    ax.plot(xs, dml_theta * xs, color="#d7301f", lw=2,
            label=f"DML slope = {dml_theta:.3f}")
    ax.set_xlabel("Residualised exposure beta")
    ax.set_ylabel("Residualised outcome beta")
    ax.set_title("(A) Cross-fitted residuals")
    ax.legend(frameon=True)

    # Panel 2: original beta scatter
    ax = axes[1]
    ax.scatter(data["beta_exp_aligned"], data["beta_out"],
               s=18, alpha=0.7, color="#41ab5d")
    ax.set_xlabel("Exposure beta (aligned)")
    ax.set_ylabel("Outcome beta")
    ax.set_title("(B) Original summary effects")

    # Panel 3: OLS vs DML point estimate
    ax = axes[2]
    labels = ["OLS", "DML (LassoCV)"]
    vals = [ols_theta, dml_theta]
    ax.bar(labels, vals, color=["#9e9ac8", "#6a51a3"], width=0.55)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel("Estimated logOR")
    ax.set_title("(C) OLS vs DML")

    plt.tight_layout()
    out = os.path.join(FIG_DIR, "Figure3_DML_panel.png")
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    return out


def plot_sensitivity_forest(estimates):
    """Forest of OR (95% CI) for the Lasso and random-forest learners."""
    fig, ax = plt.subplots(figsize=(7, 3))
    for i, row in enumerate(estimates.itertuples()):
        ax.plot([row.OR_LCI, row.OR_UCI], [i, i], color="#3182bd", lw=2)
        ax.plot(row.OR, i, "s", color="#08519c", ms=8)
    ax.axvline(1, ls="--", color="grey")
    ax.set_yticks(range(len(estimates)))
    ax.set_yticklabels(estimates["learner"])
    ax.invert_yaxis()
    ax.set_xlabel("Odds ratio (95% CI)")
    ax.set_title("DML sensitivity: LassoCV vs Random Forest")
    plt.tight_layout()
    out = os.path.join(FIG_DIR, "DML_sensitivity_forest.png")
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    return out


def main():
    # Load the saved estimates; the residual arrays are produced by
    # 01_dml_main.py and may be exported there for full reproducibility.
    estimates = pd.read_csv(os.path.join(OUT_DIR, "DML_estimates.csv"))
    plot_sensitivity_forest(estimates)
    print("Sensitivity forest saved.")


if __name__ == "__main__":
    main()
