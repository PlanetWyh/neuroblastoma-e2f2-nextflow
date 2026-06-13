import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr

parser = argparse.ArgumentParser()
parser.add_argument("--input_file", required=True)
parser.add_argument("--output_dir", required=True)
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

df = pd.read_csv(args.input_file, index_col=0)

#plot 1: E2F2 score by death from disease
plot_df = df.dropna(subset=["death from disease", "E2F2_score"]).copy()

survived = plot_df.loc[plot_df["death from disease"] == 0, "E2F2_score"]
died = plot_df.loc[plot_df["death from disease"] == 1, "E2F2_score"]

stat, pval = mannwhitneyu(survived, died, alternative="two-sided")

plt.figure(figsize=(5, 5))
plt.boxplot(
    [survived, died],
    tick_labels=["Survived", "Died"]
)
plt.ylabel("E2F2 ssGSEA score")
plt.title(f"GSE49711: E2F2 score by outcome\nMann–Whitney p={pval:.3e}")
plt.tight_layout()
plt.savefig(
    os.path.join(args.output_dir, "GSE49711_E2F2_by_death.png"),
    dpi=300,
    bbox_inches="tight"
)
plt.close()


#plot 2: E2F2 score by MYCN status
if "mycn status" in df.columns:
    mycn_df = df.dropna(subset=["mycn status", "E2F2_score"]).copy()

    non_amp = mycn_df.loc[mycn_df["mycn status"] == 0, "E2F2_score"]
    amp = mycn_df.loc[mycn_df["mycn status"] == 1, "E2F2_score"]

    if len(non_amp) > 0 and len(amp) > 0:
        stat, pval = mannwhitneyu(non_amp, amp, alternative="two-sided")

        plt.figure(figsize=(5, 5))
        plt.boxplot(
            [non_amp, amp],
            tick_labels=["MYCN not amplified", "MYCN amplified"]
        )
        plt.ylabel("E2F2 ssGSEA score")
        plt.title(f"GSE49711: E2F2 score by MYCN status\nMann–Whitney p={pval:.3e}")
        plt.tight_layout()
        plt.savefig(
            os.path.join(args.output_dir, "GSE49711_E2F2_by_MYCN.png"),
            dpi=300,
            bbox_inches="tight"
        )
        plt.close()


#plot 3: E2F2 score by risk group
if "high risk" in df.columns:
    risk_df = df.dropna(subset=["high risk", "E2F2_score"]).copy()

    low_risk = risk_df.loc[risk_df["high risk"] == 0, "E2F2_score"]
    high_risk = risk_df.loc[risk_df["high risk"] == 1, "E2F2_score"]

    if len(low_risk) > 0 and len(high_risk) > 0:
        stat, pval = mannwhitneyu(low_risk, high_risk, alternative="two-sided")

        plt.figure(figsize=(5, 5))
        plt.boxplot(
            [low_risk, high_risk],
            tick_labels=["Non-high risk", "High risk"]
        )
        plt.ylabel("E2F2 ssGSEA score")
        plt.title(f"GSE49711: E2F2 score by risk group\nMann–Whitney p={pval:.3e}")
        plt.tight_layout()
        plt.savefig(
            os.path.join(args.output_dir, "GSE49711_E2F2_by_high_risk.png"),
            dpi=300,
            bbox_inches="tight"
        )
        plt.close()


#plot 4: E2F2 score vs age
if "age at diagnosis" in df.columns:
    age_df = df.dropna(subset=["age at diagnosis", "E2F2_score"]).copy()

    rho, pval = spearmanr(age_df["age at diagnosis"], age_df["E2F2_score"])

    plt.figure(figsize=(5, 5))
    plt.scatter(age_df["age at diagnosis"], age_df["E2F2_score"], alpha=0.7)
    plt.xlabel("Age at diagnosis")
    plt.ylabel("E2F2 ssGSEA score")
    plt.title(f"GSE49711: E2F2 vs age\nSpearman r={rho:.2f}, p={pval:.3e}")
    plt.tight_layout()
    plt.savefig(
        os.path.join(args.output_dir, "GSE49711_E2F2_vs_age.png"),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()


print("Saved plots to:", args.output_dir)