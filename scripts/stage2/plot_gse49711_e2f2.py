import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr

#this script takes the merged clinical + E2F2 score table from validate_e2f2_gse49711.py and produces four boxplots
#each boxplot compares E2F2 ssGSEA scores across a different group:
#1st group: outcome (survived vs. died from disease)
#2nd group: MYCN amplification status
#3rd group: risk group (high vs. non-high risk)
#4th group: age at diagnosis (scatter + Spearman correlation)

#statistical test used throughout is the Mann-Whitney U test (Wilcoxon rank-sum test) rather than a t-test
#because E2F2 ssGSEA scores are not guaranteed to be normally distributed and Mann-Whitney tests whether one group tends to have
#higher scores than the other without any assumptions about distribution shape, 
 
parser = argparse.ArgumentParser()
parser.add_argument("--input_file", required=True)
parser.add_argument("--output_dir", required=True)
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

#index_col=0 restores the GSM accession as the row index (it was saved that way by validate_e2f2_gse49711.py)
df = pd.read_csv(args.input_file, index_col=0)

#plot 1: E2F2 score by death from disease (survived vs. died from disease)
#dropna to ensure patients with missing outcome data don't contaminate some group
#"death from disease" column is 0/1 (set to numeric in validate_e2f2_gse49711.py)
plot_df = df.dropna(subset=["death from disease", "E2F2_score"]).copy()

survived = plot_df.loc[plot_df["death from disease"] == 0, "E2F2_score"]
died = plot_df.loc[plot_df["death from disease"] == 1, "E2F2_score"]

#test for ANY difference in E2F2 score between groups, not specifically whether "died" is higher
#(one-sided test would give a lower p-value but requires pre-specifying the direction)
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


#plot 2: E2F2 score by MYCN amplification status
#MYCN amplification is the most established marker of aggressive neuroblastoma and was hypothesized to correlate with E2F2 activity
#given their functional linkage (MYCN drives E2F target gene activation)
#column is 0 = not amplified, 1 = amplified
if "mycn status" in df.columns:
    mycn_df = df.dropna(subset=["mycn status", "E2F2_score"]).copy()

    non_amp = mycn_df.loc[mycn_df["mycn status"] == 0, "E2F2_score"]
    amp = mycn_df.loc[mycn_df["mycn status"] == 1, "E2F2_score"]
    #only plot if both groups are non-empty (to guard in case the dataset happens to have no amplified or no non-amplified cases)
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
#"high risk" column values are  0 = non-high-risk, 1 = high-risk.
#(this is alternative for the COG risk stratification used in TARGET, so a consistent direction across both cohorts
#strengthens the claim that E2F2 is specifically elevated in aggressive disease
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


#plot 4: E2F2 score vs age at diagnosis (scatter + correlation)
#age is a continuous variable so it gets a scatter plot and Spearman correlation
#Spearman is used (not Pearson) because it measures monotonic association without assuming linearity or normality, appropriate
#given the skewed age distribution visible in the scatterplot.
#NOTE: this plot was not included in the final report
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
