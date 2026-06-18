from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import os

#this script takes the survival-ready table from prepare_survival_data.py and produces Kaplan-Meier survival curves 
#comparing E2F2 high vs low patients, with log-rank test p-values
#three comparisons are run: all patients + High Risk patients only + Low/Intermediate Risk patients only (too small num of patients here to divide further)
#this stratification by COG risk group tests whether the E2F2 signal is specific to the most aggressive disease subgroup rather than a generic marker across all severities

parser = argparse.ArgumentParser()

parser.add_argument("--input_file", required=True)
parser.add_argument("--output_dir", required=True)

args = parser.parse_args()

merged = pd.read_csv(args.input_file)
os.makedirs(args.output_dir, exist_ok=True)

plot_path = os.path.join(args.output_dir, "km_e2f2_survival.png")
stats_path = os.path.join(args.output_dir, "logrank_results.csv")

#accumulates one results row per plot_km() call, then saved as csv at the end of the script
results = []

#core plotting function: draws one KM panel on a given matplotlib axes object and records the log-rank test result
def plot_km(ax, data, title, comparison_name):
    """
    ax:               matplotlib xxes to draw on
    data:             DataFrame subset for this comparison (already filtered to the relevant risk group)
    title:            panel title
    comparison_name:  string key stored in the results csv
    """
    kmf = KaplanMeierFitter()

    high = data[data["E2F2_group"] == "High"]
    low = data[data["E2F2_group"] == "Low"]

    #if either group is empty (e.g. a risk subgroup has no "high" patients after filtering), skip plotting and not crash
    #(in practice this shouldn't happen with a median split but just to be on a save side)

    if len(high) == 0 or len(low) == 0:
        ax.set_title(f"{title}\nNot enough groups")
        return
        
    #fit and plot the E2F2-high survival curve.
    #KaplanMeierFitter.fit() takes:
    #durations: time to event or censoring (days here)
    #event_observed: 1 = death occurred, 0 = censored (still alive or lost to follow-up at last contact)
    #ci_show=True adds the 95% confidence interval shaded region
    kmf.fit(
        high["duration"],
        high["event"],
        label=f"E2F2 High (n={len(high)})"
    )
    kmf.plot_survival_function(ax=ax, ci_show=True)

    #fit and plot the E2F2-Low survival curve on the same axes
    kmf.fit(
        low["duration"],
        low["event"],
        label=f"E2F2 Low (n={len(low)})"
    )
    kmf.plot_survival_function(ax=ax, ci_show=True)

    #log-rank test: tests whether the two survival curves are statistically distinguishable
    #the null hypothesis is that the two groups have the same underlying survival function
    #p < 0.05 is taken as evidence of a significant difference
    lr = logrank_test(
        high["duration"],
        low["duration"],
        high["event"],
        low["event"]
    )

    ax.set_title(f"{title}\np={lr.p_value:.3g}")
    ax.set_xlabel("Days")
    ax.set_ylabel("Survival probability")
    #store full statistics for the csv output
    results.append({
        "comparison": comparison_name,
        "n_total": len(data),
        "n_high": len(high),
        "n_low": len(low),
        "events_high": int(high["event"].sum()),
        "events_low": int(low["event"].sum()),
        "p_value": lr.p_value,
        "test_statistic": lr.test_statistic
    })

#build a 1x3 figure and populate each panel
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

#plot 1: all patients regardless of risk group
plot_km(
    axes[0],
    merged,
    "All patients",
    "all_patients"
)

#plot 2: High Risk patients only
high_risk = merged[
    merged["diagnoses.cog_neuroblastoma_risk_group"] == "High Risk"
]

plot_km(
    axes[1],
    high_risk,
    "High Risk only",
    "high_risk_only"
)

#plot 3: Low + Intermediate Risk patients combined (n=33 combined already limits power)
low_intermediate = merged[
    merged["diagnoses.cog_neuroblastoma_risk_group"].isin(
        ["Low Risk", "Intermediate Risk"]
    )
]

plot_km(
    axes[2],
    low_intermediate,
    "Low + Intermediate Risk",
    "low_intermediate_risk"
)

plt.suptitle(
    "E2F2 Regulon Score — Kaplan-Meier Survival Analysis",
    fontsize=14,
    fontweight="bold"
)

plt.tight_layout()
plt.savefig(plot_path, dpi=300, bbox_inches="tight")
plt.close()

pd.DataFrame(results).to_csv(stats_path, index=False)

print(f"Saved plot to {plot_path}")
print(f"Saved statistics to {stats_path}")
