from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import os


parser = argparse.ArgumentParser()

parser.add_argument("--input_file", required=True)
parser.add_argument("--output_dir", required=True)

args = parser.parse_args()

merged = pd.read_csv(args.input_file)
os.makedirs(args.output_dir, exist_ok=True)

plot_path = os.path.join(args.output_dir, "km_e2f2_survival.png")
stats_path = os.path.join(args.output_dir, "logrank_results.csv")


results = []


def plot_km(ax, data, title, comparison_name):
    kmf = KaplanMeierFitter()

    high = data[data["E2F2_group"] == "High"]
    low = data[data["E2F2_group"] == "Low"]

    if len(high) == 0 or len(low) == 0:
        ax.set_title(f"{title}\nNot enough groups")
        return

    kmf.fit(
        high["duration"],
        high["event"],
        label=f"E2F2 High (n={len(high)})"
    )
    kmf.plot_survival_function(ax=ax, ci_show=True)

    kmf.fit(
        low["duration"],
        low["event"],
        label=f"E2F2 Low (n={len(low)})"
    )
    kmf.plot_survival_function(ax=ax, ci_show=True)

    lr = logrank_test(
        high["duration"],
        low["duration"],
        high["event"],
        low["event"]
    )

    ax.set_title(f"{title}\np={lr.p_value:.3g}")
    ax.set_xlabel("Days")
    ax.set_ylabel("Survival probability")

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


fig, axes = plt.subplots(1, 3, figsize=(18, 6))

plot_km(
    axes[0],
    merged,
    "All patients",
    "all_patients"
)

high_risk = merged[
    merged["diagnoses.cog_neuroblastoma_risk_group"] == "High Risk"
]

plot_km(
    axes[1],
    high_risk,
    "High Risk only",
    "high_risk_only"
)

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