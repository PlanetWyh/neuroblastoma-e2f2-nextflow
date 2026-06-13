import argparse
import pandas as pd


parser = argparse.ArgumentParser()

parser.add_argument("--clinical_file", required=True)
parser.add_argument("--e2f2_scores_file", required=True)
parser.add_argument("--output_file", required=True)
parser.add_argument("--summary_file", required=True)

parser.add_argument("--mycn_file", required=False, default=None)
parser.add_argument("--min_group_size", required=False, type=int, default=15)

args = parser.parse_args()

#load clinical data
clinical = pd.read_csv(args.clinical_file, sep="\t")

cols_needed = [
    "cases.submitter_id",
    "demographic.vital_status",
    "demographic.days_to_death",
    "diagnoses.days_to_last_follow_up",
    "diagnoses.age_at_diagnosis",
    "diagnoses.inss_stage",
    "diagnoses.cog_neuroblastoma_risk_group",
]

missing_cols = [c for c in cols_needed if c not in clinical.columns]
if missing_cols:
    raise ValueError(f"Missing clinical columns: {missing_cols}")

clinical_clean = clinical[cols_needed].drop_duplicates(
    subset="cases.submitter_id"
).copy()

clinical_clean["demographic.days_to_death"] = pd.to_numeric(
    clinical_clean["demographic.days_to_death"],
    errors="coerce"
)

clinical_clean["diagnoses.days_to_last_follow_up"] = pd.to_numeric(
    clinical_clean["diagnoses.days_to_last_follow_up"],
    errors="coerce"
)

clinical_clean["duration"] = clinical_clean["demographic.days_to_death"].fillna(
    clinical_clean["diagnoses.days_to_last_follow_up"]
)

clinical_clean["event"] = (
    clinical_clean["demographic.vital_status"] == "Dead"
).astype(int)

clinical_clean = clinical_clean.dropna(subset=["duration"])


#load E2F2 scores
e2f2 = pd.read_csv(args.e2f2_scores_file)

e2f2.columns = ["sample_id", "E2F2_score"]

# TARGET-30-XXXXXX-01A-01R → TARGET-30-XXXXXX
e2f2["cases.submitter_id"] = e2f2["sample_id"].str[:16]


#merge clinical + E2F2
merged = clinical_clean.merge(
    e2f2,
    on="cases.submitter_id",
    how="inner"
)


#optional MYCN status
mycn_status_available = False
mycn_usable_for_stratification = False

if args.mycn_file is not None:
    mycn_data = pd.read_csv(args.mycn_file, sep="\t")

    def get_mycn_status(row):
        for slot in [
            "follow_ups.0.molecular_tests.0",
            "follow_ups.1.molecular_tests.0",
        ]:
            gene = row.get(f"{slot}.gene_symbol")
            result = row.get(f"{slot}.test_result")

            if gene == "MYCN" and result in ["Amplified", "Not Amplified"]:
                return result

        return None

    mycn_data["mycn_status"] = mycn_data.apply(get_mycn_status, axis=1)

    mycn_clean = mycn_data[["submitter_id", "mycn_status"]].copy()
    mycn_clean = mycn_clean[mycn_clean["mycn_status"].notna()]

    merged = merged.merge(
        mycn_clean,
        left_on="cases.submitter_id",
        right_on="submitter_id",
        how="left"
    )

    amplified_n = int((merged["mycn_status"] == "Amplified").sum())
    not_amplified_n = int((merged["mycn_status"] == "Not Amplified").sum())

    if amplified_n >= args.min_group_size and not_amplified_n >= args.min_group_size:
        mycn_usable_for_stratification = True
else:
    merged["mycn_status"] = "Not available"

mycn_counts = merged["mycn_status"].value_counts(dropna=False)
print("MYCN status counts:")
print(mycn_counts)

if "Amplified" in mycn_counts.index:
    if mycn_counts["Amplified"] < args.min_group_size:
        print(
            "WARNING: MYCN-amplified group is too small for reliable stratified survival analysis. "
            "Use COG risk group instead."
        )

#Add E2F2 high/low groups
median_score = merged["E2F2_score"].median()

merged["E2F2_group"] = merged["E2F2_score"].apply(
    lambda x: "High" if x >= median_score else "Low"
)



#Save output table
merged.to_csv(args.output_file, index=False)

#summary
with open(args.summary_file, "w") as f:
    f.write("Survival input preparation summary\n")
    f.write("==================================\n\n")

    f.write(f"Clinical patients after cleaning: {clinical_clean.shape[0]}\n")
    f.write(f"Patients with E2F2 scores merged: {merged.shape[0]}\n")
    f.write(f"Events/deaths: {int(merged['event'].sum())}\n")
    f.write(f"Censored: {int((merged['event'] == 0).sum())}\n\n")

    f.write("E2F2 group counts:\n")
    f.write(str(merged["E2F2_group"].value_counts()))
    f.write("\n\n")

    f.write("COG risk group counts:\n")
    f.write(str(merged["diagnoses.cog_neuroblastoma_risk_group"].value_counts(dropna=False)))
    f.write("\n\n")

    if args.mycn_file is not None:
        f.write("MYCN status counts:\n")
        f.write(str(merged["mycn_status"].value_counts(dropna=False)))
        f.write("\n\n")

        f.write(f"Minimum group size threshold: {args.min_group_size}\n")
        f.write(f"MYCN usable for stratified survival analysis: {mycn_usable_for_stratification}\n\n")

        if not mycn_usable_for_stratification:
            f.write(
                "MYCN stratification was not used as the primary stratification variable "
                "because at least one MYCN subgroup was below the minimum group size threshold. "
                "COG neuroblastoma risk group should be used as the primary clinical stratifier instead.\n"
            )
    else:
        f.write("No MYCN file was provided.\n")

print(f"Saved survival table to {args.output_file}")
print(f"Saved summary to {args.summary_file}")