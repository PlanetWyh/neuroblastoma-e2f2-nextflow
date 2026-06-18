import argparse
import pandas as pd

#this script combines three separate data sources: clinical/survival info + E2F2 ssGSEA scores and (optionally) MYCN amplification status
#it creates one tidy per-patient table ready for Kaplan-Meier survival analysis and human-readable summary of how many patients made it through each step

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

#any non-numeric/placeholder values (e.g. "'--") become NaN to not crash the script
clinical_clean["demographic.days_to_death"] = pd.to_numeric(
    clinical_clean["demographic.days_to_death"],
    errors="coerce"
)

clinical_clean["diagnoses.days_to_last_follow_up"] = pd.to_numeric(
    clinical_clean["diagnoses.days_to_last_follow_up"],
    errors="coerce"
)

#build standard survival-analysis "duration" and "event" columns
#duration: how long the patient was followed. 
#if they died, use days_to_death
#otherwise (still alive / lost to follow-up at last contact), use days_to_last_follow_up
#event: 1 if the patient died (an observed event), 0 if they were censored (alive/unknown at last contact) 
#this is the format lifelines KaplanMeierFitter and logrank_test expect.
clinical_clean["duration"] = clinical_clean["demographic.days_to_death"].fillna(
    clinical_clean["diagnoses.days_to_last_follow_up"]
)

clinical_clean["event"] = (
    clinical_clean["demographic.vital_status"] == "Dead"
).astype(int)

#patients with neither a death time nor a follow-up time can't be used in survival analysis at all, so they're dropped
clinical_clean = clinical_clean.dropna(subset=["duration"])


#load E2F2 scores and convert sample-level IDs to patient-level IDs
#TARGET barcodes look like "TARGET-30-XXXXXX-01A-01R":
#first 16 characters identify the PATIENT and everything after that (like sample type, aliquot, batch) varies even for the same patient
#truncating to 16 characters lets these scores be merged against the clinical table's bare patient IDs
e2f2 = pd.read_csv(args.e2f2_scores_file)

e2f2.columns = ["sample_id", "E2F2_score"]

# TARGET-30-XXXXXX-01A-01R → TARGET-30-XXXXXX
e2f2["cases.submitter_id"] = e2f2["sample_id"].str[:16]

#inner join: keep only patients who have BOTH usable survival data + an E2F2 score
merged = clinical_clean.merge(
    e2f2,
    on="cases.submitter_id",
    how="inner"
)


#optional: attach MYCN amplification status, if a MYCN file was given
#GDC's data export denormalises test records into numbered flat columns ("follow_ups.0.molecular_tests.0.gene_symbol")
#this code checks the first two follow-up "slots" for a MYCN test result
#NOTE: if a patient's MYCN test happens to be recorded under a LATER follow-up index (slot 2, 3, ...),
#this function won't find it and will silently treat that patient as having no MYCN status
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
    
    #left join: keep every patient from `merged` regardless of whether MYCN data is available for them (missing = NaN)
    merged = merged.merge(
        mycn_clean,
        left_on="cases.submitter_id",
        right_on="submitter_id",
        how="left"
    )

    amplified_n = int((merged["mycn_status"] == "Amplified").sum())
    not_amplified_n = int((merged["mycn_status"] == "Not Amplified").sum())

    #don't trust a MYCN-stratified analysis if either subgroup is too small to give a statistically meaningful result
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

#stratify patients into E2F2 high/low groups by median split
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
