import pandas as pd
import numpy as np
import gseapy as gp
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--logcpm_file", required=True)
parser.add_argument("--geneset_file", required=True)
parser.add_argument("--output_file", required=True)

args = parser.parse_args()

#inputs
logcpm_file = args.logcpm_file
geneset_file = args.geneset_file
output_file = args.output_file

#load logCPM matrix
log_cpm = pd.read_csv(logcpm_file, index_col=0)

#load E2F2 gene set
with open(geneset_file, "r") as f:
    lines = f.read().strip().split("\n")

e2f2_geneset = [g for g in lines if not g.startswith("#")]

print(f"E2F2 targets before cleaning: {len(e2f2_geneset)}")

#clean gene set
e2f2_geneset_clean = [
    g.upper()
    for g in e2f2_geneset
    if g != "E2F2_TARGET_GENES"
    and not g.startswith("ENSG")
]

#uppercase matrix gene names
log_cpm.index = log_cpm.index.astype(str).str.upper()

#check overlap
overlap = [g for g in e2f2_geneset_clean if g in log_cpm.index]
print(f"Clean gene set size: {len(e2f2_geneset_clean)}")
print(f"Overlap with expression matrix: {len(overlap)}")

if len(overlap) < 10:
    raise ValueError("Too few E2F2 target genes overlap with expression matrix.")

#run ssGSEA
ssgsea_results = gp.ssgsea(
    data=log_cpm,
    gene_sets={"E2F2_TARGET_GENES": e2f2_geneset_clean},
    outdir=None,
    sample_norm_method="rank",
    no_plot=True,
    min_size=10,
    max_size=2000
)

#save scores
e2f2_scores = ssgsea_results.res2d.set_index("Name")["NES"].astype(float)
e2f2_scores.name = "E2F2_score"

print(e2f2_scores.describe())

e2f2_scores.to_csv(output_file, header=True)

print(f"Saved {output_file}")