import pandas as pd
import numpy as np
import gseapy as gp
import argparse

#this script computes for every TARGET patient a single number which summarises how active the E2F2 target gene program is in their tumor
#using single-sample Gene Set Enrichment Analysis (ssGSEA) 
#this score can later be correlated with survival


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
# MSigDB .grp files are one gene per line, "#" lines are comments and get skipped
with open(geneset_file, "r") as f:
    lines = f.read().strip().split("\n")

e2f2_geneset = [g for g in lines if not g.startswith("#")]

print(f"E2F2 targets before cleaning: {len(e2f2_geneset)}")
#clean gene set
#two cleanup rules: drop the literal string "E2F2_TARGET_GENES" if it shows up as an entry
#+ drop any Ensembl gene IDs (starting with "ENSG") since the expression matrix is indexed by gene SYMBOL, not Ensembl ID
e2f2_geneset_clean = [
    g.upper()
    for g in e2f2_geneset
    if g != "E2F2_TARGET_GENES"
    and not g.startswith("ENSG")
]

#uppercase the matrix's gene names too, so matching against the (already uppercased) gene set isn't broken by case differences
log_cpm.index = log_cpm.index.astype(str).str.upper()

#check overlap
#how many of the E2F2 target genes actually exist in this expression matrix? 
#if too few do, the resulting ssGSEA score would be unreliable, so the script deliberately crashes here
overlap = [g for g in e2f2_geneset_clean if g in log_cpm.index]
print(f"Clean gene set size: {len(e2f2_geneset_clean)}")
print(f"Overlap with expression matrix: {len(overlap)}")

if len(overlap) < 10:
    raise ValueError("Too few E2F2 target genes overlap with expression matrix.")

#run ssGSEA
#sample_norm_method="rank": within each patient, gene expression values are converted to ranks before scoring
#this makes scores comparable across patients regardless of absolute expression scale, and is the standard ssGSEA preprocessing choice
#outdir=None / no_plot=True: suppress gseapy's own file/plot output
#min_size / max_size: gene set size bounds gseapy enforces
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
#gseapy returns one row per patient in ssgsea_results.res2d, with several score columns
#NES ("Normalized Enrichment Score") is the standard ssGSEA output to use
#this score is normalised for gene set size making it comparable across analyses
#this is my per-patient "E2F2 activity score" used for the rest of Stage 2
e2f2_scores = ssgsea_results.res2d.set_index("Name")["NES"].astype(float)
e2f2_scores.name = "E2F2_score"

print(e2f2_scores.describe())

e2f2_scores.to_csv(output_file, header=True)

print(f"Saved {output_file}")
