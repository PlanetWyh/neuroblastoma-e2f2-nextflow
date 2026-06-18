import os
import gzip
import argparse
import urllib.request

import GEOparse
import pandas as pd
import gseapy as gp

#this script independently reproduces the E2F2 ssGSEA scoring pipeline on a completely separate neuroblastoma cohort (GSE49711) to check that
#the TARGET survival association isn't a cohort-specific artifact
#this script handles everything from scratch: downloading the GEO metadata and expression data, computing ssGSEA scores with the
#same gene signature used for TARGET, and producing a merged clinical+scores table ready for plotting

parser = argparse.ArgumentParser()
parser.add_argument("--gse_id", default="GSE49711")
parser.add_argument("--geneset_file", required=True)
parser.add_argument("--output_dir", required=True)
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

#download and parse GEO metadata using GEOparse: GEOparse fetches the GSE SOFT file 
#(structured flat-file format GEO uses to store dataset and sample metadata) and parses it into python object
#each individual sample holds per-sample clinical annotations under metadata["characteristics_ch1"] 
#thay are stored as a list of "key: value" strings
gse = GEOparse.get_GEO(args.gse_id, destdir=args.output_dir)

clinical_rows = []
gsm_to_title = {}

for gsm_name, gsm in gse.gsms.items():
    row = {"gsm": gsm_name}
    #"title" is a short label GEO depositors assign to each sample (e.g. "SEQC_NB001")
    #need this later to map expression matrix columns (they use these titles) back to GSM accession IDs (clinical data uses this)
    title = gsm.metadata["title"][0]
    gsm_to_title[gsm_name] = title

    #parse the free-text characteristics into a flat key:value dict.
    #GEO characteristics_ch1 entries look like "age at diagnosis: 365", "mycn status: 1"
    #splitting on semicolon gives  column names and values
    chars = gsm.metadata["characteristics_ch1"]
    for c in chars:
        if ": " in c:
            key, val = c.split(": ", 1)
            row[key.strip()] = val.strip()

    clinical_rows.append(row)

clinical = pd.DataFrame(clinical_rows)
clinical_file = os.path.join(args.output_dir, "GSE49711_clinical.csv")
clinical.to_csv(clinical_file, index=False)

print("Clinical:", clinical.shape)

#download the expression matrix directly from GEO 
#no need to re-normalise since GSE49711 provides a pre-processed log2-transformed expression matrix as a supplementary file
#the file is gzip-compressed
url = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/"
    "GSE49711/suppl/GSE49711_SEQC_NB_TUC_G_log2.txt.gz"
)

expr_gz = os.path.join(args.output_dir, "GSE49711_expression.txt.gz")
#only download if the file isn't already on disk
if not os.path.exists(expr_gz):
    print("Downloading GSE49711 expression matrix...")
    urllib.request.urlretrieve(url, expr_gz)

with gzip.open(expr_gz, "rt") as f:
    expr = pd.read_csv(f, sep="\t", index_col=0)

print("Expression:", expr.shape)

#map expression matrix columns to samples ids: the expression matrix columns use short sample titles ("SEQC_NB001") but not samples ids
#first truncate to 10 characters to strip off any suffix variation and then look up the corresponding sample id
# using the title->sample mapping built above
expr.columns = [c[:10] for c in expr.columns]  # SEQC_NB001 

title_to_gsm = {title: gsm for gsm, title in gsm_to_title.items()}
expr.columns = [title_to_gsm.get(c, c) for c in expr.columns]

#keep only samples present in BOTH the expression matrix and clinical table (inner join on sample ID)
clinical = clinical.set_index("gsm")

overlap = sorted(set(expr.columns) & set(clinical.index))
expr = expr[overlap]
clinical = clinical.loc[overlap].copy()

print("Matched samples:", len(overlap))

#uppercase gene names for consistent matching against the gene set (same step as in calculate_e2f2_ssgsea.py)
expr.index = expr.index.astype(str).str.upper()

expr_file = os.path.join(args.output_dir, "GSE49711_expression_log2.csv")
expr.to_csv(expr_file)

#load E2F2 gene set

with open(args.geneset_file, "r") as f:
    genes = f.read().strip().splitlines()

e2f2_genes = [
    g.upper()
    for g in genes
    if g
    and not g.startswith("#")
    and g != "E2F2_TARGET_GENES"
    and not g.startswith("ENSG")
]

overlap_genes = [g for g in e2f2_genes if g in expr.index]

print(f"E2F2 genes in GSE49711: {len(overlap_genes)} / {len(e2f2_genes)}")

if len(overlap_genes) < 10:
    raise ValueError("Too few E2F2 genes overlap with expression matrix.")

#ssGSEA with same parameters as for target cohort
ssgsea = gp.ssgsea(
    data=expr,
    gene_sets={"E2F2_TARGET_GENES": e2f2_genes},
    outdir=None,
    sample_norm_method="rank",
    no_plot=True,
    min_size=10,
    max_size=2000,
)

scores = ssgsea.res2d.set_index("Name")["NES"].astype(float)
scores.name = "E2F2_score"

scores_file = os.path.join(args.output_dir, "GSE49711_e2f2_ssgsea_scores.csv")
scores.to_csv(scores_file, header=True)


#merge clinical + E2F2 score + split into high/low groups by median
clinical["E2F2_score"] = scores
clinical["E2F2_group"] = (
    clinical["E2F2_score"] >= clinical["E2F2_score"].median()
).map({True: "High", False: "Low"})

#convert clinical variable columns to numeric where possible
#GSE49711 stores binary outcomes (death, progression, mycn status, risk group) as "0"/"1" strings
#and age as a number so use pd.to_numeric with errors="coerce" to handle both, convert N/A to NaN and avoid crashing
for col in ["death from disease", "progression", "age at diagnosis", "mycn status", "high risk"]:
    if col in clinical.columns:
        clinical[col] = clinical[col].replace("N/A", pd.NA)
        clinical[col] = pd.to_numeric(clinical[col], errors="coerce")

merged_file = os.path.join(args.output_dir, "GSE49711_e2f2_clinical_merged.csv")
clinical.to_csv(merged_file)

print(clinical["E2F2_group"].value_counts())
print("Saved:", merged_file)
