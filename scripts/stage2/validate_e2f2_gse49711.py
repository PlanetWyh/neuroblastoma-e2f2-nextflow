import os
import gzip
import argparse
import urllib.request

import GEOparse
import pandas as pd
import gseapy as gp


parser = argparse.ArgumentParser()
parser.add_argument("--gse_id", default="GSE49711")
parser.add_argument("--geneset_file", required=True)
parser.add_argument("--output_dir", required=True)
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

#download GEO metadata
gse = GEOparse.get_GEO(args.gse_id, destdir=args.output_dir)

clinical_rows = []
gsm_to_title = {}

for gsm_name, gsm in gse.gsms.items():
    row = {"gsm": gsm_name}

    title = gsm.metadata["title"][0]
    gsm_to_title[gsm_name] = title

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

#download expression matrix
url = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/"
    "GSE49711/suppl/GSE49711_SEQC_NB_TUC_G_log2.txt.gz"
)

expr_gz = os.path.join(args.output_dir, "GSE49711_expression.txt.gz")

if not os.path.exists(expr_gz):
    print("Downloading GSE49711 expression matrix...")
    urllib.request.urlretrieve(url, expr_gz)

with gzip.open(expr_gz, "rt") as f:
    expr = pd.read_csv(f, sep="\t", index_col=0)

print("Expression:", expr.shape)

#map expression columns to GSM IDs
expr.columns = [c[:10] for c in expr.columns]  # SEQC_NB001 etc.

title_to_gsm = {title: gsm for gsm, title in gsm_to_title.items()}
expr.columns = [title_to_gsm.get(c, c) for c in expr.columns]

clinical = clinical.set_index("gsm")

overlap = sorted(set(expr.columns) & set(clinical.index))
expr = expr[overlap]
clinical = clinical.loc[overlap].copy()

print("Matched samples:", len(overlap))

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

#ssGSEA
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


#merge clinical + E2F2 score
clinical["E2F2_score"] = scores
clinical["E2F2_group"] = (
    clinical["E2F2_score"] >= clinical["E2F2_score"].median()
).map({True: "High", False: "Low"})

for col in ["death from disease", "progression", "age at diagnosis", "mycn status", "high risk"]:
    if col in clinical.columns:
        clinical[col] = clinical[col].replace("N/A", pd.NA)
        clinical[col] = pd.to_numeric(clinical[col], errors="coerce")

merged_file = os.path.join(args.output_dir, "GSE49711_e2f2_clinical_merged.csv")
clinical.to_csv(merged_file)

print(clinical["E2F2_group"].value_counts())
print("Saved:", merged_file)