import argparse
import json
import os

import numpy as np
import pandas as pd


parser = argparse.ArgumentParser()

parser.add_argument("--metadata", required=True)
parser.add_argument("--data_dir", required=True)
parser.add_argument("--counts_out", required=True)
parser.add_argument("--logcpm_out", required=True)

args = parser.parse_args()


#load metadata and map file_id to case_id
with open(args.metadata, "r") as f:
    metadata = json.load(f)

file_to_case = {}

for entry in metadata:
    file_id = entry["file_id"]
    case_id = entry["associated_entities"][0]["entity_submitter_id"]
    file_to_case[file_id] = case_id


#build matrix
matrices = []

for file_id in os.listdir(args.data_dir):
    file_path = os.path.join(args.data_dir, file_id)

    if not os.path.isdir(file_path):
        continue

    tsv_files = [f for f in os.listdir(file_path) if f.endswith(".tsv")]

    if not tsv_files:
        continue

    tsv_path = os.path.join(file_path, tsv_files[0])

    df = pd.read_csv(tsv_path, sep="\t", skiprows=1, header=0)
    df = df[df["gene_id"].str.startswith("ENSG")]

    case_id = file_to_case.get(file_id, file_id)

    matrices.append(
        df.set_index("gene_name")["unstranded"].rename(case_id)
    )


if not matrices:
    raise ValueError(f"No expression .tsv files found in {args.data_dir}")


expression_matrix = pd.concat(matrices, axis=1)

print(f"Raw matrix shape: {expression_matrix.shape}")


#deduplicate samples
expression_matrix = expression_matrix.T
expression_matrix.index.name = "case_id"

expression_matrix = (
    expression_matrix
    .assign(total=expression_matrix.sum(axis=1))
    .sort_values("total", ascending=False)
    .groupby(level=0)
    .first()
    .drop(columns="total")
)

expression_matrix = expression_matrix.T

print(f"Matrix shape after deduplication: {expression_matrix.shape}")

expression_matrix.to_csv(args.counts_out)


#CPM + log2 normalization
cpm = expression_matrix.div(expression_matrix.sum(axis=0), axis=1) * 1e6
log_cpm = cpm.apply(lambda x: np.log2(x + 1))

print(f"logCPM matrix shape: {log_cpm.shape}")

if "MYCN" in log_cpm.index:
    print(log_cpm.loc["MYCN"].describe())

log_cpm.to_csv(args.logcpm_out)

print(f"Saved counts matrix to: {args.counts_out}")
print(f"Saved logCPM matrix to: {args.logcpm_out}")