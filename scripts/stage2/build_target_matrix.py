import argparse
import json
import os

import numpy as np
import pandas as pd

#the TARGET neuroblastoma cohort is downloaded from the GDC as one gene-quantification .tsv file per sequencing file
#each sitting in its own subfolder named after a GDC file_id (NOT a patient ID)
#separate metadata JSON file (downloaded alongside the data, e.g. "metadata.cart.<date>.json") maps each file_id to the actual patient/case ID
#this script is to: read that mapping, read every sample's gene count file, assemble everything into one gene-by-patient expression matrix,
#handle patients who have more than one file (keeping the best one) + normalize to log2(CPM+1) for downstream analysis.
 
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
#build a file_id -> case_id lookup from the GDC metadata JSON
#metadata format is a list of entries, each describing one downloaded file 
#"associated_entities" links that file back to case it came from
#only need the first associated entity here since each TARGET RNA-seq file corresponds to exactly one case
for entry in metadata:
    file_id = entry["file_id"]
    case_id = entry["associated_entities"][0]["entity_submitter_id"]
    file_to_case[file_id] = case_id


#build matrix: wallk through every downloaded file_id subfolder and read its gene-quantification .tsv file.
matrices = []

for file_id in os.listdir(args.data_dir):
    file_path = os.path.join(args.data_dir, file_id)
    #skip anything that isn't a per-file subfolder (e.g. stray files like the metadata)

    if not os.path.isdir(file_path):
        continue

    tsv_files = [f for f in os.listdir(file_path) if f.endswith(".tsv")]

    if not tsv_files:
        continue

    tsv_path = os.path.join(file_path, tsv_files[0])
    
    #skiprows=1 skips a leading comment/metadata line that GDC STAR gene-counts files include above the real header row
    df = pd.read_csv(tsv_path, sep="\t", skiprows=1, header=0)

    #keeping only rows whose gene_id starts with "ENSG" (skips GSC summary rows N_mapped, etc.)
    df = df[df["gene_id"].str.startswith("ENSG")]

    #look up which patient this file belongs to
    #if a file_id is somehow missing from the metadata fall back to using the raw file_id as a stand-in case_id rather than crashing
    case_id = file_to_case.get(file_id, file_id)
    
    #leep only the "unstranded" column which is a standard, non-strand specific raw read count per gene + label this sample's column with its patient case_id
    matrices.append(
        df.set_index("gene_name")["unstranded"].rename(case_id)
    )


if not matrices:
    raise ValueError(f"No expression .tsv files found in {args.data_dir}")

#combine every patient's gene-count Series into one gene x patient matrix (each Series becomes one column, aligned by gene name)
expression_matrix = pd.concat(matrices, axis=1)

print(f"Raw matrix shape: {expression_matrix.shape}")


#deduplicate patients with more than one sequencing file (repeat runs
#for each case_id, keep only the file with the highest total read count
expression_matrix = expression_matrix.T
expression_matrix.index.name = "case_id"

expression_matrix = (
    expression_matrix
    .assign(total=expression_matrix.sum(axis=1)) #rows = patients, cols = genes
    .sort_values("total", ascending=False) #highest-coverage file first
    .groupby(level=0) #group by case_id (the index)
    .first() #since sorted, "first" = highest total
    .drop(columns="total") #drop the helper column
)

expression_matrix = expression_matrix.T  #back to genes x patients

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
