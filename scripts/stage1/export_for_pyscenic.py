import scanpy as sc
import pandas as pd
import argparse

import os

#this script is a format converter
#pySCENIC command line tools don't read AnnData/.h5ad files directly
#they expect a plain CSV expression matrix (cells as rows, genes as
#columns) plus a separate metadata file
#this script converts the
#cleaned NE-only AnnData object into that format

parser = argparse.ArgumentParser()

#ne_clean.h5ad
parser.add_argument("--input_file", required=True)
parser.add_argument("--output_dir", required=True)

args = parser.parse_args()

#adata.X at this point in the pipeline is NOT raw or simply log-normalized expression
#this data has already been: restricted to the top 2,000 highly-variable genes + z-score scaled -- twice 
# (once across the whole tumor + within just the NE cells)
#this is not an implementation recommended by pySCENIC documentation (Z-scoring introduces negative values
#and removes absolute expression magnitude, and restricting to
#HVGs limits which genes can even be considered as TF targets)

adata = sc.read_h5ad(args.input_file)

#pySCENIC wants cells x genes
expr = adata.to_df()


os.makedirs(args.output_dir, exist_ok=True)

csv_file = os.path.join(args.output_dir, "ne_expr_matrix.csv")
metadata_file = os.path.join(args.output_dir, "ne_metadata.csv")

#main expression matrix that GRNBoost2/ctx/AUCell will all read
expr.to_csv(csv_file)

#save metadata for later
#not used by pySCENIC itself, but needed afterwards in pyscenic_plots.py to color the UMAP by cluster
#identity once regulon activity scores come back)

adata.obs[["full_leiden", "ne_leiden"]].to_csv(metadata_file)

# print(expr.shape)
