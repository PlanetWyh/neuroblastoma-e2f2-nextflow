import scanpy as sc
import pandas as pd
import argparse

import os

parser = argparse.ArgumentParser()

#ne_clean.h5ad
parser.add_argument("--input_file", required=True)
parser.add_argument("--output_dir", required=True)

args = parser.parse_args()

adata = sc.read_h5ad(args.input_file)

#pySCENIC wants cells x genes
expr = adata.to_df()


os.makedirs(args.output_dir, exist_ok=True)

csv_file = os.path.join(args.output_dir, "ne_expr_matrix.csv")
metadata_file = os.path.join(args.output_dir, "ne_metadata.csv")

expr.to_csv(csv_file)

# save metadata for later

adata.obs[["full_leiden", "ne_leiden"]].to_csv(metadata_file)

# print(expr.shape)