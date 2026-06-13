import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("--input_matrix", required=True)
parser.add_argument("--pyscenic_resources", required=True)
parser.add_argument("--input_clean", required=True)
parser.add_argument("--auc", required=True)
parser.add_argument("--output_file", required=True)

args = parser.parse_args()

#"ne_expr_matrix.csv"
expr = pd.read_csv(args.input_matrix, index_col=0)
#pyscenic_resources/allTFs_hg38.txt
tfs = pd.read_csv(args.pyscenic_resources, header=None)[0].astype(str)
#"ne_clean.h5ad"
adata = sc.read_h5ad(args.input_clean)
#pyscenic_output/auc_mtx.csv
auc = pd.read_csv(args.auc, index_col=0)

# align cells
auc = auc.loc[adata.obs_names]

# add regulon score
adata.obs["E2F2_regulon"] = auc["E2F2(+)"]

sc.pl.umap(
    adata,
    color=["E2F2_regulon", "PCLAF"],
    cmap="viridis",
    show=False
)

#figures/umap_ne_labelled.png
plt.savefig(args.output_file, dpi=300, bbox_inches="tight")
