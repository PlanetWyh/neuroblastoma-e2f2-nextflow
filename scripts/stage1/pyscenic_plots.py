import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import argparse

#this step produces key figure for stage 1 results: a UMAP of the NE cell population 
#colored by E2F2 regulon activity (from pySCENIC's AUCell step) and PCLAF expression 
#(the proliferation marker used to flag the aggressive subpopulation)
#visually overlapping high E2F2 activity with high PCLAF on the
#same embedding is the evidence that E2F2 is specifically active in the most malignant NE subpopulation

parser = argparse.ArgumentParser()

parser.add_argument("--input_matrix", required=True)
parser.add_argument("--pyscenic_resources", required=True)
parser.add_argument("--input_clean", required=True)
parser.add_argument("--auc", required=True)
parser.add_argument("--output_file", required=True)

args = parser.parse_args()

#"ne_expr_matrix.csv"
# expr = pd.read_csv(args.input_matrix, index_col=0)
#pyscenic_resources/allTFs_hg38.txt
# tfs = pd.read_csv(args.pyscenic_resources, header=None)[0].astype(str)
#"ne_clean.h5ad"
adata = sc.read_h5ad(args.input_clean)
#pyscenic_output/auc_mtx.csv
auc = pd.read_csv(args.auc, index_col=0)

#align AUCell scores to the NE object's cell order
#IMPORTANT: this assumes every single cell in `adata` has a matching row in `auc` 
#if pySCENIC internal steps ever drop a cell this line will raise a KeyError rather than silently misaligning anything
auc = auc.loc[adata.obs_names]

#take the E2F2 regulon's per-cell activity score and attach it to the AnnData object as a new column
#so scanpy's plotting functions can use it
adata.obs["E2F2_regulon"] = auc["E2F2(+)"]

#plot two panels side by side on the same NE UMAP
sc.pl.umap(
    adata,
    color=["E2F2_regulon", "PCLAF"],
    cmap="viridis",
    show=False
)

#figures/umap_ne_labelled.png
plt.savefig(args.output_file, dpi=300, bbox_inches="tight")
