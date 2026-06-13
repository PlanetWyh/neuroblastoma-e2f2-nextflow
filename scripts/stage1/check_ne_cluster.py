import os
import argparse
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--input_file", required=True)
parser.add_argument("--output_dir", required=True)
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

adata = sc.read_h5ad(args.input_file)

# Fix possible log1p/base issue
if "log1p" in adata.uns and "base" in adata.uns["log1p"]:
    del adata.uns["log1p"]["base"]

# Marker sets
ne_markers = ["GAL", "DDX1", "CCND1", "STMN2", "MEG3"]
immune_markers = ["NKG7", "GNLY", "CCL5", "LYZ", "S100A8", "CD79A", "IGKC"]
prolif_markers = ["PCLAF", "BIRC5"]

all_markers = ne_markers + immune_markers + prolif_markers
all_markers = [g for g in all_markers if g in adata.var_names]

present_ne = [g for g in ne_markers if g in adata.var_names]
print("NE markers found:", present_ne)

# Add NE score
sc.tl.score_genes(
    adata,
    gene_list=present_ne,
    score_name="NE_score"
)

# UMAP: leiden
sc.pl.umap(
    adata,
    color="leiden",
    legend_loc="on data",
    title="Leiden clusters",
    show=False
)
plt.savefig(
    os.path.join(args.output_dir, "check_01_umap_leiden.png"),
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# UMAP: NE score
sc.pl.umap(
    adata,
    color="NE_score",
    cmap="viridis",
    title="NE marker score",
    show=False
)
plt.savefig(
    os.path.join(args.output_dir, "check_02_umap_ne_score.png"),
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# UMAP: individual NE markers
sc.pl.umap(
    adata,
    color=present_ne,
    cmap="viridis",
    ncols=3,
    show=False
)
plt.savefig(
    os.path.join(args.output_dir, "check_03_umap_ne_markers.png"),
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# Dotplot
sc.pl.dotplot(
    adata,
    var_names=all_markers,
    groupby="leiden",
    standard_scale="var",
    show=False
)
plt.savefig(
    os.path.join(args.output_dir, "check_04_dotplot_markers.png"),
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# Cluster-level summary
summary = (
    adata.obs
    .groupby("leiden", observed=True)["NE_score"]
    .agg(["mean", "median", "count"])
    .sort_values("mean", ascending=False)
)

summary.to_csv(os.path.join(args.output_dir, "check_05_ne_score_by_cluster.csv"))

print(summary)
print("\nBest NE candidate cluster:", summary.index[0])