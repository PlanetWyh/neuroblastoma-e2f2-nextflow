import scanpy as sc
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--input_file", required=True)
parser.add_argument("--output_file", required=True)

args = parser.parse_args()

adata = sc.read_h5ad(args.input_file)

sc.tl.pca(adata, random_state=777)
sc.pp.neighbors(adata, random_state=777)
sc.tl.umap(adata, random_state=777)
sc.tl.leiden(adata, resolution=0.5, random_state=777)

adata.obs_names_make_unique()

#suspected NE cluster
ne_clusters = ["9"]

ne = adata[adata.obs["leiden"].isin(ne_clusters)].copy()
ne.obs["full_leiden"] = ne.obs["leiden"].copy()

sc.pp.filter_genes(ne, min_cells=3)
sc.pp.scale(ne, max_value=10)
sc.tl.pca(ne, svd_solver="arpack", random_state=777)
sc.pp.neighbors(ne, random_state=777)
sc.tl.umap(ne, random_state=777)
sc.tl.leiden(ne, resolution=0.4, random_state=777)

ne.obs["ne_leiden"] = ne.obs["leiden"].copy()

sc.pl.umap(
    ne,
    color=["full_leiden", "ne_leiden", "PCLAF"],
    cmap="viridis",
    show=False
)

ne.write(args.output_file)