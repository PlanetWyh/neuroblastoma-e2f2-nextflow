import scanpy as sc
import numpy as np
import argparse

#this script takes fully processed clustered dataset from preprocess_nb.py
#end pulls out only the cells belonging to NE (identified manually as Leiden cluster "9"
#it then reclusters just the NE popolation 

parser = argparse.ArgumentParser()

parser.add_argument("--input_file", required=True)
parser.add_argument("--output_file", required=True)

args = parser.parse_args()

adata = sc.read_h5ad(args.input_file)

#recomputing same things as in preprocess_nb.py
sc.tl.pca(adata, random_state=777)
sc.pp.neighbors(adata, random_state=777)
sc.tl.umap(adata, random_state=777)
sc.tl.leiden(adata, resolution=0.5, random_state=777)

adata.obs_names_make_unique()

#cluster "9" was identified as the neuroendocrine population based
# on marker gene expression (PCLAF, RGS4, PAGE2, BIRC5)
#suspected NE cluster
ne_clusters = ["9"]
#save from which cluster cells came from ( previous clustering step)
ne = adata[adata.obs["leiden"].isin(ne_clusters)].copy()
ne.obs["full_leiden"] = ne.obs["leiden"].copy()

#reclustering within NE cells only: mainly drop genes that have no significant variability inside NE cluster
sc.pp.filter_genes(ne, min_cells=3)
#rescale again 
sc.pp.scale(ne, max_value=10)
#PCA again to be specific within NE population, with solver suitable to small number of cells
sc.tl.pca(ne, svd_solver="arpack", random_state=777)
sc.pp.neighbors(ne, random_state=777)
sc.tl.umap(ne, random_state=777)
#smaller, heterogeneous population requires smaller resolution
sc.tl.leiden(ne, resolution=0.4, random_state=777)

ne.obs["ne_leiden"] = ne.obs["leiden"].copy()

sc.pl.umap(
    ne,
    color=["full_leiden", "ne_leiden", "PCLAF"],
    cmap="viridis",
    show=False
)

#save the NE object only with original clusters + new reclustering labels
ne.write(args.output_file)
