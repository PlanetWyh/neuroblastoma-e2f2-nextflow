import scanpy as sc
import matplotlib.pyplot as plt
import argparse

#this script takes a AnnData obect from load_nb.py and turns it into clean dataset
#Scanpy QC -> normalization -> dimensionality reduction -> clustering
#workflow, with thresholds chosen to match the original Sun et al. 2024 paper being replicated

parser = argparse.ArgumentParser()

parser.add_argument("--input_file", required=True)
parser.add_argument("--output_file", required=True)

args = parser.parse_args()

adata = sc.read_h5ad(args.input_file)

#apply quality control thresholds
#flag mitochondrial genes whose var_names start with "MT-" in human gene annotation 
# high mitochondrial RNA content in a cell indicates a dying/lysed cell

adata.var["mt"] = adata.var_names.str.startswith("MT-")
#compute per cell: total counts, number of genes detected, and percentage of those counts coming from mitochondrial genes
# stored as adata.obs["pct_counts_mt"]
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

#apply thresholds as in paper
#cells failing any condition are dropped
adata = adata[
    (adata.obs.n_genes_by_counts > 300) &
    (adata.obs.n_genes_by_counts < 7500) &
    (adata.obs.total_counts > 500) &
    (adata.obs.total_counts < 10000) &
    (adata.obs.pct_counts_mt < 25)
].copy()

print(adata)

#normalisation of cell count to 10000 (not all cells get same amount of seqs)
sc.pp.normalize_total(adata, target_sum=1e4)

#log transform needed for PCA to work better
sc.pp.log1p(adata)

#identify top 2000 genes with the most cell-to-cell variability
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000
)

adata = adata[:, adata.var.highly_variable]

#z-score each gene + clip extreme values at 10
sc.pp.scale(adata, max_value=10)

#PCA
sc.tl.pca(adata) #reduce to top principal components

#neighbors
sc.pp.neighbors(adata) #build cell-to-cell similarity graph from PCA space

#UMAP - 2D visualisation
sc.tl.umap(adata)

#Leiden clustering (community detection)
sc.tl.leiden(adata, resolution=0.5)

#Plot
sc.pl.umap(
    adata,
    color=["leiden", "sample"],
    show = False
)

#Save
adata.write(args.output_file)
