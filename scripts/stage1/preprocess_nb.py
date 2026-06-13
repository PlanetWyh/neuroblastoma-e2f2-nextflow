import scanpy as sc
import matplotlib.pyplot as plt
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

#mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    inplace=True
)

#QC plots
# sc.pl.violin(
#     adata,
#     ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
#     multi_panel=True
# )

#apply thresholds as in paper
adata = adata[
    (adata.obs.n_genes_by_counts > 300) &
    (adata.obs.n_genes_by_counts < 7500) &
    (adata.obs.total_counts > 500) &
    (adata.obs.total_counts < 10000) &
    (adata.obs.pct_counts_mt < 25)
].copy()

print(adata)

#normalize
sc.pp.normalize_total(adata, target_sum=1e4)

#log transform
sc.pp.log1p(adata)

#highly variable genes
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000
)

adata = adata[:, adata.var.highly_variable]

#scale
sc.pp.scale(adata, max_value=10)

#PCA
sc.tl.pca(adata)

#neighbors
sc.pp.neighbors(adata)

#UMAP
sc.tl.umap(adata)

#Leiden clustering
sc.tl.leiden(adata, resolution=0.5)

#Plot
sc.pl.umap(
    adata,
    color=["leiden", "sample"],
    show = False
)

#Save
adata.write(args.output_file)