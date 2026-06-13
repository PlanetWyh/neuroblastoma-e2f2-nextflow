import scanpy as sc

ne = sc.read_h5ad("ne_cells.h5ad")

sc.tl.rank_genes_groups(
    ne,
    groupby="leiden",
    method="wilcoxon"
)

sc.pl.rank_genes_groups(
    ne,
    n_genes=10,
    sharey=False
)

sc.pl.dotplot(
    ne,
    var_names=["PCLAF", "RGS4", "PAGE2", "BIRC5"],
    groupby="leiden",
    standard_scale="var"
)

ne.write("ne_cells_markers.h5ad")