import scanpy as sc
#explorative step to identify ne pclaf+ cells
ne = sc.read_h5ad("ne_cells.h5ad")

#statistically rank which genes best distinguish each NE subcluster from the others (Wilcoxon rank-sum test)
#visualise the top distinguishing genes per subcluster
#check expression of the four curated NE proliferation marker genes (Table 1 in the report) across subclusters,
#to assign a biological identity to each one
sc.tl.rank_genes_groups(
    ne,
    groupby="leiden",
    method="wilcoxon"
)
#visualize the top 10 ranked marker genes for every subcluster
sc.pl.rank_genes_groups(
    ne,
    n_genes=10,
    sharey=False
)

#dotplot shows, per gene (rows) and per subcluster (columns), both
#the fraction of cells expressing the gene (dot size) and the
#average expression level (dot color)
sc.pl.dotplot(
    ne,
    var_names=["PCLAF", "RGS4", "PAGE2", "BIRC5"],
    groupby="leiden",
    standard_scale="var"
)

ne.write("ne_cells_markers.h5ad")
