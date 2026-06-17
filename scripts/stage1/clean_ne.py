import scanpy as sc
import argparse

#this script is just a check to print a couple of sanity checks

parser = argparse.ArgumentParser()
parser.add_argument("--input_file", required=True)
parser.add_argument("--output_file", required=True)
args = parser.parse_args()

ne = sc.read_h5ad(args.input_file)

#print the object summary (number of cells/genes, available obs/var/obsm fields)
print(ne)
# sanity-check: how many NE cells came from each ORIGINAL whole-tumor cluster 
# (should mostly/all be cluster "9", since that's what was
# selected in subset_ne.py).
print("Full dataset clusters:", ne.obs["full_leiden"].value_counts())
#how many cells fall into each of the new NE-specific
#subclusters (e.g. to confirm the PCLAF+ subpopulation is present
#and of a reasonable size)
print("NE subclusters:", ne.obs["ne_leiden"].value_counts())

ne.write(args.output_file)
