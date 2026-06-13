import scanpy as sc
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input_file", required=True)
parser.add_argument("--output_file", required=True)
args = parser.parse_args()

ne = sc.read_h5ad(args.input_file)

# Do not recluster here.
# Do not filter subclusters here.
# Just save the selected NE population.

print(ne)
print("Full dataset clusters:", ne.obs["full_leiden"].value_counts())
print("NE subclusters:", ne.obs["ne_leiden"].value_counts())

ne.write(args.output_file)