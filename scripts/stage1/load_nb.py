import argparse
import glob
import os

import anndata as ad
import scanpy as sc


parser = argparse.ArgumentParser()

parser.add_argument("--input_dir", required=True)

parser.add_argument("--output_file", required=True)

args = parser.parse_args()


#find all filtered 10x h5 files recursively
pattern = os.path.join(args.input_dir, "**", "*filtered_feature_bc_matrix.h5")
files = sorted(glob.glob(pattern, recursive=True))

if len(files) == 0:
    raise FileNotFoundError(
        f"no files matching found in {args.input_dir}"
    )

adatas = []

for f in files:
    print(f"Loading {f}")

    sample_name = os.path.basename(f).split("_filtered")[0]

    adata = sc.read_10x_h5(f)
    adata.var_names_make_unique()
    adata.obs["sample"] = sample_name

    adatas.append(adata)


#merge all samples
adata = ad.concat(
    adatas,
    join="outer",
    label="batch",
    keys=[a.obs["sample"].iloc[0] for a in adatas],
    index_unique="-"
)

print(adata)


#make sure output folder exists
output_parent = os.path.dirname(args.output_file)
if output_parent:
    os.makedirs(output_parent, exist_ok=True)


#save merged object
adata.write(args.output_file)

print(f"saved merged AnnData to {args.output_file}")