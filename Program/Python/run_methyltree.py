import argparse
import random
import re
import tempfile
import anndata
import numpy
import pandas
import methyltree
import step00


def region_chromosome(region: str) -> str:
    return re.split(r"[:_]", region, maxsplit=1)[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input H5AD file from Step 01", type=str)
    parser.add_argument("output", help="Output H5AD file", type=str)
    parser.add_argument("--individual", help="Individual to analyze", type=str, required=True)
    parser.add_argument("--clone", help="Clone column in obs (truth for accuracy)", type=str, default="large_clone_id")
    parser.add_argument("--exclude", help="Chromosomes to exclude", type=str, nargs="*", default=list())
    parser.add_argument("--method", help="Similarity method", choices=["correlation", "correlation_fast"], default="correlation_fast")
    parser.add_argument("--seed", help="Random seed (MethylTree random baseline)", type=int, default=42)

    args = parser.parse_args()

    step00.check_suffix(args.input, {".h5ad"})
    step00.check_suffix(args.output, {".h5ad"})

    random.seed(args.seed)
    numpy.random.seed(args.seed)

    adata = anndata.read_h5ad(args.input)
    assert (step00.sample_column in adata.obs.columns), f"{step00.sample_column} is not in obs!!"
    assert (args.clone in adata.obs.columns), f"{args.clone} is not in obs!!"

    adata = adata[adata.obs[step00.sample_column] == args.individual].copy()
    assert (adata.n_obs > 0), f"{args.individual} is not found!!"
    print(adata)

    covered = (~numpy.isnan(numpy.asarray(adata.X))).sum(axis=0)
    adata = adata[:, (covered > 0)].copy()
    print("Regions covered in this individual:", adata.n_vars)

    chromosomes = adata.var_names.map(region_chromosome)
    print(pandas.Series(chromosomes).value_counts())
    if args.exclude:
        adata = adata[:, ~(chromosomes.isin(args.exclude))].copy()
    print("Excluded chromosomes:", args.exclude, "/ remaining regions:", adata.n_vars)

    clone_counts = adata.obs[args.clone].value_counts(dropna=False)
    print(clone_counts)
    print("Clones with >= 2 cells:", (clone_counts.loc[clone_counts.index.notna()] >= 2).sum())

    with tempfile.TemporaryDirectory() as directory:
        adata, accuracy = methyltree.analysis.comprehensive_lineage_analysis(directory, directory, args.individual, args.clone, adata_orig=adata, compute_similarity=True, similarity_method=args.method, remove_celltype_signal=False, similarity_normalize=False, similarity_correction=True, perform_coarse_graining=False, heatmap_vmax_percentile=99.5, heatmap_vmin_percentile=40)

        accuracy_data = pandas.read_csv(f"{directory}/accuracy_{args.clone}_{args.individual}.csv")
        with open(f"{directory}/lineage_tree_{args.clone}_{args.individual}.txt", "r") as f:
            lineage_tree = f.read()

    print(accuracy_data)
    for key, value in accuracy.items():
        print(key, value)

    order = numpy.asarray(adata.uns["order_x"], dtype=int)
    leaf_order = numpy.empty(adata.n_obs, dtype=int)
    leaf_order[order] = numpy.arange(adata.n_obs)
    adata.obs["MethylTree_order"] = leaf_order

    similarity_keys = [f"X_similarity_{args.method}_raw", "X_similarity", "X_shared_sites"]
    output_data = anndata.AnnData(obs=adata.obs.copy(), var=adata.var.copy(), obsm={key: numpy.asarray(adata.obsm[key], dtype=float) for key in similarity_keys})
    output_data.uns["accuracy"] = {key: float(value) for key, value in accuracy.items()}
    output_data.uns["accuracy_table"] = accuracy_data
    output_data.uns["lineage_tree"] = lineage_tree
    output_data.uns["ordered_clone_array"] = numpy.asarray(adata.uns["ordered_clone_array"], dtype=str)
    output_data.uns["parameters"] = {"Individual": args.individual, "clone": args.clone, "exclude": ",".join(args.exclude) if args.exclude else "none", "method": args.method, "seed": args.seed, "n_cells": adata.n_obs, "n_regions": adata.n_vars}
    print(output_data)

    output_data.write(args.output, **step00.anndata_compressions)
