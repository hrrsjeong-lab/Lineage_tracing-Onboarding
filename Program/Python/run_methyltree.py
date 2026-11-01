import argparse
import random
import re
import tempfile
import typing
import anndata
import numpy
import pandas
import tqdm
import methyltree
import step00


def region_chromosome(region: str) -> str:
    return re.split(r"[:_]", region, maxsplit=1)[0]


def run_individual(adata: anndata.AnnData, individual: str, clone: str, exclude: typing.List[str], method: str, seed: int) -> anndata.AnnData:
    random.seed(seed)
    numpy.random.seed(seed)

    adata = adata[adata.obs[step00.sample_column] == individual].copy()
    assert (adata.n_obs > 0), f"{individual} is not found!!"
    print(adata)

    covered = (~numpy.isnan(numpy.asarray(adata.X))).sum(axis=0)
    adata = adata[:, (covered > 0)].copy()
    print(individual, "regions covered in this individual:", adata.n_vars)

    chromosomes = adata.var_names.map(region_chromosome)
    print(pandas.Series(chromosomes).value_counts())
    if exclude:
        adata = adata[:, ~(chromosomes.isin(exclude))].copy()
    print(individual, "excluded chromosomes:", exclude, "/ remaining regions:", adata.n_vars)

    clone_counts = adata.obs[clone].astype(str).value_counts()
    print(clone_counts)
    print(individual, "clones with >= 2 cells:", (clone_counts.loc[~clone_counts.index.str.contains("nan", regex=False)] >= 2).sum())

    with tempfile.TemporaryDirectory() as directory:
        adata, accuracy = methyltree.analysis.comprehensive_lineage_analysis(directory, directory, individual, clone, adata_orig=adata, compute_similarity=True, similarity_method=method, remove_celltype_signal=False, similarity_normalize=False, similarity_correction=True, perform_coarse_graining=False, heatmap_vmax_percentile=99.5, heatmap_vmin_percentile=40)

        accuracy_data = pandas.read_csv(f"{directory}/accuracy_{clone}_{individual}.csv")
        with open(f"{directory}/lineage_tree_{clone}_{individual}.txt", "r") as f:
            lineage_tree = f.read()

    print(accuracy_data)
    for key, value in accuracy.items():
        print(individual, key, value)

    order = numpy.asarray(adata.uns["order_x"], dtype=int)
    leaf_order = numpy.empty(adata.n_obs, dtype=int)
    leaf_order[order] = numpy.arange(adata.n_obs)
    adata.obs["MethylTree_order"] = leaf_order

    similarity_keys = [f"X_similarity_{method}_raw", "X_similarity", "X_shared_sites"]
    output_data = anndata.AnnData(obs=adata.obs.copy(), obsm={key: numpy.asarray(adata.obsm[key], dtype=float) for key in similarity_keys})
    output_data.uns["accuracy"] = {key: float(value) for key, value in accuracy.items()}
    output_data.uns["accuracy_table"] = accuracy_data
    output_data.uns["lineage_tree"] = lineage_tree
    output_data.uns["ordered_clone_array"] = numpy.asarray(adata.uns["ordered_clone_array"], dtype=str)
    output_data.uns["parameters"] = {"Individual": individual, "clone": clone, "exclude": ",".join(exclude) if exclude else "none", "method": method, "seed": seed, "n_cells": adata.n_obs, "n_regions": adata.n_vars}
    return output_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input H5AD file from Step 01", type=str)
    parser.add_argument("output", help="Output H5AD file", type=str)
    parser.add_argument("--individuals", help="Individuals to analyze (default: all in the input)", type=str, nargs="*", default=list())
    parser.add_argument("--clone", help="Clone column in obs (truth for accuracy)", type=str, default="large_clone_id")
    parser.add_argument("--exclude", help="Chromosomes to exclude", type=str, nargs="*", default=list())
    parser.add_argument("--method", help="Similarity method", choices=["correlation", "correlation_fast"], default="correlation_fast")
    parser.add_argument("--seed", help="Random seed (MethylTree random baseline)", type=int, default=42)

    args = parser.parse_args()

    step00.check_suffix(args.input, {".h5ad"})
    step00.check_suffix(args.output, {".h5ad"})

    input_data = anndata.read_h5ad(args.input)
    assert (step00.sample_column in input_data.obs.columns), f"{step00.sample_column} is not in obs!!"
    assert (args.clone in input_data.obs.columns), f"{args.clone} is not in obs!!"

    available_individuals = sorted(set(input_data.obs[step00.sample_column].astype(str)))
    individual_list = args.individuals if args.individuals else available_individuals
    assert (set(individual_list) <= set(available_individuals)), f"Unknown individuals: {sorted(set(individual_list) - set(available_individuals))}!!"
    print("Individuals:", individual_list)

    results = [run_individual(input_data, individual, args.clone, args.exclude, args.method, args.seed) for individual in tqdm.tqdm(individual_list)]

    obs_data = pandas.concat([result.obs for result in results])
    n_cells = len(obs_data)
    obsp_data: typing.Dict[str, numpy.ndarray] = dict()
    for key in results[0].obsm.keys():
        matrix = numpy.full((n_cells, n_cells), numpy.nan)
        start = 0
        for result in results:
            end = start + result.n_obs
            matrix[start:end, start:end] = result.obsm[key]
            start = end
        obsp_data[key] = matrix

    output_data = anndata.AnnData(obs=obs_data, obsp=obsp_data)
    output_data.uns["MethylTree"] = {individual: dict(result.uns) for individual, result in zip(individual_list, results)}
    output_data.uns["parameters"] = {"individuals": individual_list, "clone": args.clone, "exclude": ",".join(args.exclude) if args.exclude else "none", "method": args.method, "seed": args.seed}
    print(output_data)

    output_data.write(args.output, **step00.anndata_compressions)
