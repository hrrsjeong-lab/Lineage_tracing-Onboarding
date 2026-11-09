import argparse
import random
import tempfile
import typing
import anndata
import matplotlib.pyplot
import numpy
import pandas
import methyltree
import methyltree.metric
import tqdm
import step00


def mask_coverage(matrix: numpy.ndarray, fraction: float, generator: numpy.random.Generator) -> numpy.ndarray:
    masked = matrix.copy()
    observed = ~numpy.isnan(masked)
    masked[observed & (generator.random(size=masked.shape) >= fraction)] = numpy.nan
    return masked


def shared_region_summary(matrix: numpy.ndarray) -> typing.Dict[str, float]:
    observed = (~numpy.isnan(matrix)).astype(dtype=float)
    shared = (observed @ observed.T)[numpy.triu_indices(n=matrix.shape[0], k=1)]
    return {"median_covered_regions": float(numpy.median(a=observed.sum(axis=1))), "median_shared_regions": float(numpy.median(a=shared)), "min_shared_regions": float(shared.min()), "zero_shared_pairs": int((shared == 0).sum())}


def permutation_null(ordered_clone_array: numpy.ndarray, permutations: int, generator: numpy.random.Generator) -> pandas.DataFrame:
    rows = list()
    for _ in tqdm.trange(permutations, postfix=step00.main_file()):
        table = methyltree.metric.lineage_accuracy_from_leaf_order(ordered_clone_array=generator.permutation(x=ordered_clone_array))
        rows.append({"accuracy": table["accuracy"].astype(dtype=float).mean(), "continuity": table["continuity"].astype(dtype=float).mean()})
    return pandas.DataFrame(data=rows)


def run_individual(adata: anndata.AnnData, individual: str, keep: int, clone: str, replicates: int, permutations: int, method: str, seed: int) -> typing.List[typing.Dict[str, typing.Any]]:
    matrix = numpy.asarray(a=adata.X)
    fraction = keep / 100
    rows: typing.List[typing.Dict[str, typing.Any]] = list()

    for replicate in tqdm.tqdm(iterable=range(1 if (keep == 100) else replicates), leave=False, postfix=step00.main_file()):
        replicate_seed = seed + replicate
        random.seed(a=replicate_seed)
        numpy.random.seed(seed=replicate_seed)
        generator = numpy.random.default_rng(seed=replicate_seed)

        masked = matrix if (keep == 100) else mask_coverage(matrix=matrix, fraction=fraction, generator=generator)
        covered = (~numpy.isnan(masked)).sum(axis=0) > 0
        row: typing.Dict[str, typing.Any] = {"Individual": individual, "keep_percent": keep, "replicate": replicate, "seed": replicate_seed, "n_cells": adata.n_obs, "n_regions": int(covered.sum())} | shared_region_summary(matrix=masked[:, covered])
        print(row)

        if row["zero_shared_pairs"] > 0:
            rows.append(row | {"status": "zero_shared_pairs"})
            continue

        masked_data = anndata.AnnData(X=masked[:, covered], obs=adata.obs.copy(), var=adata.var.loc[covered].copy())
        with tempfile.TemporaryDirectory() as directory:
            masked_data, accuracy = methyltree.analysis.comprehensive_lineage_analysis(out_dir=directory, data_path=directory, save_data_des=individual, clone_key=clone, adata_orig=masked_data, compute_similarity=True, similarity_method=method, remove_celltype_signal=False, similarity_normalize=False, similarity_correction=True, perform_coarse_graining=False, heatmap_vmax_percentile=99.5, heatmap_vmin_percentile=40)
            accuracy_data = pandas.read_csv(filepath_or_buffer=f"{directory}/accuracy_{clone}_{individual}.csv")
        matplotlib.pyplot.close(fig="all")

        null_data = permutation_null(ordered_clone_array=numpy.asarray(a=masked_data.uns["ordered_clone_array"]), permutations=permutations, generator=generator)
        row |= {"status": "ok", "n_clones": len(accuracy_data), "methyltree_random_accuracy": accuracy["mean_random_accuracy"]}
        for metric in ["accuracy", "continuity"]:
            observed = float(accuracy[f"mean_{metric}"])
            row |= {metric: observed, f"null_{metric}_mean": null_data[metric].mean(), f"null_{metric}_q95": null_data[metric].quantile(q=0.95), f"p_{metric}": (1 + (null_data[metric] >= observed).sum()) / (1 + len(null_data))}
        print(row)
        rows.append(row)

    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input H5AD file from Step 01", type=str)
    parser.add_argument("output", help="Output TSV(.gz) file", type=str)
    parser.add_argument("--individuals", help="Individuals to analyze (default: all individuals in the input)", type=str, nargs="+", default=None)
    parser.add_argument("--keep", help="Percents of covered regions to keep per cell", type=int, nargs="+", default=[100, 80, 50, 30, 20, 10, 5])
    parser.add_argument("--clone", help="Clone column in obs (truth for accuracy)", type=str, default="large_clone_id")
    parser.add_argument("--replicates", help="Number of random masks (ignored at 100 percent)", type=int, default=5)
    parser.add_argument("--permutations", help="Number of leaf-order permutations for the null", type=int, default=1000)
    parser.add_argument("--method", help="Similarity method", choices=["correlation", "correlation_fast"], default="correlation_fast")
    parser.add_argument("--seed", help="Random seed", type=int, default=42)

    args = parser.parse_args()

    step00.check_suffix(filename=args.input, suffixes={".h5ad"})
    step00.check_suffix(filename=args.output, suffixes={".tsv", ".tsv.gz"})
    assert all((0 < keep <= 100) for keep in args.keep), f"Keep must be in (0, 100]: {args.keep}!!"
    assert (len(set(args.keep)) == len(args.keep)), f"Duplicated keep: {args.keep}!!"
    assert (args.replicates > 0), "Replicates must be positive!!"
    assert (args.permutations > 0), "Permutations must be positive!!"

    input_data = anndata.read_h5ad(filename=args.input)
    assert (step00.sample_column in input_data.obs.columns), f"{step00.sample_column} is not in obs!!"
    assert (args.clone in input_data.obs.columns), f"{args.clone} is not in obs!!"

    available_individuals = sorted(set(input_data.obs[step00.sample_column].astype(dtype=str)))
    individual_list = available_individuals if (args.individuals is None) else args.individuals
    missing = sorted(set(individual_list) - set(available_individuals))
    assert (not missing), f"{missing} are not found in {available_individuals}!!"
    keep_list = sorted(args.keep, reverse=True)
    print("Individuals:", individual_list)
    print("Keep percents:", keep_list)

    rows: typing.List[typing.Dict[str, typing.Any]] = list()
    for individual in tqdm.tqdm(iterable=individual_list, postfix=step00.main_file()):
        adata = input_data[input_data.obs[step00.sample_column].astype(dtype=str) == individual].copy()
        adata = adata[:, (~numpy.isnan(numpy.asarray(a=adata.X))).sum(axis=0) > 0].copy()
        print(individual, adata)

        for keep in tqdm.tqdm(iterable=keep_list, leave=False, postfix=step00.main_file()):
            rows.extend(run_individual(adata=adata, individual=individual, keep=keep, clone=args.clone, replicates=args.replicates, permutations=args.permutations, method=args.method, seed=args.seed))

    output_data = pandas.DataFrame(data=rows)
    with pandas.option_context("display.max_columns", None, "display.width", 250):
        print(output_data)
    output_data.to_csv(path_or_buf=args.output, sep="\t", index=False)
