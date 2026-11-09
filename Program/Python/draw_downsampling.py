import argparse
import itertools
import math
import os
import tempfile
import typing
import zipfile
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import tqdm
import step00

metric_list = ["accuracy", "continuity"]
x_columns = {"median_shared_regions": "Median shared regions per cell pair", "keep_percent": "Covered regions kept per cell (%)"}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV(.gz) files from downsample_methyltree.py", type=str, nargs="+")
    parser.add_argument("output", help="Output ZIP file", type=str)

    args = parser.parse_args()

    step00.check_suffixes(args.input, {".tsv", ".tsv.gz"})
    step00.check_suffix(args.output, {".zip"})

    input_data = pandas.concat([pandas.read_csv(file, sep="\t") for file in tqdm.tqdm(args.input)], ignore_index=True)
    ok_data = input_data.loc[input_data["status"] == "ok"]

    individual_list = sorted(set(input_data[step00.sample_column]))
    keep_list = sorted(set(input_data["keep_percent"]))
    palette = dict(zip(individual_list, itertools.cycle(matplotlib.colors.TABLEAU_COLORS)))
    failed_counts = input_data.loc[input_data["status"] != "ok"].groupby(step00.sample_column).size()

    summary_data = input_data.groupby([step00.sample_column, "keep_percent"]).agg(runs=("status", "size"), ok_runs=("status", lambda x: int((x == "ok").sum())), median_shared_regions=("median_shared_regions", "median"), min_shared_regions=("min_shared_regions", "min"), accuracy=("accuracy", "mean"), null_accuracy_q95=("null_accuracy_q95", "mean"), significant_runs=("p_accuracy", lambda x: int((x < 0.05).sum())))
    with pandas.option_context("display.max_columns", None, "display.width", 250, "display.max_rows", None):
        print(summary_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    ncols = min(2, len(individual_list))
    nrows = math.ceil(len(individual_list) / ncols)

    with tempfile.TemporaryDirectory() as directory:
        figure_list: typing.List[str] = list()

        for (x_column, x_label), metric in tqdm.tqdm(list(itertools.product(x_columns.items(), metric_list))):
            fig, axes = matplotlib.pyplot.subplots(nrows=nrows, ncols=ncols, figsize=(14 * ncols, 12 * nrows), sharey=True, squeeze=False)
            for ax, individual in zip(axes.flat, individual_list):
                drawing_data = ok_data.loc[ok_data[step00.sample_column] == individual]
                line_data = drawing_data.groupby("keep_percent").agg(x=(x_column, "median"), y=(metric, "mean"), q95=(f"null_{metric}_q95", "mean"), null=(f"null_{metric}_mean", "mean")).sort_values("x")

                ax.plot(line_data["x"], line_data["q95"], color="black", linestyle="--", linewidth=2, label="Null 95th percentile")
                ax.plot(line_data["x"], line_data["null"], color="dimgray", linestyle=":", linewidth=2, label="Null mean")
                ax.plot(line_data["x"], line_data["y"], color=palette[individual], linewidth=4, label="Observed mean")
                ax.scatter(drawing_data[x_column], drawing_data[metric], color=palette[individual], s=200, alpha=0.6, edgecolor="white", linewidth=1, zorder=3, label="Replicate")

                ax.set_xscale("log")
                if x_column == "keep_percent":
                    ax.set_xticks(keep_list, [str(keep) for keep in keep_list], fontsize="x-small", rotation=90)
                    ax.minorticks_off()
                ax.set_ylim(-0.05, 1.05)
                ax.set_xlabel(x_label, fontsize="small")
                ax.set_ylabel(f"Mean per-clone {metric}", fontsize="small")
                failed = int(failed_counts.get(individual, 0))
                ax.set_title(f"{individual} ({failed} runs with zero-shared pairs)" if failed else individual, fontsize="small")
                ax.legend(loc="best", fontsize="xx-small")

            for ax in list(axes.flat)[len(individual_list):]:
                ax.set_visible(False)

            matplotlib.pyplot.tight_layout()
            step00.save_figure(fig, f"{directory}/{metric}-{x_column}", figure_list)

        fig, axes = matplotlib.pyplot.subplots(nrows=nrows, ncols=ncols, figsize=(14 * ncols, 12 * nrows), sharey=True, squeeze=False)
        for ax, individual in zip(axes.flat, individual_list):
            drawing_data = ok_data.loc[ok_data[step00.sample_column] == individual]
            ax.scatter(drawing_data["median_shared_regions"], -numpy.log10(drawing_data["p_accuracy"]), color=palette[individual], s=200, alpha=0.6, edgecolor="white", linewidth=1, zorder=3)
            ax.axhline(-numpy.log10(0.05), color="black", linestyle="--", linewidth=2)
            ax.set_xscale("log")
            ax.set_xlabel(x_columns["median_shared_regions"], fontsize="small")
            ax.set_ylabel("-log10(permutation p), accuracy", fontsize="small")
            ax.set_title(individual, fontsize="small")

        for ax in list(axes.flat)[len(individual_list):]:
            ax.set_visible(False)

        matplotlib.pyplot.tight_layout()
        step00.save_figure(fig, f"{directory}/pvalue-median_shared_regions", figure_list)

        with zipfile.ZipFile(args.output, mode="w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zip_file:
            for figure in tqdm.tqdm(figure_list):
                zip_file.write(figure, arcname=os.path.basename(figure))
