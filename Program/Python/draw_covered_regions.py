import argparse
import itertools
import os
import tempfile
import typing
import zipfile
import matplotlib
import matplotlib.colors
import matplotlib.figure
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import tqdm
import step00

cell_test_list = ["n_covered_regions", "log1p_n_covered_regions", "mean_methylation", "pct_covered_regions_chrX", "pct_covered_regions_chrY", "pct_covered_regions_chrM"]
cell_pair_list = [("log1p_n_covered_regions", "mean_methylation"), ("log1p_n_covered_regions", "pct_covered_regions_chrX"), ("log1p_n_covered_regions", "pct_covered_regions_chrY")]
region_test_list = ["n_covered_cells", "pct_covered_cells", "mean_methylation"]


def label_facet(x, color, label):
    ax = matplotlib.pyplot.gca()
    ax.text(0.01, 0.85, f"{label} (n={len(x)})", color="black", fontsize="small", horizontalalignment="left", verticalalignment="top", transform=ax.transAxes, path_effects=step00.path_effects)


def is_small_integer(values: pandas.Series, max_range: int = 200) -> bool:
    values = values.dropna()
    return bool(numpy.allclose(values, numpy.round(values)) and ((values.max() - values.min()) <= max_range))


def draw_facet_histogram(data: pandas.DataFrame, test: str, individual_list: typing.List[str], palette: typing.Dict[str, str], ylabel: str) -> matplotlib.figure.Figure:
    height = max(6.0, 18.0 / len(individual_list))
    g = seaborn.FacetGrid(data=data, row=step00.sample_column, row_order=individual_list, hue=step00.sample_column, hue_order=individual_list, palette=palette, aspect=24.0 / height, height=height, sharey=False)

    if is_small_integer(data[test]):
        g.map(seaborn.histplot, test, stat="proportion", discrete=True, linewidth=1, edgecolor="white")
    else:
        g.map(seaborn.histplot, test, stat="proportion", bins=numpy.histogram_bin_edges(data[test].dropna(), bins=50), linewidth=1, edgecolor="white")

    g.map(label_facet, test)
    g.set_titles("")
    g.set_axis_labels("", "")
    g.figure.supxlabel(test.replace("_", " "))
    g.figure.supylabel(ylabel)
    g.tight_layout()
    return g.figure


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("cell", help="Input cell QC TSV(.gz) file", type=str)
    parser.add_argument("region", help="Input region QC TSV(.gz) file", type=str)
    parser.add_argument("output", help="Output ZIP file", type=str)

    args = parser.parse_args()

    step00.check_suffix(args.cell, {".tsv", ".tsv.gz"})
    step00.check_suffix(args.region, {".tsv", ".tsv.gz"})
    step00.check_suffix(args.output, {".zip"})

    cell_data = pandas.read_csv(args.cell, sep="\t", index_col=0)
    print(cell_data)

    region_data = pandas.read_csv(args.region, sep="\t", index_col=0)
    print(region_data)

    individual_list = sorted(set(cell_data[step00.sample_column]))
    assert (individual_list == sorted(set(region_data[step00.sample_column]))), "Individuals differ between cell and region tables!!"

    palette = dict(zip(individual_list, itertools.cycle(matplotlib.colors.TABLEAU_COLORS)))
    print("Individual:", len(individual_list), palette)

    assert all(column in cell_data.columns for column in cell_test_list), "Missing cell QC columns!!"
    assert all(column in region_data.columns for column in region_test_list), "Missing region QC columns!!"
    print(cell_data.groupby(step00.sample_column)[cell_test_list].describe().T)
    print(region_data.groupby(step00.sample_column)[region_test_list].describe().T)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    with tempfile.TemporaryDirectory() as directory:
        figure_list: typing.List[str] = list()

        for cell_test in tqdm.tqdm(cell_test_list):
            fig, ax = matplotlib.pyplot.subplots(figsize=(24, 18))
            seaborn.histplot(data=cell_data, x=cell_test, hue=step00.sample_column, hue_order=individual_list, palette=palette, stat="count", bins=50, multiple="stack", linewidth=1, edgecolor="white", ax=ax)
            seaborn.move_legend(ax, loc="upper right", fontsize="x-small", title_fontsize="x-small")
            matplotlib.pyplot.xlabel(cell_test.replace("_", " "))
            matplotlib.pyplot.ylabel("Number of cells")
            matplotlib.pyplot.tight_layout()
            step00.save_figure(fig, f"{directory}/Cell-{cell_test}-Hist", figure_list)

            fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))
            seaborn.violinplot(data=cell_data, x=step00.sample_column, order=individual_list, y=cell_test, hue=step00.sample_column, hue_order=individual_list, palette=palette, inner="box", cut=0, linewidth=2, legend=False, ax=ax)
            counts = cell_data[step00.sample_column].value_counts()
            ax.set_xticks(range(len(individual_list)), [f"{individual}\n(n={counts[individual]})" for individual in individual_list], fontsize="x-small")
            matplotlib.pyplot.xlabel("")
            matplotlib.pyplot.ylabel(cell_test.replace("_", " "))
            matplotlib.pyplot.tight_layout()
            step00.save_figure(fig, f"{directory}/Cell-{cell_test}-Violin", figure_list)

            if len(individual_list) < 2:
                continue

            p_data = pandas.DataFrame(0.0, index=individual_list, columns=individual_list)
            for individual_a in individual_list:
                for individual_b in individual_list:
                    if individual_a == individual_b:
                        continue
                    values_a = cell_data.loc[cell_data[step00.sample_column] == individual_a, cell_test].dropna()
                    values_b = cell_data.loc[cell_data[step00.sample_column] == individual_b, cell_test].dropna()
                    p_value = scipy.stats.mannwhitneyu(values_a, values_b)[1]
                    p_data.loc[individual_a, individual_b] = 0.0 if numpy.isnan(p_value) else -numpy.log10(p_value)
            print(cell_test, "-log10(p), Mann-Whitney U, uncorrected:")
            print(p_data)

            fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))
            seaborn.heatmap(data=p_data, cmap="Reds", vmin=0, annot=True, fmt=".1f", cbar=False, square=True, xticklabels=True, yticklabels=True, ax=ax)
            matplotlib.pyplot.title(cell_test.replace("_", " "))
            matplotlib.pyplot.xlabel("")
            matplotlib.pyplot.ylabel("")
            matplotlib.pyplot.xticks(rotation="vertical", fontsize="x-small")
            matplotlib.pyplot.yticks(rotation="horizontal", fontsize="x-small")
            matplotlib.pyplot.tight_layout()
            step00.save_figure(fig, f"{directory}/Cell-{cell_test}-Heatmap", figure_list)

        for cell_test_1, cell_test_2 in tqdm.tqdm(cell_pair_list):
            g = seaborn.FacetGrid(data=cell_data, col=step00.sample_column, col_order=individual_list, hue=step00.sample_column, hue_order=individual_list, palette=palette, col_wrap=min(2, len(individual_list)), height=12, aspect=1)
            g.map(seaborn.scatterplot, cell_test_1, cell_test_2, s=80, edgecolor="white", linewidth=1)
            g.set_titles("{col_name}")
            g.set_axis_labels(cell_test_1.replace("_", " "), cell_test_2.replace("_", " "))
            g.tight_layout()
            step00.save_figure(g.figure, f"{directory}/Cell-{cell_test_1}-{cell_test_2}-Scatter", figure_list)

        for region_test in tqdm.tqdm(region_test_list):
            step00.save_figure(draw_facet_histogram(region_data, region_test, individual_list, palette, "Proportion of regions"), f"{directory}/Region-{region_test}-HistMultiple", figure_list)

            fig, ax = matplotlib.pyplot.subplots(figsize=(24, 18))
            seaborn.ecdfplot(data=region_data, x=region_test, hue=step00.sample_column, hue_order=individual_list, palette=palette, linewidth=4, ax=ax)
            seaborn.move_legend(ax, loc="lower right", fontsize="x-small", title_fontsize="x-small")
            matplotlib.pyplot.xlabel(region_test.replace("_", " "))
            matplotlib.pyplot.ylabel("Cumulative proportion of regions")
            matplotlib.pyplot.tight_layout()
            step00.save_figure(fig, f"{directory}/Region-{region_test}-ECDF", figure_list)

        x_edges = numpy.arange(region_data["n_covered_cells"].min() - 0.5, region_data["n_covered_cells"].max() + 1.5, 1.0)
        y_edges = numpy.histogram_bin_edges(region_data["mean_methylation"].dropna(), bins=50)
        count_max = max(numpy.histogram2d(data["n_covered_cells"], data["mean_methylation"], bins=(x_edges, y_edges))[0].max() for _, data in region_data.dropna(subset=["mean_methylation"]).groupby(step00.sample_column))
        print("2DHist bins:", len(x_edges) - 1, "x", len(y_edges) - 1, "/ color range: 0 -", count_max)

        g = seaborn.FacetGrid(data=region_data, col=step00.sample_column, col_order=individual_list, col_wrap=min(2, len(individual_list)), height=12, aspect=1)
        g.map_dataframe(seaborn.histplot, x="n_covered_cells", y="mean_methylation", bins=(x_edges, y_edges), cmap="Blues", vmin=0, vmax=count_max, cbar=True)
        g.set_titles("{col_name}")
        g.set_axis_labels("n covered cells", "mean methylation")
        g.tight_layout()
        step00.save_figure(g.figure, f"{directory}/Region-n_covered_cells-mean_methylation-2DHist", figure_list)

        chromosome_data = pandas.crosstab(region_data["Chromosome"], region_data[step00.sample_column]).reindex(columns=individual_list)
        chromosome_data = chromosome_data.loc[sorted(chromosome_data.index, key=step00.chromosome_order)]
        print(chromosome_data)

        fig, ax = matplotlib.pyplot.subplots(figsize=(6 + 4 * len(individual_list), 24))
        seaborn.heatmap(data=chromosome_data, cmap="Blues", annot=True, fmt="d", annot_kws={"fontsize": "xx-small"}, cbar=False, xticklabels=True, yticklabels=True, ax=ax)
        matplotlib.pyplot.xlabel("")
        matplotlib.pyplot.ylabel("")
        matplotlib.pyplot.xticks(rotation="vertical", fontsize="x-small")
        matplotlib.pyplot.yticks(rotation="horizontal", fontsize="x-small")
        matplotlib.pyplot.tight_layout()
        step00.save_figure(fig, f"{directory}/Region-Chromosome-Heatmap", figure_list)

        with zipfile.ZipFile(args.output, mode="w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zip_file:
            for figure in tqdm.tqdm(figure_list):
                zip_file.write(figure, arcname=os.path.basename(figure))
