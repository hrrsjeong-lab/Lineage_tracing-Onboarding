import argparse
import itertools
import os
import tempfile
import typing
import zipfile
import anndata
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

metric_list = ["accuracy", "continuity"]


def exclude_label(exclude: str) -> str:
    return "All chromosomes" if (exclude.lower() == "none") else f"Without {exclude}"


def clone_codes(clones: pandas.Series) -> typing.Tuple[numpy.ndarray, typing.List[str]]:
    clones = clones.astype(object).where(clones.notna(), None)
    counts = clones.value_counts()
    clones = clones.where(clones.map(lambda x: (x is not None) and (counts.get(x, 0) >= 2)), None)
    clone_order = list(dict.fromkeys(clones.dropna()))
    codes = numpy.array([0 if (x is None) else clone_order.index(x) + 1 for x in clones])
    return codes, clone_order


def draw_ordered_heatmap(matrix: numpy.ndarray, codes: numpy.ndarray, n_clones: int, cmap: str, vmin: float, vmax: float, title: str, colorbar_label: str) -> matplotlib.figure.Figure:
    clone_colors = [step00.missing_color] + [color for _, color in zip(range(n_clones), itertools.cycle(matplotlib.colormaps["tab20"].colors))]
    clone_cmap = matplotlib.colors.ListedColormap(clone_colors)

    fig = matplotlib.pyplot.figure(figsize=(26, 25))
    grid = fig.add_gridspec(2, 3, width_ratios=[1, 24, 1], height_ratios=[1, 24], wspace=0.02, hspace=0.02)
    ax_top = fig.add_subplot(grid[0, 1])
    ax_left = fig.add_subplot(grid[1, 0])
    ax = fig.add_subplot(grid[1, 1])
    cax = fig.add_subplot(grid[1, 2])

    image = ax.imshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto", interpolation="none")
    ax_top.imshow(codes[numpy.newaxis, :], cmap=clone_cmap, vmin=0, vmax=n_clones, aspect="auto", interpolation="none")
    ax_left.imshow(codes[:, numpy.newaxis], cmap=clone_cmap, vmin=0, vmax=n_clones, aspect="auto", interpolation="none")
    fig.colorbar(image, cax=cax, label=colorbar_label)

    for a in tqdm.tqdm([ax, ax_top, ax_left]):
        a.grid(False)
        a.set_xticks([])
        a.set_yticks([])
    ax.set_xlabel("Cells (MethylTree leaf order)")
    ax_left.set_ylabel("Cells (MethylTree leaf order)")
    ax_top.set_title(title, fontsize="small")
    fig.subplots_adjust(left=0.06, right=0.90, top=0.90, bottom=0.05)
    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input MethylTree H5AD files from run_methyltree.py", type=str, nargs="+")
    parser.add_argument("output", help="Output ZIP file", type=str)

    args = parser.parse_args()

    step00.check_suffixes(args.input, {".h5ad"})
    step00.check_suffix(args.output, {".zip"})

    adatas: typing.List[anndata.AnnData] = list()
    clone_tables: typing.List[pandas.DataFrame] = list()

    for file in tqdm.tqdm(args.input):
        adata = anndata.read_h5ad(file)
        adatas.append(adata)

        clone_data = adata.uns["accuracy_table"].copy()
        clone_data[step00.sample_column] = adata.uns["parameters"]["Individual"]
        clone_data["Exclude"] = exclude_label(adata.uns["parameters"]["exclude"])
        clone_tables.append(clone_data)

    clone_data = pandas.concat(clone_tables, ignore_index=True)
    clone_data["clone"] = clone_data["clone"].astype(str)
    print(clone_data)

    individual_list = sorted(set(clone_data[step00.sample_column]))
    exclude_list = sorted(set(clone_data["Exclude"]), key=lambda x: (x != exclude_label("none"), x))
    palette = dict(zip(individual_list, itertools.cycle(matplotlib.colors.TABLEAU_COLORS)))
    print("Individual:", individual_list, "/ Exclude:", exclude_list)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    with tempfile.TemporaryDirectory() as directory:
        figure_list: typing.List[str] = list()

        for adata in tqdm.tqdm(adatas):
            parameters = adata.uns["parameters"]
            individual = parameters["Individual"]
            name = f"{individual}-{exclude_label(parameters['exclude']).replace(' ', '_').replace(',', '_')}"
            accuracy = adata.uns["accuracy"]

            order = numpy.argsort(adata.obs["MethylTree_order"].to_numpy())
            codes, clone_order = clone_codes(adata.obs[parameters["clone"]].iloc[order].reset_index(drop=True))
            title = f"{individual} / {exclude_label(parameters['exclude'])} / {parameters['n_cells']} cells / {parameters['n_regions']} regions\naccuracy {accuracy['mean_accuracy']:.3f} (random {accuracy['mean_random_accuracy']:.3f}) / continuity {accuracy['mean_continuity']:.3f} (random {accuracy['mean_random_continuity']:.3f})"

            similarity = adata.obsm["X_similarity"][order][:, order]
            vmin = numpy.percentile(similarity.flatten(), 40)
            vmax = numpy.percentile(numpy.triu(similarity, k=1), 99.5)
            print(name, "similarity color range:", vmin, vmax)
            step00.save_figure(draw_ordered_heatmap(similarity, codes, len(clone_order), "viridis", vmin, vmax, title, "Corrected similarity"), f"{directory}/{name}-Similarity-Heatmap", figure_list)

            shared_sites = adata.obsm["X_shared_sites"][order][:, order]
            upper = numpy.triu_indices_from(shared_sites, k=1)
            step00.save_figure(draw_ordered_heatmap(shared_sites, codes, len(clone_order), "Blues", 0, numpy.percentile(shared_sites[upper], 99.5), title, "Shared covered regions"), f"{directory}/{name}-SharedSites-Heatmap", figure_list)

            pair_data = pandas.DataFrame({"Shared regions": shared_sites[upper], "Corrected similarity": similarity[upper], "Pair": numpy.where((codes[upper[0]] > 0) & (codes[upper[0]] == codes[upper[1]]), "Same clone", "Different clone")})
            correlation = scipy.stats.spearmanr(pair_data["Shared regions"], pair_data["Corrected similarity"])[0]
            print(name, "pairs:", pair_data["Pair"].value_counts().to_dict(), "/ Spearman r (shared regions, similarity):", correlation)

            fig, ax = matplotlib.pyplot.subplots(figsize=(24, 18))
            seaborn.scatterplot(data=pair_data, x="Shared regions", y="Corrected similarity", hue="Pair", hue_order=["Different clone", "Same clone"], palette={"Different clone": step00.missing_color, "Same clone": palette[individual]}, s=40, linewidth=0, rasterized=True, ax=ax)
            seaborn.move_legend(ax, loc="upper right", fontsize="x-small", title_fontsize="x-small")
            ax.text(0.01, 0.99, f"Spearman r = {correlation:.3f} (all {len(pair_data)} pairs)", fontsize="x-small", horizontalalignment="left", verticalalignment="top", transform=ax.transAxes, path_effects=step00.path_effects)
            matplotlib.pyplot.title(f"{individual} / {exclude_label(parameters['exclude'])}", fontsize="small")
            matplotlib.pyplot.tight_layout()
            step00.save_figure(fig, f"{directory}/{name}-SharedSites-Similarity-Scatter", figure_list)

        for metric in tqdm.tqdm(metric_list):
            observed_data = clone_data[[step00.sample_column, "Exclude", "clone", "clone_size", metric]].assign(Type="Observed")
            random_data = clone_data[[step00.sample_column, "Exclude", "clone", "clone_size", f"random_{metric}"]].rename(columns={f"random_{metric}": metric}).assign(Type="Random")
            long_data = pandas.concat([observed_data, random_data], ignore_index=True)

            fig, axes = matplotlib.pyplot.subplots(ncols=len(exclude_list), figsize=(6 + 4 * len(individual_list) * len(exclude_list), 18), sharey=True, squeeze=False)
            for ax, exclude in zip(axes.flat, exclude_list):
                drawing_data = long_data.loc[(long_data["Exclude"] == exclude)]
                seaborn.boxplot(data=drawing_data, x=step00.sample_column, order=individual_list, y=metric, hue="Type", hue_order=["Observed", "Random"], palette={"Observed": "tab:gray", "Random": "white"}, linewidth=2, fliersize=0, legend=(ax is axes.flat[0]), ax=ax)
                seaborn.stripplot(data=drawing_data, x=step00.sample_column, order=individual_list, y=metric, hue="Type", hue_order=["Observed", "Random"], palette={"Observed": "black", "Random": "darkgray"}, dodge=True, size=8, jitter=0.2, legend=False, ax=ax)
                counts = drawing_data.loc[drawing_data["Type"] == "Observed"].groupby(step00.sample_column).size()
                ax.set_xticks(range(len(individual_list)), [f"{individual}\n({counts.get(individual, 0)} clones)" for individual in individual_list], fontsize="x-small", rotation=90)
                ax.set_title(exclude, fontsize="small")
                ax.set_xlabel("")
                ax.set_ylabel(f"Per-clone {metric}")
                ax.set_ylim(-0.05, 1.05)
            seaborn.move_legend(axes.flat[0], loc="center left", fontsize="x-small", title_fontsize="x-small")
            matplotlib.pyplot.tight_layout()
            step00.save_figure(fig, f"{directory}/Clone-{metric}-Box", figure_list)

            g = seaborn.relplot(data=clone_data, x="clone_size", y=metric, hue=step00.sample_column, hue_order=individual_list, palette=palette, col="Exclude", col_order=exclude_list, s=400, alpha=0.6, edgecolor="white", linewidth=1, height=18, aspect=1.0, facet_kws={"legend_out": False})
            g.set_titles("{col_name}")
            g.set_axis_labels("Clone size (cells)", f"Per-clone {metric}")
            g.set(ylim=(-0.05, 1.05))
            g.tight_layout()
            step00.save_figure(g.figure, f"{directory}/Clone-size-{metric}-Scatter", figure_list)

            if len(exclude_list) < 2:
                continue

            pair_data = clone_data.pivot_table(index=[step00.sample_column, "clone"], columns="Exclude", values=metric).dropna().reset_index()
            print(metric, exclude_list[0], "vs", exclude_list[1], "clones:", len(pair_data))

            g = seaborn.relplot(data=pair_data, x=exclude_list[0], y=exclude_list[1], hue=step00.sample_column, hue_order=individual_list, palette=palette, col=step00.sample_column, col_order=individual_list, col_wrap=min(2, len(individual_list)), s=400, alpha=0.6, edgecolor="white", linewidth=1, height=12, aspect=1.0, legend=False)
            counts = pair_data.groupby(step00.sample_column).size()
            for ax, individual in zip(g.axes.flat, individual_list):
                ax.plot([0, 1], [0, 1], color="black", linestyle="--", linewidth=2, zorder=0)
                ax.set_title(f"{individual} ({counts.get(individual, 0)} clones)")
            g.set_axis_labels(f"{exclude_list[0]}: {metric}", f"{exclude_list[1]}: {metric}")
            g.set(xlim=(-0.05, 1.05), ylim=(-0.05, 1.05))
            g.tight_layout()
            step00.save_figure(g.figure, f"{directory}/Clone-{metric}-Exclude-Scatter", figure_list)

        with zipfile.ZipFile(args.output, mode="w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zip_file:
            for figure in tqdm.tqdm(figure_list):
                zip_file.write(figure, arcname=os.path.basename(figure))
