import argparse
import itertools
import os
import tempfile
import typing
import zipfile
import anndata
import matplotlib
import matplotlib.colors
import matplotlib.patches
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import tqdm
import step00

metric_list = ["accuracy", "continuity"]


def exclude_label(exclude: str) -> str:
    return "All chromosomes" if (exclude in step00.missing_values) else f"Without {exclude}"


def cell_clone_labels(clones: pandas.Series, min_clone_size: int = 2) -> typing.List[str]:
    labels = clones.astype(str)
    labels = labels.where(~labels.str.contains("nan", regex=False), "nan")
    counts = labels.value_counts()
    labels = labels.where((labels == "nan") | (labels.map(counts) >= min_clone_size), "nan")
    return labels.tolist()


def mantel_spearman_pvalue(matrix_a: numpy.ndarray, matrix_b: numpy.ndarray, permutations: int, seed: int) -> float:
    upper = numpy.triu_indices_from(matrix_a, k=1)
    rank_a = numpy.zeros(matrix_a.shape, dtype=float)
    rank_a[upper] = scipy.stats.rankdata(matrix_a[upper])
    rank_a = rank_a + rank_a.T
    rank_b = scipy.stats.rankdata(matrix_b[upper])
    observed = abs(numpy.corrcoef(rank_a[upper], rank_b)[0, 1])

    generator = numpy.random.default_rng(seed)
    count = 0
    for _ in tqdm.trange(permutations, postfix=step00.main_file()):
        permutation = generator.permutation(matrix_a.shape[0])
        if abs(numpy.corrcoef(rank_a[permutation][:, permutation][upper], rank_b)[0, 1]) >= observed:
            count += 1
    return (count + 1) / (permutations + 1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input MethylTree H5AD files from run_methyltree.py", type=str, nargs="+")
    parser.add_argument("output", help="Output ZIP file", type=str)
    parser.add_argument("--permutations", help="Number of Mantel permutations", type=int, default=10 ** 4)

    args = parser.parse_args()

    step00.check_suffixes(args.input, {".h5ad", ".h5", ".hdf5"})
    step00.check_suffix(args.output, {".zip"})
    assert (args.permutations > 0), "Permutations must be positive!!"

    adatas: typing.List[anndata.AnnData] = list()
    clone_tables: typing.List[pandas.DataFrame] = list()

    for file in tqdm.tqdm(args.input, postfix=step00.main_file()):
        for adata in tqdm.tqdm(step00.split_methyltree(anndata.read_h5ad(file)), leave=False, postfix=step00.main_file()):
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
    print("Individual:", individual_list)
    print("Exclude:", exclude_list)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    with tempfile.TemporaryDirectory() as directory:
        figure_list: typing.List[str] = list()

        for adata in tqdm.tqdm(adatas, postfix=step00.main_file()):
            parameters = adata.uns["parameters"]
            individual = parameters["Individual"]
            accuracy = adata.uns["accuracy"]

            order = numpy.argsort(adata.obs["MethylTree_order"].to_numpy())
            cell_names: typing.List[str] = adata.obs_names[order].tolist()
            clone_labels: typing.List[str] = cell_clone_labels(adata.obs[parameters["clone"]].iloc[order])
            assert (clone_labels == cell_clone_labels(pandas.Series(adata.uns["ordered_clone_array"]))), f"{individual}: clone labels differ from MethylTree!!"
            cell_clone: typing.Dict[str, str] = dict(zip(cell_names, clone_labels))
            clone_order = sorted(dict.fromkeys(x for x in clone_labels if (x != "nan")))
            codes = numpy.array([(0 if (x in step00.missing_values) else clone_order.index(x) + 1) for x in clone_labels])
            print("Clone order:", len(clone_order), clone_order)

            similarity = pandas.DataFrame(adata.obsm["X_similarity"][order][:, order], index=order, columns=order, dtype=numpy.float64)
            print(similarity)

            local_palette = dict(zip(sorted(set(clone_order)), itertools.cycle(matplotlib.colormaps["tab20"].colors)))
            clone_cmap = matplotlib.colors.ListedColormap([step00.missing_color] + [local_palette[x] for x in clone_order])

            width_ratios, height_ratios = [6, 1, 24, 1], [1, 24]
            fig, ax = matplotlib.pyplot.subplot_mosaic([["Legend", ".", "Horizontal_Color", "."], ["Legend", "Vertical_Color", "Heatmap", "Colorbar"]], figsize=(35, 25), width_ratios=width_ratios, height_ratios=height_ratios)
            seaborn.heatmap(data=similarity, vmin=numpy.percentile(similarity.to_numpy().flatten(), 40), vmax=numpy.percentile(numpy.triu(similarity.to_numpy(), k=1), 99), cmap="viridis", cbar=True, cbar_ax=ax["Colorbar"], xticklabels=False, yticklabels=False, rasterized=True, ax=ax["Heatmap"])
            ax["Heatmap"].set_xlabel(f"Cells (n={len(order)})")
            ax["Vertical_Color"].set_ylabel(f"Cells (n={len(order)})")
            ax["Vertical_Color"].imshow(codes[:, numpy.newaxis], cmap=clone_cmap, vmin=0, vmax=len(clone_order), aspect="auto", interpolation="none", rasterized=True)
            ax["Horizontal_Color"].imshow(codes[numpy.newaxis, :], cmap=clone_cmap, vmin=0, vmax=len(clone_order), aspect="auto", interpolation="none", rasterized=True)
            for name in tqdm.tqdm(["Vertical_Color", "Horizontal_Color"], postfix=step00.main_file()):
                ax[name].grid(False)
                ax[name].set_xticks([])
                ax[name].set_yticks([])
            clone_counts = pandas.Series(clone_labels).value_counts()
            legend_labels = [f"Unassigned (n={clone_counts.get('nan', 0)})"] + [f"{clone} (n={clone_counts[clone]})" for clone in clone_order]
            legend_handles = [matplotlib.patches.Patch(facecolor=clone_cmap(code), edgecolor="none", label=label) for code, label in enumerate(legend_labels)]
            ax["Legend"].set_axis_off()
            ax["Legend"].legend(handles=legend_handles, title="Clone", loc="upper left", bbox_to_anchor=(0, 1), borderaxespad=0, frameon=False)
            ax["Colorbar"].tick_params(axis="y", labelsize="xx-small")
            ax["Colorbar"].set_ylabel("Corrected similarity", fontsize="x-small")
            fig.suptitle(f"{individual} / {exclude_label(parameters['exclude'])}\ncell n={parameters['n_cells']} / region n={parameters['n_regions']}")
            fig.tight_layout()
            step00.save_figure(fig=fig, name=f"{directory}/{individual}-{step00.safe_name(exclude_label(parameters['exclude']))}-Similarity-Heatmap", figure_list=figure_list)

            shared_sites = pandas.DataFrame(adata.obsm["X_shared_sites"][order][:, order], index=order, columns=order, dtype=int)
            upper = numpy.triu_indices_from(shared_sites.to_numpy(), k=1)

            fig, ax = matplotlib.pyplot.subplot_mosaic([["Legend", ".", "Horizontal_Color", "."], ["Legend", "Vertical_Color", "Heatmap", "Colorbar"]], figsize=(35, 25), width_ratios=width_ratios, height_ratios=height_ratios)
            seaborn.heatmap(data=shared_sites, vmin=0, vmax=numpy.percentile(shared_sites.to_numpy()[upper], 99), cmap="Blues", cbar=True, cbar_ax=ax["Colorbar"], xticklabels=False, yticklabels=False, rasterized=True, ax=ax["Heatmap"])
            ax["Heatmap"].set_xlabel(f"Cells (n={len(order)})")
            ax["Vertical_Color"].set_ylabel(f"Cells (n={len(order)})")
            ax["Vertical_Color"].imshow(codes[:, numpy.newaxis], cmap=clone_cmap, vmin=0, vmax=len(clone_order), aspect="auto", interpolation="none", rasterized=True)
            ax["Horizontal_Color"].imshow(codes[numpy.newaxis, :], cmap=clone_cmap, vmin=0, vmax=len(clone_order), aspect="auto", interpolation="none", rasterized=True)
            for name in tqdm.tqdm(["Vertical_Color", "Horizontal_Color"], postfix=step00.main_file()):
                ax[name].grid(False)
                ax[name].set_xticks([])
                ax[name].set_yticks([])
            clone_counts = pandas.Series(clone_labels).value_counts()
            legend_labels = [f"Unassigned (n={clone_counts.get('nan', 0)})"] + [f"{clone} (n={clone_counts[clone]})" for clone in clone_order]
            legend_handles = [matplotlib.patches.Patch(facecolor=clone_cmap(code), edgecolor="none", label=label) for code, label in enumerate(legend_labels)]
            ax["Legend"].set_axis_off()
            ax["Legend"].legend(handles=legend_handles, title="Clone", loc="upper left", bbox_to_anchor=(0, 1), borderaxespad=0, frameon=False)
            ax["Colorbar"].tick_params(axis="y", labelsize="xx-small")
            ax["Colorbar"].set_ylabel("Corrected similarity", fontsize="x-small")
            fig.suptitle(f"{individual} / {exclude_label(parameters['exclude'])}\ncell n={parameters['n_cells']} / region n={parameters['n_regions']}")
            fig.tight_layout()
            step00.save_figure(fig=fig, name=f"{directory}/{individual}-{step00.safe_name(exclude_label(parameters['exclude']))}-SharedSites-Heatmap", figure_list=figure_list)

            pair_data = pandas.DataFrame({"Shared regions": shared_sites.to_numpy()[upper], "Corrected similarity": similarity.to_numpy()[upper], "Pair": numpy.where((codes[upper[0]] > 0) & (codes[upper[0]] == codes[upper[1]]), "Same clone", "Different clone")})
            correlation, p_value = scipy.stats.spearmanr(pair_data["Shared regions"], pair_data["Corrected similarity"])
            mantel_p_value = mantel_spearman_pvalue(shared_sites.to_numpy(), similarity.to_numpy(), args.permutations, 42)
            print(individual, "pairs:", pair_data["Pair"].value_counts().to_dict(), "/ Spearman r:", correlation, "/ p-value", step00.pvalue_format(p_value=p_value), "/ Mantel p:", step00.pvalue_format(p_value=mantel_p_value))

            fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))
            seaborn.regplot(data=pair_data, x="Shared regions", y="Corrected similarity", color=palette[individual], fit_reg=True, robust=True, scatter_kws=dict(s=40, linewidth=0, rasterized=True), ax=ax)
            matplotlib.pyplot.text(numpy.median(pair_data["Shared regions"]), numpy.median(pair_data["Corrected similarity"]), f"Spearman r={correlation:.3f}, p={step00.pvalue_format(p_value=p_value)}", color="black", path_effects=step00.path_effects, fontsize="x-small", horizontalalignment="center", verticalalignment="center")
            matplotlib.pyplot.tight_layout()
            step00.save_figure(fig=fig, name=f"{directory}/{individual}-{step00.safe_name(exclude_label(parameters['exclude']))}-SharedSites-Similarity-Scatter", figure_list=figure_list)

            g = seaborn.lmplot(data=pair_data, x="Shared regions", y="Corrected similarity", hue="Pair", hue_order=["Different clone", "Same clone"], palette={"Different clone": step00.missing_color, "Same clone": palette[individual]}, height=18, aspect=1, legend=True, legend_out=False, scatter=True, fit_reg=True, robust=True, scatter_kws=dict(s=5))
            g.ax.text(0.5, 0.5, f"Spearman r={correlation:.3f}, p={step00.pvalue_format(p_value=p_value)}", color="black", path_effects=step00.path_effects, fontsize="x-small", horizontalalignment="center", verticalalignment="center")
            step00.save_figure(fig=g.figure, name=f"{directory}/{individual}-{step00.safe_name(exclude_label(parameters['exclude']))}-SharedSites-Similarity-Regression", figure_list=figure_list)

            """
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
            seaborn.move_legend(axes.flat[0], loc="lower left", fontsize="x-small", title_fontsize="x-small")
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
            """

        with zipfile.ZipFile(args.output, mode="w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zip_file:
            for figure in tqdm.tqdm(figure_list):
                zip_file.write(figure, arcname=os.path.basename(figure))
