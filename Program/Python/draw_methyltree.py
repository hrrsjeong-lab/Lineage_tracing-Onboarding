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
import statannotations.Annotator
import tqdm
import notebook_pdf
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

    step00.check_suffixes(filenames=args.input, suffixes={".h5ad", ".h5", ".hdf5"})
    step00.check_suffix(filename=args.output, suffixes={".zip"})
    assert (args.permutations > 0), "Permutations must be positive!!"

    notebook = notebook_pdf.Notebook(title="Draw MethylTree", script=__file__, packages=("anndata", "matplotlib", "numpy", "pandas", "scipy", "seaborn", "statsmodels"))

    with notebook.cell(title="Read MethylTree results"):
        adatas: typing.List[anndata.AnnData] = list()
        clone_tables: typing.List[pandas.DataFrame] = list()

        for file in tqdm.tqdm(iterable=args.input, postfix=step00.main_file()):
            for adata in tqdm.tqdm(iterable=step00.split_methyltree(adata=anndata.read_h5ad(filename=file)), leave=False, postfix=step00.main_file()):
                adatas.append(adata)

                clone_data = adata.uns["accuracy_table"].copy()
                clone_data[step00.sample_column] = adata.uns["parameters"]["Individual"]
                clone_data["Exclude"] = exclude_label(exclude=adata.uns["parameters"]["exclude"])
                clone_tables.append(clone_data)

        clone_data = pandas.concat(objs=clone_tables, ignore_index=True)
        clone_data["clone"] = clone_data["clone"].astype(dtype=str)
        print(clone_data)

        individual_list = sorted(set(clone_data[step00.sample_column]))
        exclude_list = sorted(set(clone_data["Exclude"]), key=lambda x: (x != exclude_label(exclude="none"), x))
        palette = dict(zip(individual_list, itertools.cycle(matplotlib.colors.TABLEAU_COLORS)))
        print("Individual:", individual_list)
        print("Exclude:", exclude_list)

    with notebook.cell(title="Plot settings"):
        matplotlib.use(backend="Agg")
        matplotlib.rcParams.update(step00.matplotlib_parameters)
        seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    with tempfile.TemporaryDirectory() as directory:
        figure_list: typing.List[str] = list()
        pair_inputs: typing.List[typing.Tuple[str, str, pandas.DataFrame, numpy.ndarray, numpy.ndarray]] = list()

        with notebook.cell(title="Heatmaps in MethylTree leaf order", note="Unassigned: cells without a clone label or in a clone with only one cell (the rule MethylTree uses for accuracy)."):
            for adata in tqdm.tqdm(iterable=adatas, postfix=step00.main_file()):
                parameters = adata.uns["parameters"]
                individual = parameters["Individual"]
                name = f"{individual}-{step00.safe_name(name=exclude_label(exclude=parameters['exclude']))}"

                order = numpy.argsort(a=adata.obs["MethylTree_order"].to_numpy())
                clone_labels: typing.List[str] = cell_clone_labels(clones=adata.obs[parameters["clone"]].iloc[order])
                assert (clone_labels == cell_clone_labels(clones=pandas.Series(data=adata.uns["ordered_clone_array"]))), f"{individual}: clone labels differ from MethylTree!!"
                cell_clone: typing.Dict[str, str] = dict(zip(adata.obs_names[order].tolist(), clone_labels))
                clone_order = sorted(dict.fromkeys(x for x in clone_labels if (x != "nan")))
                codes = numpy.array(object=[(0 if (x in step00.missing_values) else clone_order.index(x) + 1) for x in clone_labels])
                print(name, "clone order:", len(clone_order), clone_order)

                local_palette = dict(zip(sorted(set(clone_order)), itertools.cycle(matplotlib.colormaps["tab20"].colors)))
                clone_cmap = matplotlib.colors.ListedColormap(colors=[step00.missing_color] + [local_palette[x] for x in clone_order])
                clone_counts = pandas.Series(data=clone_labels).value_counts()
                legend_labels = [f"Unassigned (n={clone_counts.get(key='nan', default=0)})"] + [f"{clone} (n={clone_counts[clone]})" for clone in clone_order]
                legend_handles = [matplotlib.patches.Patch(facecolor=clone_cmap(X=code), edgecolor="none", label=label) for code, label in enumerate(iterable=legend_labels)]

                similarity = pandas.DataFrame(data=adata.obsm["X_similarity"][order][:, order], index=order, columns=order, dtype=numpy.float64)
                shared_sites = pandas.DataFrame(data=adata.obsm["X_shared_sites"][order][:, order], index=order, columns=order, dtype=int)
                upper = numpy.triu_indices_from(arr=shared_sites.to_numpy(), k=1)

                heatmaps = [("Similarity", similarity, "viridis", numpy.percentile(a=similarity.to_numpy().flatten(), q=40), numpy.percentile(a=numpy.triu(m=similarity.to_numpy(), k=1), q=99), "Corrected similarity"), ("SharedSites", shared_sites, "Blues", 0, numpy.percentile(a=shared_sites.to_numpy()[upper], q=99), "Shared covered regions")]
                for heatmap_name, heatmap_data, heatmap_cmap, heatmap_vmin, heatmap_vmax, colorbar_label in tqdm.tqdm(heatmaps, postfix=step00.main_file()):
                    fig, ax = matplotlib.pyplot.subplot_mosaic(mosaic=[["Legend", ".", "Horizontal_Color", "."], ["Legend", "Vertical_Color", "Heatmap", "Colorbar"]], figsize=(35, 25), width_ratios=[6, 1, 24, 1], height_ratios=[1, 24])
                    seaborn.heatmap(data=heatmap_data, vmin=heatmap_vmin, vmax=heatmap_vmax, cmap=heatmap_cmap, cbar=True, cbar_ax=ax["Colorbar"], xticklabels=False, yticklabels=False, rasterized=True, ax=ax["Heatmap"])
                    ax["Heatmap"].set_xlabel(xlabel=f"Cells (n={len(order)})")
                    ax["Vertical_Color"].set_ylabel(ylabel=f"Cells (n={len(order)})")
                    ax["Vertical_Color"].imshow(X=codes[:, numpy.newaxis], cmap=clone_cmap, vmin=0, vmax=len(clone_order), aspect="auto", interpolation="none", rasterized=True)
                    ax["Horizontal_Color"].imshow(X=codes[numpy.newaxis, :], cmap=clone_cmap, vmin=0, vmax=len(clone_order), aspect="auto", interpolation="none", rasterized=True)
                    for strip in ["Vertical_Color", "Horizontal_Color"]:
                        ax[strip].grid(visible=False)
                        ax[strip].set_xticks(ticks=[])
                        ax[strip].set_yticks(ticks=[])
                    ax["Legend"].set_axis_off()
                    ax["Legend"].legend(handles=legend_handles, title="Clone", loc="upper left", bbox_to_anchor=(0, 1), borderaxespad=0, frameon=False)
                    ax["Colorbar"].tick_params(axis="y", labelsize="xx-small")
                    ax["Colorbar"].set_ylabel(ylabel=colorbar_label, fontsize="x-small")
                    fig.suptitle(t=f"{individual} / {exclude_label(exclude=parameters['exclude'])}\ncell n={parameters['n_cells']} / region n={parameters['n_regions']}")
                    fig.tight_layout()
                    step00.save_figure(fig=fig, name=f"{directory}/{name}-{heatmap_name}-Heatmap", figure_list=figure_list)

                pair_data = pandas.DataFrame(data={"Shared regions": shared_sites.to_numpy()[upper], "Corrected similarity": similarity.to_numpy()[upper], "Pair": numpy.where((codes[upper[0]] > 0) & (codes[upper[0]] == codes[upper[1]]), "Same clone", "Different clone")})
                pair_inputs.append((individual, name, pair_data, shared_sites.to_numpy(), similarity.to_numpy()))

        with notebook.cell(title="Shared covered regions vs corrected similarity", note=f"Spearman p-value treats cell pairs as independent; Mantel p-value permutes cells ({args.permutations} permutations)."):
            for individual, name, pair_data, shared_matrix, similarity_matrix in tqdm.tqdm(iterable=pair_inputs, postfix=step00.main_file()):
                correlation, p_value = scipy.stats.spearmanr(a=pair_data["Shared regions"], b=pair_data["Corrected similarity"])
                mantel_p_value = mantel_spearman_pvalue(matrix_a=shared_matrix, matrix_b=similarity_matrix, permutations=args.permutations, seed=42)
                print(name, "pairs:", pair_data["Pair"].value_counts().to_dict(), "/ Spearman r:", step00.pvalue_format(correlation), "/ p-value", step00.pvalue_format(p_value=p_value), "/ Mantel p:", step00.pvalue_format(p_value=mantel_p_value))

                fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))
                seaborn.regplot(data=pair_data, x="Shared regions", y="Corrected similarity", color=palette[individual], fit_reg=True, robust=True, scatter_kws=dict(s=40, linewidth=0, rasterized=True), ax=ax)
                matplotlib.pyplot.text(x=numpy.median(a=pair_data["Shared regions"]), y=numpy.median(a=pair_data["Corrected similarity"]), s=f"Spearman r={correlation:.3f}, p={step00.pvalue_format(p_value=p_value)}", color="black", path_effects=step00.path_effects, fontsize="x-small", horizontalalignment="center", verticalalignment="center")
                matplotlib.pyplot.tight_layout()
                step00.save_figure(fig=fig, name=f"{directory}/{name}-SharedSites-Similarity-Scatter", figure_list=figure_list)

                g = seaborn.lmplot(data=pair_data, x="Shared regions", y="Corrected similarity", hue="Pair", hue_order=["Different clone", "Same clone"], palette={"Different clone": step00.missing_color, "Same clone": palette[individual]}, height=18, aspect=1, legend=True, facet_kws={"legend_out": False}, scatter=True, fit_reg=True, robust=True, scatter_kws=dict(s=5))
                g.ax.text(x=0.5, y=0.5, s=f"Spearman r={correlation:.3f}, p={step00.pvalue_format(p_value=p_value)}", transform=g.ax.transAxes, color="black", path_effects=step00.path_effects, fontsize="x-small", horizontalalignment="center", verticalalignment="center")
                step00.save_figure(fig=g.figure, name=f"{directory}/{name}-SharedSites-Similarity-Regression", figure_list=figure_list)

        with notebook.cell(title="Observed vs Random"):
            for metric in tqdm.tqdm(metric_list, postfix=step00.main_file()):
                observed_data = clone_data[[step00.sample_column, "Exclude", "clone", "clone_size", metric]].assign(Type="Observed")
                random_data = clone_data[[step00.sample_column, "Exclude", "clone", "clone_size", f"random_{metric}"]].rename(columns={f"random_{metric}": metric}).assign(Type="Random")
                long_data = pandas.concat([observed_data, random_data], ignore_index=True)

                fig, axes = matplotlib.pyplot.subplots(ncols=len(exclude_list), figsize=(18 * len(exclude_list), 18), sharey=True, squeeze=False)
                for ax, exclude in zip(axes.flat, exclude_list):
                    drawing_data = long_data.loc[(long_data["Exclude"] == exclude)]
                    seaborn.violinplot(data=drawing_data, x=step00.sample_column, order=individual_list, y=metric, hue="Type", hue_order=["Observed", "Random"], linewidth=2.5, cut=0, inner="box", legend=(ax is axes.flat[0]), ax=ax)
                    statannotations.Annotator.Annotator(ax, [((individual, "Observed"), (individual, "Random")) for individual in individual_list], data=drawing_data, x=step00.sample_column, order=individual_list, y=metric, hue="Type", hue_order=["Observed", "Random"]).configure(test="Mann-Whitney", text_format="full", hide_non_significant=True, comparisons_correction=None, loc="inside", verbose=1).apply_and_annotate()
                    counts = drawing_data.loc[(drawing_data["Type"] == "Observed")].groupby(step00.sample_column).size()
                    ax.set_xticks(range(len(individual_list)), [f"{individual}\n({counts.get(individual, 0)} clones)" for individual in individual_list], fontsize="x-small", rotation=90)
                    ax.set_title(exclude, fontsize="small")
                    ax.set_xlabel("")
                    ax.set_ylabel(f"Per-clone {metric}")
                matplotlib.pyplot.tight_layout()
                step00.save_figure(fig, f"{directory}/Clone-{metric}-Violin", figure_list)

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

                g = seaborn.relplot(data=pair_data, x=exclude_list[0], y=exclude_list[1], hue=step00.sample_column, hue_order=individual_list, palette=palette, col=step00.sample_column, col_order=individual_list, col_wrap=min(2, len(individual_list)), s=400, alpha=0.6, edgecolor="white", linewidth=1, height=24, aspect=1.0, legend=False)
                counts = pair_data.groupby(step00.sample_column).size()
                for ax, individual in zip(g.axes.flat, individual_list):
                    ax.plot([0, 1], [0, 1], color="black", linestyle="--", linewidth=2, zorder=0)
                    ax.set_title(f"{individual} ({counts.get(individual, 0)} clones)")
                g.set_axis_labels(f"{exclude_list[0]}: {metric}", f"{exclude_list[1]}: {metric}")
                g.set(xlim=(-0.05, 1.05), ylim=(-0.05, 1.05))
                g.tight_layout()
                step00.save_figure(g.figure, f"{directory}/Clone-{metric}-Exclude-Scatter", figure_list)

        with notebook.cell(title="Save"):
            with zipfile.ZipFile(file=args.output, mode="w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zip_file:
                for figure in tqdm.tqdm(iterable=figure_list, postfix=step00.main_file()):
                    zip_file.write(filename=figure, arcname=os.path.basename(p=figure))

        notebook.write_pdf(path=args.output)
