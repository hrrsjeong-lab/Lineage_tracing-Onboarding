import datetime
import os
import re
import shlex
import sys
import typing
import warnings
import anndata
import matplotlib.patches
import matplotlib.patheffects
import matplotlib.pyplot
import matplotlib.transforms
import numpy
import tqdm

cpus_error_message = "CPUs must be positive!!"
default_error_message = "Something went wrong!!"

matplotlib_parameters: typing.Dict[str, typing.Any] = {"font.size": 50, "axes.labelsize": 50, "axes.titlesize": 60, "figure.titlesize": 60, "xtick.labelsize": 40, "ytick.labelsize": 40, "legend.fontsize": 20, "legend.title_fontsize": 25, "figure.dpi": 600, "text.color": "black", "font.family": "sans-serif", "pdf.fonttype": 42, "ps.fonttype": 42, "pdf.compression": 9}
suffixes = ["pdf", "png"]

epsilon = sys.float_info.epsilon
float_minimum = sys.float_info.min

arrowprops = {"arrowstyle": "-", "color": "silver", "linewidth": 0.5}
path_effects = [matplotlib.patheffects.withStroke(linewidth=5, foreground="white")]

anndata_compressions = {"compression": "gzip", "compression_opts": 4}
missing_values = {"", "nan", "nan0", "NaN", "None", "none"}
missing_color = "lightgray"

sample_column = "Individual"
covered_column = "Covered"
methylation_column = "Methylation"


def main_file(line: bool = True) -> str:
    path = getattr(sys.modules.get("__main__"), "__file__", None)
    if path is None:
        return "<interactive>"
    path = path if path.startswith("<") else os.path.abspath(path)
    if not line:
        return path
    frame = sys._getframe(1)
    while (frame is not None) and (frame.f_globals.get("__name__") != "__main__"):
        frame = frame.f_back
    return path if (frame is None) else f"{path}:{frame.f_lineno}"


def check_suffix(filename: str, suffixes: typing.Set[str]) -> None:
    filename = filename.lower()
    for suffix in tqdm.tqdm(suffixes, postfix=main_file()):
        if filename.endswith(suffix.lower()):
            return
    raise ValueError(f"{filename} must be ended with one of {sorted(suffixes)}!!")


def check_suffixes(filenames: typing.List[str], suffixes: typing.Set[str]) -> None:
    for filename in tqdm.tqdm(filenames, postfix=main_file()):
        check_suffix(filename, suffixes)
    return


def check_cpus(cpus: int) -> None:
    assert (cpus >= 1), cpus_error_message


def confidence_ellipse(x: typing.List[float], y: typing.List[float], ax, n_std: float = 2.0, facecolor: str = "none", **kwargs) -> matplotlib.patches.Patch:
    assert (len(x) == len(y)), "x and y must be the same size!!"

    if len(x) < 3:
        return matplotlib.patches.Ellipse((0, 0), width=0, height=0, facecolor=facecolor, **kwargs)

    cov = numpy.cov(x, y)
    pearson = cov[0, 1] / numpy.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = numpy.sqrt(1 + pearson)
    ell_radius_y = numpy.sqrt(1 - pearson)
    ellipse = matplotlib.patches.Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2, facecolor=facecolor, **kwargs)

    scale_x = numpy.sqrt(cov[0, 0]) * n_std
    mean_x = numpy.mean(x)

    scale_y = numpy.sqrt(cov[1, 1]) * n_std
    mean_y = numpy.mean(y)

    transf = matplotlib.transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(mean_x, mean_y)
    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


def pvalue_format(p_value: float) -> str:
    assert (0.0 <= p_value <= 1.0), f"p={p_value} is not a valid p-value!!"

    if p_value >= 1e-3:
        return f"{p_value:.3f}"
    elif p_value < 1e-3:
        return f"{p_value:.1e}"

    raise ValueError(default_error_message)


def figure_metadata(extra: typing.Optional[typing.Mapping[str, typing.Any]] = None) -> typing.Dict[str, str]:
    information = {"Script": main_file(), "Command": shlex.join(sys.argv), "Timestamp": datetime.datetime.now().astimezone().isoformat(timespec="seconds")}
    if extra is not None:
        information.update({str(key): str(value) for key, value in extra.items()})
    return information


def save_figure(fig: matplotlib.figure.Figure, name: str, figure_list: typing.List[str], extra_metadata: typing.Optional[typing.Mapping[str, typing.Any]] = None) -> None:
    information = figure_metadata(extra_metadata)
    for suffix in tqdm.tqdm(suffixes, postfix=main_file()):
        assert (f"{name}.{suffix}" not in figure_list), f"{name}.{suffix} is duplicated!!"
        figure_list.append(f"{name}.{suffix}")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Unknown infodict keyword", category=UserWarning)
            fig.savefig(figure_list[-1], metadata=information)
    matplotlib.pyplot.close(fig)


def chromosome_order(chromosome: str) -> typing.Tuple[int, str]:
    name = re.sub(r"^chr", "", chromosome)
    return (int(name), "") if name.isdigit() else (1000, name)


def safe_name(name: str) -> str:
    return name.replace(' ', '_').replace(',', '_')


def split_methyltree(adata: anndata.AnnData) -> typing.List[anndata.AnnData]:
    results: typing.List[anndata.AnnData] = list()
    for individual, result in tqdm.tqdm(adata.uns["MethylTree"].items(), total=len(adata.uns["MethylTree"]), postfix=main_file()):
        mask = (adata.obs[sample_column] == individual).to_numpy()
        results.append(anndata.AnnData(obs=adata.obs.loc[mask].copy(), obsm={key: numpy.asarray(adata.obsp[key])[numpy.ix_(mask, mask)] for key in adata.obsp.keys()}, uns=dict(result)))
    return results
