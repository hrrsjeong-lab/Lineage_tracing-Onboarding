import argparse
import re
import numpy
import pandas
import scanpy
import scipy.sparse
import tqdm
import step00

chromosome_columns = ["chrX", "chrY", "chrM"]


def get_chromosome(region: str) -> str:
    return re.split(r"[:_]", region, maxsplit=1)[0]


def mean_methylation(methylation_sum: pandas.Series, covered: pandas.Series) -> pandas.Series:
    return (methylation_sum / covered).where(covered > 0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input H5AD file", type=str)
    parser.add_argument("output", help="Output TSV(.gz) file", type=str)
    parser.add_argument("--target", help="Target QC", choices=["cell", "region"], required=True)

    args = parser.parse_args()

    step00.check_suffix(args.input, {".h5ad", ".h5", ".hdf5"})
    step00.check_suffix(args.output, {".tsv", ".tsv.gz"})

    input_adata = scanpy.read_h5ad(args.input)
    print(input_adata)

    assert (not scipy.sparse.issparse(input_adata.X)), "X must be dense: uncovered regions are expected as NaN!!"
    assert (step00.sample_column in input_adata.obs.columns), f"{step00.sample_column} is not in obs!!"
    print("Value range:", float(numpy.nanmin(input_adata.X)), float(numpy.nanmax(input_adata.X)))

    chromosomes = pandas.Series(list(map(get_chromosome, input_adata.var_names)), index=input_adata.var_names)
    print(chromosomes.value_counts())
    for column in tqdm.tqdm(chromosome_columns):
        input_adata.var[column] = (chromosomes == column).to_numpy()

    input_adata.layers[step00.covered_column] = (~numpy.isnan(input_adata.X)).astype(numpy.float32)
    input_adata.layers[step00.methylation_column] = numpy.nan_to_num(input_adata.X, nan=0.0)

    if args.target == "cell":
        covered_data, _ = scanpy.pp.calculate_qc_metrics(input_adata, layer=step00.covered_column, qc_vars=chromosome_columns, percent_top=None, log1p=True, inplace=False)
        methylation_data, _ = scanpy.pp.calculate_qc_metrics(input_adata, layer=step00.methylation_column, percent_top=None, log1p=False, inplace=False)

        output_data = pandas.DataFrame(index=input_adata.obs_names)
        output_data[step00.sample_column] = input_adata.obs[step00.sample_column].astype(str)
        output_data["n_covered_regions"] = covered_data["n_genes_by_counts"]
        output_data["log1p_n_covered_regions"] = covered_data["log1p_n_genes_by_counts"]
        output_data["mean_methylation"] = mean_methylation(methylation_data["total_counts"], covered_data["n_genes_by_counts"])
        for column in tqdm.tqdm(chromosome_columns):
            output_data[f"n_covered_regions_{column}"] = covered_data[f"total_counts_{column}"].astype(int)
            output_data[f"pct_covered_regions_{column}"] = covered_data[f"pct_counts_{column}"]

        assert (output_data["n_covered_regions"] == numpy.sum(input_adata.layers[step00.covered_column], axis=1)).all(), step00.default_error_message
        print(output_data.groupby(step00.sample_column).describe().T)
    elif args.target == "region":
        region_list = list()
        for name in tqdm.tqdm(sorted(input_adata.obs[step00.sample_column].unique())):
            adata = input_adata[input_adata.obs[step00.sample_column] == name]
            _, covered_data = scanpy.pp.calculate_qc_metrics(adata, layer=step00.covered_column, percent_top=None, log1p=True, inplace=False)
            _, methylation_data = scanpy.pp.calculate_qc_metrics(adata, layer=step00.methylation_column, percent_top=None, log1p=False, inplace=False)

            region_data = pandas.DataFrame(index=input_adata.var_names)
            region_data[step00.sample_column] = name
            region_data["Chromosome"] = chromosomes
            region_data["n_covered_cells"] = covered_data["n_cells_by_counts"]
            region_data["pct_covered_cells"] = 100.0 - covered_data["pct_dropout_by_counts"]
            region_data["mean_methylation"] = mean_methylation(methylation_data["total_counts"], covered_data["n_cells_by_counts"])

            print(name, f"{(region_data['n_covered_cells'] == 0).sum()} regions without coverage are removed")
            region_list.append(region_data[region_data["n_covered_cells"] > 0])

        output_data = pandas.concat(region_list)
        output_data.index.name = "Region"
        print(output_data.groupby(step00.sample_column).describe().T)
    else:
        raise ValueError(step00.default_error_message)

    print(output_data)
    output_data.to_csv(args.output, sep="\t")
