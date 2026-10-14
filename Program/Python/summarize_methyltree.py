import argparse
import typing
import anndata
import pandas
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input MethylTree H5AD files from run_methyltree.py", type=str, nargs="+")
    parser.add_argument("output", help="Output TSV(.gz) file", type=str)

    args = parser.parse_args()

    step00.check_suffixes(args.input, {".h5ad"})
    step00.check_suffix(args.output, {".tsv", ".tsv.gz"})

    rows: typing.List[typing.Dict[str, typing.Any]] = list()
    for file in tqdm.tqdm(args.input):
        adata = anndata.read_h5ad(file)
        accuracy_data = adata.uns["accuracy_table"]
        row = {"File": file.split("/")[-1]} | dict(adata.uns["parameters"]) | {"n_clones": len(accuracy_data), "n_clonal_cells": int(accuracy_data["clone_size"].sum())} | dict(adata.uns["accuracy"])
        rows.append(row)

    output_data = pandas.DataFrame(rows).sort_values(["Individual", "exclude"])
    print(output_data)
    output_data.to_csv(args.output, sep="\t", index=False)
