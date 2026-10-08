import argparse
import numpy
import pandas
import scanpy
import scipy.sparse
import tqdm
import step00


def as_writable(series: pandas.Series) -> pandas.Series:
    if series.dtype != object:
        return series
    if series.dropna().map(lambda x: isinstance(x, (bool, numpy.bool_))).all():
        return series.astype("boolean")
    return series.map(lambda x: x if pandas.isna(x) else str(x))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input H5AD file (MethylTree output)", type=str)
    parser.add_argument("sample", help="Input sample sheet TSV(.gz) file", type=str)
    parser.add_argument("output", help="Output HDF5 file", type=str)
    parser.add_argument("--name", help="Individual (mouse) name", type=str, required=True)
    parser.add_argument("--clone", help="Clone column in obs", type=str, default="large_clone_id")

    args = parser.parse_args()

    step00.check_suffix(args.input, {".h5ad", ".h5", ".hdf5"})
    step00.check_suffix(args.sample, {".tsv", ".tsv.gz"})
    step00.check_suffix(args.output, {".h5ad", ".h5", ".hdf5"})

    input_adata = scanpy.read_h5ad(args.input)
    print(input_adata)

    assert (not scipy.sparse.issparse(input_adata.X)), "X must be dense: uncovered regions are expected as NaN!!"

    covered = ~numpy.isnan(input_adata.X)
    print("dtype:", input_adata.X.dtype)
    print(f"Covered fraction: {covered.mean():.4f}")
    print("Covered regions per cell:", pandas.Series(covered.sum(axis=1)).describe().to_dict())
    print("Covered cells per region:", pandas.Series(covered.sum(axis=0)).describe().to_dict())
    print("Value range:", float(numpy.nanmin(input_adata.X)), float(numpy.nanmax(input_adata.X)))

    sample_data = pandas.read_csv(args.sample, sep="\t", index_col="sample")
    print(sample_data)

    assert (not sample_data.index.duplicated().any()), "Duplicated sample in sample sheet!!"

    missing = list(input_adata.obs_names.difference(sample_data.index))
    assert (not missing), f"{len(missing)} cells are missing from sample sheet: {list(missing)}!!"
    sample_data = sample_data.loc[input_adata.obs_names]

    for column in tqdm.tqdm(list(sample_data.columns)):
        if column in input_adata.obs.columns:
            mismatch = (input_adata.obs[column].astype(str) != sample_data[column].astype(str)).sum()
            print(f"obs['{column}'] kept; {mismatch} cells differ from sample sheet")
        else:
            input_adata.obs[column] = as_writable(sample_data[column])
            print(f"obs['{column}'] added from sample sheet")

    assert (args.clone in input_adata.obs.columns), f"{args.clone} is not in obs!!"
    print(input_adata.obs[args.clone].value_counts(dropna=False))

    input_adata.obs["Individual"] = args.name
    print(input_adata)
    input_adata.write_h5ad(args.output, **step00.anndata_compressions)
