import argparse
import typing
import anndata
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


def read_individual(filename: str, name: str, sample_data: pandas.DataFrame, clone: str) -> anndata.AnnData:
    print(">", name, filename)
    adata = scanpy.read_h5ad(filename)
    print(adata)

    assert (not scipy.sparse.issparse(adata.X)), f"{name}: X must be dense: uncovered regions are expected as NaN!!"

    covered = ~(numpy.isnan(adata.X))
    print("dtype:", adata.X.dtype)
    print(f"Covered fraction: {covered.mean():.4f}")
    print("Covered regions per cell:", pandas.Series(covered.sum(axis=1)).describe().to_dict())
    print("Covered cells per region:", pandas.Series(covered.sum(axis=0)).describe().to_dict())
    print("Value range:", float(numpy.nanmin(adata.X)), float(numpy.nanmax(adata.X)))

    missing = list(adata.obs_names.difference(sample_data.index))
    assert (not missing), f"{name}: {len(missing)} cells are missing from sample sheet: {missing}!!"
    sheet = sample_data.loc[adata.obs_names]

    for column in tqdm.tqdm(list(sheet.columns)):
        if column in adata.obs.columns:
            mismatch = (adata.obs[column].astype(str) != sheet[column].astype(str)).sum()
            print(f"obs['{column}'] kept; {mismatch} cells differ from sample sheet")
        else:
            adata.obs[column] = as_writable(sheet[column])
            print(f"obs['{column}'] added from sample sheet")

    assert (clone in adata.obs.columns), f"{name}: {clone} is not in obs!!"
    print(adata.obs[clone].value_counts(dropna=False))

    for key in tqdm.tqdm(list(adata.obsp.keys())):
        print(f"obsp['{key}'] dropped (cell-by-cell within {name})")
        del adata.obsp[key]

    for key in tqdm.tqdm(list(adata.obsm.keys())):
        if adata.obsm[key].shape[1] == adata.n_obs:
            print(f"obsm['{key}'] dropped (cell-by-cell within {name})")
            del adata.obsm[key]

    if adata.uns:
        print(f"uns {sorted(adata.uns.keys())} dropped")
        adata.uns.clear()

    return adata


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input H5AD files (MethylTree output)", type=str, nargs="+")
    parser.add_argument("sample", help="Input sample sheet TSV(.gz) file", type=str)
    parser.add_argument("output", help="Output H5AD file", type=str)
    parser.add_argument("--names", help="Individual (mouse) names in the same order as input", type=str, nargs="+", required=True)
    parser.add_argument("--clone", help="Clone column in obs", type=str, default="large_clone_id")
    parser.add_argument("--join", help="Region join across individuals", choices=["outer", "inner"], default="outer")

    args = parser.parse_args()

    step00.check_suffixes(args.input, {".h5ad", ".h5", ".hdf5"})
    step00.check_suffix(args.sample, {".tsv", ".tsv.gz"})
    step00.check_suffix(args.output, {".h5ad", ".h5", ".hdf5"})

    assert (len(args.input) == len(args.names)), f"{len(args.input)} inputs vs {len(args.names)} names!!"
    assert (len(set(args.names)) == len(args.names)), f"Duplicated names: {args.names}!!"

    sample_data = pandas.read_csv(args.sample, sep="\t", index_col="sample")
    print(sample_data)

    assert (not sample_data.index.duplicated().any()), "Duplicated sample in sample sheet!!"

    adatas: typing.Dict[str, anndata.AnnData] = dict()
    for filename, name in tqdm.tqdm(list(zip(args.input, args.names))):
        adatas[name] = read_individual(filename, name, sample_data, args.clone)

    regions = {name: set(adata.var_names) for name, adata in adatas.items()}
    shared = pandas.DataFrame({a: {b: len(regions[a] & regions[b]) for b in args.names} for a in args.names})
    print("Shared regions between individuals:")
    print(shared)
    union, common = set.union(*regions.values()), set.intersection(*regions.values())
    print(f"Regions: union {len(union)}, shared by all {len(common)}, join {args.join}")

    for key in sorted(set().union(*[set(adata.obsm.keys()) for adata in adatas.values()])):
        widths = {adata.obsm[key].shape[1] if key in adata.obsm else -1 for adata in adatas.values()}
        if (len(widths) > 1) or (-1 in widths):
            print(f"obsm['{key}'] dropped (not shared with the same width)")
            for adata in adatas.values():
                if key in adata.obsm:
                    del adata.obsm[key]

    n_obs = sum(adata.n_obs for adata in adatas.values())
    n_vars = len(union) if args.join == "outer" else len(common)
    print(f"Expected X: {n_obs} x {n_vars}, {n_obs * n_vars * adatas[args.names[0]].X.dtype.itemsize / 1024 ** 3:.2f} GiB")

    output_adata = anndata.concat(adatas, join=args.join, merge="same", label="Individual", index_unique=None, fill_value=numpy.nan)
    assert (not output_adata.obs_names.duplicated().any()), "Duplicated cells across individuals!!"

    for column in list(output_adata.obs.columns):
        output_adata.obs[column] = as_writable(output_adata.obs[column])

    clone_values = output_adata.obs[args.clone].astype(str)
    is_missing = output_adata.obs[args.clone].isna() | clone_values.isin(step00.missing_values)
    output_adata.obs["Individual_clone"] = pandas.Series(numpy.where(is_missing, numpy.nan, output_adata.obs["Individual"].astype(str) + ":" + clone_values), index=output_adata.obs_names, dtype=object)
    print(pandas.crosstab(output_adata.obs[args.clone].astype(str), output_adata.obs["Individual"]).head(20))
    print("Cells without clone:", int(is_missing.sum()))
    print("Clones (global):", output_adata.obs["Individual_clone"].nunique())

    covered = ~(numpy.isnan(output_adata.X))
    print(f"Covered fraction (merged): {covered.mean():.4f}")
    print(output_adata)
    output_adata.write_h5ad(args.output, **step00.anndata_compressions)
