#!/usr/bin/env python3
"""Preprocess H5AD inputs for Lucida.

This script owns the H5AD preprocessing step before `fit <- lucida()`:
contrast subsetting, Scanpy-like QC metrics from raw counts, counts_per_gene, and
canonical covariate recoding. It writes a preprocessed H5AD for the R-only
lucida_contrast.R script.
"""

import argparse
from pathlib import Path

import h5py
import numpy as np
import scanpy as sc
import SCTools as sct

QC_VARS = ["mt", "ribo", "hb"]


def str2bool(value):
    if isinstance(value, bool):
        return value
    value = str(value).strip().lower()
    if value in {"true", "t", "yes", "y", "1", "on"}:
        return True
    if value in {"false", "f", "no", "n", "0", "off"}:
        return False
    raise argparse.ArgumentTypeError(f"Expected boolean value, got {value!r}")


def first_present(obs, candidates, required=True):
    candidates = [c for c in candidates if c]
    for col in candidates:
        if col in obs.columns:
            return col
    if required:
        raise ValueError(f"Missing required obs column. Tried: {', '.join(candidates)}")
    return None


def normalize_sex(values):
    out = np.full(len(values), np.nan, dtype=float)
    for i, value in enumerate(values):
        if value is None:
            continue
        if isinstance(value, (int, float, np.number)):
            val = float(value)
            out[i] = val if np.isfinite(val) else np.nan
            continue
        text = str(value).strip().lower()
        if text in {"male", "m", "1", "true"}:
            out[i] = 1.0
        elif text in {"female", "f", "0", "false"}:
            out[i] = 0.0
    return out


def default_output(args):
    if args.output is not None:
        return args.output
    if args.cohort and args.out_dir:
        return args.out_dir / f"{args.cohort}_preproc.h5ad"
    return args.file.with_name(f"{args.file.stem}_preproc.h5ad")


def add_qc_var_flags(adata):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]", regex=True)


def add_qc_metrics_in_memory(adata):
    add_qc_var_flags(adata)
    if adata.raw is None:
        adata.raw = adata
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=QC_VARS,
        percent_top=False,
        use_raw=True,
        inplace=True,
        log1p=False,
    )
    adata.obs["counts_per_gene"] = (
        adata.obs["total_counts"] / np.maximum(adata.obs["n_genes_by_counts"], 1)
    )


def sparse_encoding(h5_group):
    encoding = h5_group.attrs.get("encoding-type", "")
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    return encoding


def matrix_shape(h5_group):
    shape = h5_group.attrs["shape"]
    return int(shape[0]), int(shape[1])


def csr_row_sums(values, indptr, gene_mask=None):
    n_rows = len(indptr) - 1
    row_nnz = np.diff(indptr)
    out = np.zeros(n_rows, dtype=float)
    nonempty = row_nnz > 0
    if not nonempty.any():
        return out

    starts = indptr[:-1][nonempty]
    if gene_mask is None:
        row_values = values
    else:
        row_values = np.where(gene_mask, values, 0)
    out[nonempty] = np.add.reduceat(row_values, starts)
    return out


def csr_row_n_genes(values, indptr):
    n_rows = len(indptr) - 1
    row_nnz = np.diff(indptr)
    out = np.zeros(n_rows, dtype=float)
    nonempty = row_nnz > 0
    if not nonempty.any():
        return out

    starts = indptr[:-1][nonempty]
    out[nonempty] = np.add.reduceat((values != 0).astype(np.int8, copy=False), starts)
    return out


def calculate_qc_metrics_ondisk(h5ad_path, adata, raw=True, chunk_size=5000):
    """Compute the same QC metrics as Scanpy without materializing full AnnData."""
    add_qc_var_flags(adata)
    with h5py.File(h5ad_path, "r") as f:
        x_path = "raw/X" if raw and "raw/X" in f else "X"
        group = f[x_path]
        if sparse_encoding(group) != "csr_matrix":
            raise ValueError(f"{h5ad_path}:{x_path} must be CSR for on-disk preprocessing")

        n_obs, n_vars = matrix_shape(group)
        if adata.n_obs != n_obs:
            raise ValueError(
                f"{x_path} has {n_obs:,} rows but metadata has {adata.n_obs:,} obs"
            )
        if adata.n_vars != n_vars:
            raise ValueError(
                f"{x_path} has {n_vars:,} columns but metadata has {adata.n_vars:,} vars"
            )

        obs_metrics = {
            "n_genes_by_counts": np.zeros(n_obs, dtype=float),
            "total_counts": np.zeros(n_obs, dtype=float),
        }

        qc_masks = {
            qc_var: adata.var[qc_var].to_numpy(dtype=bool)
            for qc_var in QC_VARS
        }
        for qc_var in QC_VARS:
            obs_metrics[f"total_counts_{qc_var}"] = np.zeros(n_obs, dtype=float)

        var_total_counts = np.zeros(n_vars, dtype=float)
        var_n_cells = np.zeros(n_vars, dtype=float)

        indptr = group["indptr"][:]
        data = group["data"]
        indices = group["indices"]

        for row_start in range(0, n_obs, chunk_size):
            row_end = min(row_start + chunk_size, n_obs)
            data_start = int(indptr[row_start])
            data_end = int(indptr[row_end])
            rel_indptr = indptr[row_start : row_end + 1] - data_start
            values = np.asarray(data[data_start:data_end])
            gene_indices = np.asarray(indices[data_start:data_end])

            obs_slice = slice(row_start, row_end)
            obs_metrics["total_counts"][obs_slice] = csr_row_sums(values, rel_indptr)
            obs_metrics["n_genes_by_counts"][obs_slice] = csr_row_n_genes(
                values, rel_indptr
            )

            for qc_var, mask in qc_masks.items():
                obs_metrics[f"total_counts_{qc_var}"][obs_slice] = csr_row_sums(
                    values, rel_indptr, mask[gene_indices]
                )

            if len(values):
                np.add.at(var_total_counts, gene_indices, values)
                np.add.at(var_n_cells, gene_indices, values != 0)

    total = obs_metrics["total_counts"]
    for qc_var in QC_VARS:
        key = f"total_counts_{qc_var}"
        obs_metrics[f"pct_counts_{qc_var}"] = obs_metrics[key] / total * 100

    obs_metrics["counts_per_gene"] = total / np.maximum(
        obs_metrics["n_genes_by_counts"], 1
    )

    var_metrics = {
        "n_cells_by_counts": var_n_cells,
        "mean_counts": var_total_counts / n_obs,
        "pct_dropout_by_counts": (1 - var_n_cells / n_obs) * 100,
        "total_counts": var_total_counts,
    }
    return obs_metrics, var_metrics


def add_obs_recodes(adata, disease_col, subject_col, age_col, gender_col, pmi_col, chemistry_col, disease):
    disease_values = adata.obs[disease_col].astype(str)
    adata.obs["Disease"] = disease_values
    adata.obs["dx"] = (disease_values == disease).astype(float)
    adata.obs["Sex"] = normalize_sex(adata.obs[gender_col].to_numpy())
    adata.obs["Age"] = adata.obs[age_col].astype(float)
    adata.obs["PMI"] = adata.obs[pmi_col].astype(float)
    adata.obs["Donor"] = adata.obs[subject_col].astype(str).astype("category")
    if chemistry_col is not None:
        adata.obs["Chemistry"] = adata.obs[chemistry_col].astype(str).astype("category")


def resolve_columns(obs, args):
    return {
        "disease": first_present(obs, [args.disease_col, "Disease", "disease"]),
        "subject": first_present(obs, [args.subject_id_col, "individualID", "SubID", "donor_id", "subject_id"]),
        "age": first_present(obs, [args.age_col, "Age", "age"]),
        "gender": first_present(obs, [args.gender_col, "Gender", "gender", "Sex", "sex"]),
        "pmi": first_present(obs, [args.pmi_col, "PMI", "pmi"]),
        "chemistry": first_present(obs, [args.chemistry_col, "Chemistry", "chemistry"], required=False),
    }


def contrast_mask(obs, disease_col, args):
    return obs[disease_col].astype(str).isin([args.control, args.disease]).to_numpy()


def add_preproc_annotations(adata, columns, args, qc_metrics=None, var_metrics=None):
    if qc_metrics is None:
        add_qc_metrics_in_memory(adata)
    else:
        for key, value in qc_metrics.items():
            adata.obs[key] = value
        for key, value in var_metrics.items():
            adata.var[key] = value
    add_obs_recodes(
        adata,
        columns["disease"],
        columns["subject"],
        columns["age"],
        columns["gender"],
        columns["pmi"],
        columns["chemistry"],
        args.disease,
    )
    sct.io.clean_unused_categories(adata)
    return adata


def subset_output_path(args, output):
    if args.cohort:
        return output.parent / f"{args.cohort}.h5ad"
    return output.with_name(f"{output.stem}_subset.h5ad")


def preprocess_ondisk(args, output):
    subset_output = subset_output_path(args, output)
    tmp_subset = subset_output.with_name(f"{subset_output.stem}.tmp{subset_output.suffix}")
    tmp_output = output.with_name(f"{output.stem}.tmp{output.suffix}")
    for path in (tmp_subset, tmp_output):
        if path.exists():
            path.unlink()
    if subset_output.exists():
        print(f"Reusing existing on-disk contrast subset {subset_output}", flush=True)
    else:
        print(f"Reading metadata from {args.file}", flush=True)
        metadata = sct.io.read_everything_but_X(args.file)
        columns = resolve_columns(metadata.obs, args)
        keep = contrast_mask(metadata.obs, columns["disease"], args)
        print(f"Subsetting contrast {args.disease} vs {args.control}: {keep.sum():,}/{metadata.n_obs:,} cells", flush=True)

        # remove metadata to save memory
        del metadata

        if keep.sum() == 0:
            raise ValueError(f"No cells remain for contrast {args.disease} vs {args.control}")

        print(f"Writing on-disk contrast subset {subset_output}", flush=True)
        sct.io.ondisk_subset(
            args.file,
            tmp_subset,
            keep,
            raw=True,
            adata=None,
            chunk_size=args.subset_chunk_size,
        )
        tmp_subset.replace(subset_output)

    subset_metadata = sct.io.read_everything_but_X(subset_output)
    subset_columns = resolve_columns(subset_metadata.obs, args)

    print("Computing QC metrics from subset raw counts in CSR chunks", flush=True)
    qc_metrics, var_metrics = calculate_qc_metrics_ondisk(
        subset_output,
        subset_metadata,
        raw=True,
        chunk_size=args.subset_chunk_size,
    )
    add_preproc_annotations(
        subset_metadata,
        subset_columns,
        args,
        qc_metrics=qc_metrics,
        var_metrics=var_metrics,
    )

    print(f"Writing preprocessed on-disk subset {output}", flush=True)
    sct.io.ondisk_subset(
        subset_output,
        tmp_output,
        None,
        raw=True,
        adata=subset_metadata,
        chunk_size=args.subset_chunk_size,
    )
    tmp_output.replace(output)
    print(f"Saved preprocessed H5AD: {output}", flush=True)


def preprocess_in_memory(args, output):
    print(f"Reading {args.file}", flush=True)
    adata = sc.read_h5ad(args.file)

    columns = resolve_columns(adata.obs, args)
    keep = contrast_mask(adata.obs, columns["disease"], args)
    print(f"Subsetting contrast {args.disease} vs {args.control}: {keep.sum():,}/{adata.n_obs:,} cells", flush=True)
    if keep.sum() == 0:
        raise ValueError(f"No cells remain for contrast {args.disease} vs {args.control}")

    adata = adata[keep].copy()
    add_preproc_annotations(adata, columns, args)

    print(f"Writing {output}", flush=True)
    adata.write_h5ad(output)
    print(f"Saved preprocessed H5AD: {output}", flush=True)


def main():
    parser = argparse.ArgumentParser(description="Preprocess an H5AD for Lucida")
    parser.add_argument("--file", required=True, type=Path, help="Input H5AD")
    parser.add_argument("--output", type=Path, default=None, help="Output preprocessed H5AD")
    parser.add_argument("--out_dir", type=Path, default=None, help="Output directory for {cohort}_preproc.h5ad")
    parser.add_argument("--cohort", default=None, help="Cohort name for default output")
    parser.add_argument("--disease", required=True)
    parser.add_argument("--control", default="CTRL")
    parser.add_argument("--disease_col", default=None)
    parser.add_argument("--subject_id_col", default=None)
    parser.add_argument("--age_col", default=None)
    parser.add_argument("--gender_col", default=None)
    parser.add_argument("--pmi_col", default=None)
    parser.add_argument("--chemistry_col", default=None)
    parser.add_argument("--overwrite", nargs="?", const=True, default=False, type=str2bool)
    parser.add_argument("--ondisk", nargs="?", const=True, default=False, type=str2bool,
                        help="Use SCTools on-disk subsetting before adding preprocessing annotations")
    parser.add_argument("--subset_chunk_size", type=int, default=5000,
                        help="Rows per chunk for SCTools on-disk subsetting/QC")
    args = parser.parse_args()

    if not args.file.is_file():
        raise FileNotFoundError(args.file)
    output = default_output(args)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        if not args.overwrite:
            print(f"Preprocessed H5AD already exists: {output}", flush=True)
            return
        output.unlink()

    if args.ondisk:
        preprocess_ondisk(args, output)
        return

    preprocess_in_memory(args, output)


if __name__ == "__main__":
    main()
