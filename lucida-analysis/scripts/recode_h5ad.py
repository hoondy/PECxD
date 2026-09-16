#!/usr/bin/env python3

import sys
import shutil
import argparse
import anndata as ad
import h5py
import numpy as np
import scipy.sparse as sp
from pathlib import Path

# Cross-version anndata IO imports. Modern (>=0.10): anndata.io.
# Older: anndata.experimental. Fall back to private path if needed.
try:
  from anndata.io import write_elem, sparse_dataset
except ImportError:
  try:
    from anndata.experimental import write_elem, sparse_dataset
  except ImportError:
    from anndata._io.specs import write_elem
    from anndata._core.sparse_dataset import sparse_dataset


# -- streaming sparse format conversion (no full materialization) -------------

def _matrix_shape(h5_group):
  """Return (n_rows, n_cols) as a Python int tuple from an anndata sparse group.
  Handles np.int64/list/tuple shape attr variations across anndata versions."""
  shape = h5_group.attrs["shape"]
  return (int(shape[0]), int(shape[1]))


def _sparse_format(h5_group):
  """Return 'csr', 'csc', or None for a sparse-encoded h5 group.
  None if the group is not sparse-encoded."""
  enc = h5_group.attrs.get("encoding-type", b"")
  if isinstance(enc, bytes):
    enc = enc.decode()
  if enc == "csr_matrix":
    return "csr"
  if enc == "csc_matrix":
    return "csc"
  return None


def _stream_swap_layout(parent_out, name, in_group, target_format,
                        col_block=64, value_chunk=8_000_000, verbose=True):
  """Convert a sparse h5 group between CSR and CSC layouts without ever holding
  the full source matrix in memory.

  Adopts the same on-disk-write pattern as SCTools.io.ondisk_subset:
  seed an empty sparse matrix at parent_out[name] via write_elem, then append
  blocks of the converted matrix using sparse_dataset(...).append(...).

  Algorithm (CSR -> CSC; CSC -> CSR is symmetric, axes swapped):
    For each block of `col_block` destination-major axis (cols if dst is CSC,
    rows if dst is CSR):
      Walk the source in row-chunks (CSR) or col-chunks (CSC).
      For each source chunk, mask nonzeros whose minor coord falls in the
      current destination block.
      Assemble the kept triples into an in-memory sparse block of shape
      (n_minor_dst, block_width) for CSC dst, or (block_width, n_minor_dst)
      for CSR dst.
      Append the block to the destination via sparse_dataset.append.

  Memory is bounded by `block_size_nnz × ~12 bytes` for buffers plus the
  block's scipy sparse matrix (≈ block_size_nnz × 16 bytes). For a uniform
  matrix, block_size_nnz ≈ nnz_total × col_block / n_major_dst.

  Args:
    parent_out: h5py.Group that will contain the output sparse matrix
    name:       child group name within parent_out
    in_group:   h5py.Group containing the source sparse matrix
                (`data`, `indices`, `indptr`, attr `shape`, attr `encoding-type`)
    target_format: 'csr' or 'csc'
    col_block:  number of destination-major positions per pass
    value_chunk: hint for source-chunk size (in nonzeros) during the inner walk
  """
  src_fmt = _sparse_format(in_group)
  if src_fmt is None:
    raise ValueError(
      f"`{in_group.name}` is not a sparse-encoded group "
      f"(encoding-type={in_group.attrs.get('encoding-type')!r}). "
      f"Expected csr_matrix or csc_matrix.")
  if src_fmt == target_format:
    raise ValueError(f"`{in_group.name}` is already {src_fmt}; no swap needed")
  if {src_fmt, target_format} != {"csr", "csc"}:
    raise ValueError(f"can only swap csr<->csc, got {src_fmt}->{target_format}")

  shape = _matrix_shape(in_group)
  n_rows, n_cols = shape
  src_data = in_group["data"]
  src_indices = in_group["indices"]
  src_indptr = in_group["indptr"][...]   # small, load once
  nnz = int(src_data.shape[0])
  data_dtype = src_data.dtype

  # In CSR, src_indptr has len n_rows+1 (one entry per row); src_indices holds columns.
  # In CSC, src_indptr has len n_cols+1; src_indices holds rows.
  # "src_major" = the axis indexed by src_indptr (rows for CSR, cols for CSC).
  # "dst_major" = the axis indexed by dst_indptr (cols for CSC, rows for CSR).
  # By construction src_minor == dst_major (we're swapping the layout's role).
  if src_fmt == "csr":
    n_major_src = n_rows
    n_major_dst = n_cols
  else:
    n_major_src = n_cols
    n_major_dst = n_rows

  if verbose:
    print(f"  Streaming {src_fmt}->{target_format}: shape={shape}, nnz={nnz:,}")

  # Pick a source-chunk size targeted to value_chunk nonzeros worth of data
  avg_nnz_per_major = max(1, nnz // max(1, n_major_src))
  src_chunk_size = max(1024, value_chunk // avg_nnz_per_major)

  n_blocks = (n_major_dst + col_block - 1) // col_block
  if verbose:
    print(f"  {n_blocks} block(s) of {col_block} dst-major positions; "
          f"src_chunk_size={src_chunk_size:,}")

  # Seed the destination with an empty sparse matrix of the correct shape and
  # dtype. write_elem encodes encoding-type, encoding-version, and shape attrs
  # the way anndata expects.
  #
  # We always use int64 for indptr. scipy defaults to int32 for empty-shape
  # constructions, but if cumulative nnz crosses 2**31 during append,
  # sparse_dataset.append raises OverflowError. The disk overhead is ~5% on
  # gzip-compressed h5ads, which we accept for simpler logic. indices stays
  # int32 — that would only overflow with > 2.1B rows/cols, which doesn't
  # happen in real single-cell data.
  if target_format == "csc":
    seed = sp.csc_matrix((n_rows, 0), dtype=data_dtype)
  else:
    seed = sp.csr_matrix((0, n_cols), dtype=data_dtype)
  seed.indptr = seed.indptr.astype(np.int64)
  write_elem(parent_out, name, seed)
  out_sd = sparse_dataset(parent_out[name])

  # Main loop: produce one block of dst at a time, append it.
  log_every = max(1, n_blocks // 10)
  for block_idx, c0 in enumerate(range(0, n_major_dst, col_block)):
    c1 = min(c0 + col_block, n_major_dst)
    block_width = c1 - c0

    # Triple buffers (filled across source chunks, then assembled into sparse block)
    blk_src_major = []   # source-major coord of each kept nonzero (becomes dst.indices)
    blk_dst_major = []   # destination-major coord (within the block, 0..block_width-1)
    blk_data = []

    for r0 in range(0, n_major_src, src_chunk_size):
      r1 = min(r0 + src_chunk_size, n_major_src)
      d0 = int(src_indptr[r0])
      d1 = int(src_indptr[r1])
      if d1 == d0:
        continue
      idx_block = src_indices[d0:d1]
      val_block = src_data[d0:d1]
      mask = (idx_block >= c0) & (idx_block < c1)
      if not mask.any():
        continue
      sel_idx = idx_block[mask]            # absolute dst-major coord
      sel_val = val_block[mask]
      # recover source-major coord per kept nonzero from local indptr
      local_indptr = src_indptr[r0:r1+1] - d0
      row_lens = np.diff(local_indptr)
      src_majors_chunk = np.repeat(np.arange(r0, r1, dtype=np.int64), row_lens)
      sel_src_major = src_majors_chunk[mask]

      blk_src_major.append(sel_src_major)
      blk_dst_major.append((sel_idx - c0).astype(np.int64, copy=False))
      blk_data.append(sel_val)

    if not blk_data:
      # empty block: append an empty sparse matrix of the right shape
      if target_format == "csc":
        empty = sp.csc_matrix((n_rows, block_width), dtype=data_dtype)
      else:
        empty = sp.csr_matrix((block_width, n_cols), dtype=data_dtype)
      out_sd.append(empty)
      if verbose and block_idx % log_every == 0:
        print(f"    block {block_idx+1}/{n_blocks}: dst-major [{c0}:{c1}) -> empty")
      continue

    rows_arr = np.concatenate(blk_src_major)
    cols_arr = np.concatenate(blk_dst_major)
    data_arr = np.concatenate(blk_data)

    # Build the block in the target format. scipy sorts indices canonically.
    # Block indptr/indices stay int32 (sparse_dataset.append widens indptr to
    # int64 when appending to the int64 seed; indices stays int32 throughout).
    if target_format == "csc":
      block = sp.csc_matrix(
        (data_arr, (rows_arr, cols_arr)),
        shape=(n_rows, block_width),
      )
    else:
      # dst is CSR: dst-major is row, src-major (originally a column) is the column
      block = sp.csr_matrix(
        (data_arr, (cols_arr, rows_arr)),
        shape=(block_width, n_cols),
      )

    block.sort_indices()
    out_sd.append(block)

    if verbose and block_idx % log_every == 0:
      print(f"    block {block_idx+1}/{n_blocks}: dst-major [{c0}:{c1}) "
            f"-> {block.nnz:,} nnz")

  if verbose:
    print(f"  Done: final shape {out_sd.shape}, nnz {out_sd.group['data'].shape[0]:,}")


def _streaming_convert(in_path, out_path, target_format, libsize_values=None,
                       skip_libsize=False, col_block=64, value_chunk=8_000_000,
                       verbose=True):
  """Build a new h5ad at out_path with /X and /raw/X converted to target_format,
  copying all other groups verbatim. Output file size ≈ input file size (no
  HDF5 free-space waste).

  Uses _stream_swap_layout to do the conversion without loading the full matrix.
  libsize_values, if provided, is written to /obs/libSize.
  """
  target_format = target_format.lower()
  if target_format not in ("csr", "csc"):
    raise ValueError(f"target_format must be csr or csc, got {target_format}")

  with h5py.File(in_path, "r") as fin, h5py.File(out_path, "w") as fout:
    # Copy file-level attrs
    for k, v in fin.attrs.items():
      fout.attrs[k] = v

    # First pass: copy everything except /X and /raw/X verbatim.
    for name in fin:
      if name == "X":
        # we'll handle below
        continue
      if name == "raw":
        # copy raw verbatim except raw/X
        if verbose:
          print(f"Copying /raw (excluding /raw/X)...")
        raw_in = fin["raw"]
        raw_out = fout.create_group("raw")
        for k, v in raw_in.attrs.items():
          raw_out.attrs[k] = v
        for sub in raw_in:
          if sub == "X":
            continue
          raw_in.copy(sub, raw_out)
        continue
      if verbose:
        print(f"Copying /{name}...")
      fin.copy(name, fout)

    # Now write converted /X and /raw/X
    for x_path in ("X", "raw/X"):
      if x_path not in fin:
        continue
      g_in = fin[x_path]
      src_fmt = _sparse_format(g_in)
      parent_path, base = x_path.rsplit("/", 1) if "/" in x_path else ("", x_path)
      parent_out = fout[parent_path] if parent_path else fout

      if src_fmt is None:
        if verbose:
          print(f"/{x_path} not sparse "
                f"(encoding={g_in.attrs.get('encoding-type')!r}), copying verbatim")
        g_in.parent.copy(base, parent_out)
        continue
      if src_fmt == target_format:
        if verbose:
          print(f"/{x_path} already {target_format}, copying verbatim")
        g_in.parent.copy(base, parent_out)
        continue

      if verbose:
        print(f"Converting /{x_path} ({src_fmt} -> {target_format})...")
      _stream_swap_layout(parent_out, base, g_in, target_format,
                          col_block=col_block, value_chunk=value_chunk,
                          verbose=verbose)

    # libSize: write to /obs/libSize
    if not skip_libsize and libsize_values is not None:
      if verbose:
        print(f"Writing libSize to /obs/libSize...")
      obs = fout["obs"]
      if "libSize" in obs:
        del obs["libSize"]
      # Use write_elem so encoding-type/version/etc. match the anndata version.
      write_elem(obs, "libSize", np.asarray(libsize_values, dtype=np.float64))
      # update column-order
      if "column-order" in obs.attrs:
        order = list(obs.attrs["column-order"])
        order = [c.decode() if isinstance(c, bytes) else c for c in order]
        if "libSize" not in order:
          order.append("libSize")
          obs.attrs["column-order"] = np.asarray(order, dtype=h5py.string_dtype())


def _row_sum(M, chunk_size=10000):
  """Sum each row, without materializing the full matrix for backed datasets.

  For backed CSR (the common case for libSize-from-raw), we read directly from
  the HDF5 `data` + `indptr` arrays — `data[indptr[i]:indptr[i+1]]` is row i's
  nonzeros. No indices, no full load.

  For backed CSC, row sums require column-walking, so we chunk by reading
  `data`/`indices` slices keyed off `indptr` rather than relying on the
  backed object's __getitem__ (which has rough edges across anndata versions).
  """
  fmt = getattr(M, "format", None)

  # Backed CSR: stream rows from data + indptr
  if fmt == "csr" and hasattr(M, "group"):
    g = M.group
    indptr = g["indptr"][...]
    data = g["data"]
    out = np.zeros(M.shape[0], dtype=np.float64)
    # chunk row-ranges so we read `data` in reasonable HDF5 slabs
    for r0 in range(0, M.shape[0], chunk_size):
      r1 = min(r0 + chunk_size, M.shape[0])
      d0, d1 = int(indptr[r0]), int(indptr[r1])
      if d1 == d0:
        continue
      block = data[d0:d1]
      # row lengths within this chunk
      row_lens = np.diff(indptr[r0:r1+1])
      # Build local row ids in [0, r1-r0), then bincount with weights.
      # Robust to empty rows anywhere (including last row of chunk).
      local_row_ids = np.repeat(np.arange(r1 - r0, dtype=np.intp), row_lens)
      out[r0:r1] = np.bincount(local_row_ids, weights=block, minlength=r1 - r0)
    return out

  # Backed CSC: walk columns, accumulate into row sums
  if fmt == "csc" and hasattr(M, "group"):
    g = M.group
    indptr = g["indptr"][...]
    data = g["data"]
    indices = g["indices"]
    out = np.zeros(M.shape[0], dtype=np.float64)
    n_cols = M.shape[1]
    for c0 in range(0, n_cols, chunk_size):
      c1 = min(c0 + chunk_size, n_cols)
      d0, d1 = int(indptr[c0]), int(indptr[c1])
      if d1 == d0:
        continue
      np.add.at(out, indices[d0:d1], data[d0:d1])
    return out

  # Backed dense (h5py.Dataset, no .format attr, no .sum)
  if hasattr(M, "shape") and not (sp.issparse(M) or isinstance(M, np.ndarray)):
    out = np.zeros(M.shape[0], dtype=np.float64)
    for r0 in range(0, M.shape[0], chunk_size):
      r1 = min(r0 + chunk_size, M.shape[0])
      out[r0:r1] = np.asarray(M[r0:r1]).sum(axis=1)
    return out

  # In-memory sparse / dense
  return np.asarray(M.sum(axis=1)).ravel()


def _is_format(M, fmt):
  """True if M is already in the requested format ('csr' or 'csc'), whether
  it's an in-memory scipy sparse matrix or a backed sparse dataset."""
  fmt = fmt.lower()
  # backed sparse datasets expose a .format attribute ('csr' or 'csc')
  if hasattr(M, "format") and not sp.issparse(M):
    return M.format == fmt
  if fmt == "csr":
    return sp.isspmatrix_csr(M)
  return sp.isspmatrix_csc(M)


def main():

  # Define arguments
  #-------------------

  parser = argparse.ArgumentParser(
      description="Convert an AnnData .h5ad file so that X (and raw/X) is stored in CSC sparse format (v1.1)"
  )

  parser.add_argument("--input",
      type=Path, 
      required=True,
      help="Input .h5ad file")

  parser.add_argument("--output",      
      type=Path, 
      required=True,
      help="Output .h5ad file")

  parser.add_argument("--ondisk", 
      action="store_true",
      help="Use file-backed mode to reduce memory usage")

  parser.add_argument("--sortBy", 
      default=None,
      help="Cols to sort by in _decreasing_ order of importance")

  parser.add_argument(
      "--compression",
      default=None,
      choices=[None, "None", "gzip", "lzf"],
      help="Optional compression for output file (gzip or lzf). Default: None")

  parser.add_argument(
      "--format",
      default="CSR",
      choices=["CSR", "CSC"],
      help="Store sparse count matrix (X and raw/X) in CSR or CSC format.  CSR allows faster access to cells, CSC gives faster access to genes. Default: CSR")

  parser.add_argument("--noLibSize", 
      action="store_true",
      help="Skip computing libSize for each cell")

  parser.add_argument("--colBlock",
      type=int, default=64,
      help="(streaming mode only) Number of destination-major columns processed "
           "per pass. Lower = less peak memory, more passes. For very large "
           "h5ads tune down (e.g. 16 or 8). Default: 64")

  parser.add_argument("--valueChunk",
      type=int, default=8_000_000,
      help="(streaming mode only) Number of nonzeros per HDF5 read in pass 1. "
           "Default: 8000000")

  # parse args
  args = parser.parse_args()

  # if --ondisk set backed to "r"
  backed = None
  if args.ondisk:
    backed = "r"

  if args.compression == "None":
    args.compression = None

  # read file
  print("Read file...") 
  adata = ad.read_h5ad(args.input, backed=backed) 

  # Find counts entry
  if adata.X is not None:
      print("Using AnnData X matrix...")
  else:
    if 'X' in adata.layers: 
      print("Using AnnData layers/X matrix...")
      adata.X = adata.layers['X']
    elif 'counts' in adata.layers: 
      print("Using AnnData layers/counts matrix...")
      adata.X = adata.layers['counts']

  # sort cells by type  
  if args.sortBy != None:

    print("Sorting...") 
    fields = args.sortBy.split(",")

    # check if all fields are present
    missing = set(fields) - set(adata.obs.columns)

    if missing:
      print(f"Missing columns: {missing}")
      sys.exit(2)     

    # get sorted order
    # reverse since 1st sorted index is the last one for lexsort
    idx = np.lexsort( 
      keys = tuple(adata.obs[c].to_numpy() for c in fields[::-1]) )

    # apply reordering
    adata = adata[idx,:]

  # ---- streaming conversion path: --ondisk + format change needed ----
  # We have to decide this BEFORE doing format-change in memory below.
  needs_x_swap = (adata.X is not None) and (not _is_format(adata.X, args.format.lower()))
  needs_raw_swap = (adata.raw is not None) and (not _is_format(adata.raw.X, args.format.lower()))
  use_streaming = backed is not None and (needs_x_swap or needs_raw_swap)

  if use_streaming:
    if args.sortBy is not None:
      sys.exit("ERROR: --sortBy is not supported with --ondisk + format conversion. "
               "Re-run without --ondisk, or pre-sort and re-run.")
    if args.compression is not None:
      print(f"NOTE: --compression={args.compression} ignored in streaming mode "
            f"(output preserves input file's compression).")

    if args.output.is_file():
      args.output.unlink(missing_ok=True)

    libsize_values = None
    if not args.noLibSize:
      print("Compute libSize...")
      if adata.raw is not None:
        libsize_values = _row_sum(adata.raw.X)
      else:
        libsize_values = _row_sum(adata.X)

    # Close the backed input before copying the file (avoid HDF5 lock conflicts)
    in_path = str(args.input)
    out_path = str(args.output)
    adata.file.close()
    del adata

    print(f"Streaming format conversion to {args.format} (no full matrix load)...")
    _streaming_convert(in_path, out_path, args.format,
                       libsize_values=libsize_values,
                       skip_libsize=args.noLibSize,
                       col_block=args.colBlock,
                       value_chunk=args.valueChunk,
                       verbose=True)
    print("Done.")
    return

  # compute libSize (works in backed mode via chunked row-sum)
  if not args.noLibSize:
    print("Compute libSize...")
    if adata.raw is not None:
      adata.obs['libSize'] = _row_sum(adata.raw.X)
    else:
      adata.obs['libSize'] = _row_sum(adata.X)

  if args.format == "CSC":

    if not _is_format(adata.X, "csc"):
      if backed is None:
        adata = adata.copy()
      else:
        sys.exit("ERROR: --ondisk mode is strict: converting .X to CSC requires loading into memory. Re-run without --ondisk.")

    if not _is_format(adata.X, "csc"):
      print("Converting .X to CSC sparse format...")
      # convert matrix type
      adata.X = sp.csc_matrix(adata.X)

      if adata.raw is not None and not _is_format(adata.raw.X, "csc"):
        print("Converting .raw.X to CSC sparse format...")
        raw = adata.raw.to_adata()
        raw.X = sp.csc_matrix(raw.X)
        adata.raw = raw

  if args.format == "CSR":
    # X
    if not _is_format(adata.X, "csr"):
      print("Converting .X to CSR sparse format...")
      if backed is None:
        adata = adata.copy()
      else:
        sys.exit("ERROR: --ondisk mode is strict: converting .X to CSR requires loading into memory. Re-run without --ondisk.")
      
      # convert matrix type
      adata.X = sp.csr_matrix(adata.X)

      if adata.raw is not None and not _is_format(adata.raw.X, "csr"):
        print("Converting .raw.X to CSR sparse format...")
        raw = adata.raw.to_adata()
        raw.X = sp.csr_matrix(raw.X)
        adata.raw = raw

  if args.output.is_file():
    args.output.unlink(missing_ok=True)

  print("Writing H5AD...") 

  adata.write_h5ad( args.output, compression=args.compression )



if __name__ == "__main__":
  main()
