# Lucida cross-disorder analysis

Codebase documenting the Lucida workflow used for 10 disease-versus-control contrasts. Run each step manually.

Analysis settings are listed in `manifests/lucida_contrast_runs.csv`. Input files go in `data/` and outputs in `results/`. Data and results are not distributed.

## Setup

Run commands from this repository's root. Use Python 3.12 and an R environment with Lucida and GenomicDataStream installed.

```bash
python -m pip install -r requirements.txt
mkdir -p data results
```

Python dependencies are listed in `requirements.txt`. SCTools is installed from `hoondy/SCTools`.

R setup (in R; system build dependencies may also be required):

```r
install.packages(c("BiocManager", "remotes", "dplyr", "ggplot2", "ggrepel", "arrow", "metafor"))
BiocManager::install(c("SingleCellExperiment", "S4Vectors", "DelayedArray", "rhdf5"))
remotes::install_github("GabrielHoffman/GenomicDataStream@896288b2a83b5e49047b53d76a2c5b5fc38cda69")
remotes::install_github("GabrielHoffman/lucida@19909ed9cd38baeea810c7375575025460dabb7d")
```

The commands pin Lucida and GenomicDataStream to specific commits. `arrow` enables Parquet output; without it, the scripts write `.csv.gz`. `zellkonverter` is needed only for optional residualized output, disabled here. Python and R may live in separate environments: replace `python` and `Rscript` below with their respective executable paths.

## Inputs and contrasts

Place the original `*_finalmeta.h5ad` files in `data/` (or symlink individual files). Do not use previously subsetted/CSC files as step 1 inputs. On-disk preprocessing expects sparse CSR raw counts, matching raw and main gene axes, and gene-symbol variable names for mitochondrial/ribosomal/hemoglobin QC flags. Supply raw counts in `raw/X` (or `X` when no raw matrix exists).

Required observation columns for these runs: `Disease`, `individualID`, `Age`, `Gender`, `PMI`, `Chemistry`, and `subclass`. Age and PMI must be numeric-convertible. Sex is recoded as male=1, female=0; unknown values become missing. Each donor should belong to one diagnosis. The code does not replace missing covariates or perform an upstream QC/filtering workflow.

| Cohort | Input filename | Disease |
| --- | --- | --- |
| CMC_SCZ | CMC_finalmeta.h5ad | SCZ |
| SZBD_SCZ | SZBD_finalmeta.h5ad | SCZ |
| SZBD_BD | SZBD_finalmeta.h5ad | BD |
| UCLA_ASD | UCLA_finalmeta.h5ad | ASD |
| Velmeshev_ASD | Velmeshev_finalmeta.h5ad | ASD |
| MultiomeBrain_SCZ | MultiomeBrain_finalmeta.h5ad | SCZ |
| MultiomeBrain_BD | MultiomeBrain_finalmeta.h5ad | BD |
| PTSDBrainomics_PTSD | PTSDBrainomics_finalmeta.h5ad | PTSD |
| PTSDBrainomics_MDD | PTSDBrainomics_finalmeta.h5ad | MDD |
| Turecki_MDD | Turecki_finalmeta.h5ad | MDD |

All contrasts use `CTRL`. The full settings are in `manifests/lucida_contrast_runs.csv`. This is a reference manifest, not an automatically executed job list. Its `file` column names the original input; Lucida itself uses the derived CSC file. `celltype_col` maps to `--cell_type_col`; `dx_col`, `Individual_col`, and `Sex_col` map to preprocessing's `--disease_col`, `--subject_id_col`, and `--gender_col`.

## Run one contrast manually

Set these variables in a Bash shell, then run steps 1–3 in that same shell. Repeat with each row above. The remaining settings are identical across all 10 manifest rows.

```bash
cohort=CMC_SCZ
disease=SCZ
input=data/CMC_finalmeta.h5ad
```

### 1. Subset and preprocess

Subset disease and control cells, compute QC covariates from counts, and create `dx` (disease=1, control=0), `Donor`, `Age`, `Sex`, `PMI`, and `Chemistry` annotations.

```bash
python scripts/lucida_preproc.py \
  --file "$input" --cohort "$cohort" --out_dir data \
  --disease "$disease" --control CTRL \
  --disease_col Disease --subject_id_col individualID \
  --age_col Age --gender_col Gender --pmi_col PMI --chemistry_col Chemistry \
  --ondisk TRUE --subset_chunk_size 5000
```

Writes `data/${cohort}.h5ad` (intermediate subset) and `data/${cohort}_preproc.h5ad`. Existing preprocessed output is skipped; an existing intermediate subset is reused. When changing the input or contrast, use a new cohort/output location to avoid reusing stale files. `--overwrite TRUE` replaces the preprocessed output but does **not** regenerate an existing intermediate subset.

### 2. Convert CSR to CSC

Lucida streams genes, so convert the smaller preprocessed file rather than the original cohort file. This also computes `obs/libSize` from raw counts.

```bash
python scripts/recode_h5ad.py \
  --input "data/${cohort}_preproc.h5ad" \
  --output "data/${cohort}_preproc_CSC.h5ad" \
  --format CSC --ondisk
```

Both `X` and `raw/X` are converted. Streaming reduces memory usage but repeatedly scans the input and can be slow. Keep input and output paths distinct. This step replaces an existing output; do not run it against a file being used by R.

### 3. Fit Lucida and variance partitioning

```bash
Rscript scripts/lucida_contrast.R \
  --file "data/${cohort}_preproc_CSC.h5ad" \
  --cohort "$cohort" --out_dir results \
  --disease "$disease" --control CTRL \
  --formula '~ dx + (1|Donor) + Age + Sex + PMI + Chemistry + counts_per_gene + pct_counts_ribo' \
  --null_formula '~ (1|Donor)' \
  --cell_type_col subclass --coef dx \
  --nthreads 8 --chunk_size 1000 \
  --write_varpart TRUE --write_residualized FALSE --varpart_genes '' \
  --plot_ncol 6 --plot_dpi 300 --plot_width 20 --plot_height 16
```

Fits a donor random-intercept model separately by `subclass`. The tested coefficient is `dx`; donor-only null fits are used for variance partitioning. Label permutation is disabled by default. Step 3 refits and replaces same-named outputs; use a fresh results directory for a different analysis.

Outputs in `results/`, prefixed by cohort:

- `_deg.parquet` (or `_deg.csv.gz`): differential-expression table.
- `_lucida_fit.rds`, `_lucida_fit_null.rds`, `_varpart.rds`: fitted models and variance partitioning.
- `_MA.png`, `_volcano.png`, `_dispEsts.png`, `_varpart_violin.png`: diagnostics.
- `_run_parameters.tsv`, `_celltype_summary.csv`: settings and cell/donor counts.

No residualized H5AD is written for these manifest settings. Retain the CSC inputs alongside backed model objects. Plan disk space for the original data plus subset, preprocessed, and CSC copies.

## Optional step 4: cross-cohort meta-analysis

After completing all relevant contrasts, run a fixed-effect meta-analysis with `metafor::rma(method = "FE")`. For SCZ:

```bash
Rscript scripts/lucida_meta.R \
  --out_dir results/meta --disease SCZ --control CTRL --method FE \
  --input results/CMC_SCZ_deg.parquet \
  --input results/SZBD_SCZ_deg.parquet \
  --input results/MultiomeBrain_SCZ_deg.parquet \
  --nthreads 8 --plot_ncol 6 --plot_dpi 300 --plot_width 20 --plot_height 16
```

Use `.csv.gz` input paths instead when Arrow was unavailable. Repeat with `--disease` and `--input` arguments for these groups:

- BD: `SZBD_BD`, `MultiomeBrain_BD`.
- ASD: `UCLA_ASD`, `Velmeshev_ASD`.
- MDD: `PTSDBrainomics_MDD`, `Turecki_MDD`.
- PTSD has one cohort and is not meta-analyzed.

The script follows [this example](https://gist.github.com/GabrielHoffman/d01e1b778fc0ca9f41a8232fab14478f): for each gene and cell type it fits `rma(yi = logFC, sei = se, method = "FE")` using reported standard errors (or `abs(logFC / t)`), then applies one BH correction across all gene × cell-type rows (`adj.P.Val`). The within-cell-type BH is kept as `adj.P.Val.assay`, and `Q`, `Q.p.value`, and `I2` report heterogeneity between cohorts. Check `n.studies`: genes need not pass the expression filter in every cohort, and a row with `n.studies = 1` is that cohort's own estimate re-tested with a z-test, not a combined result. Outputs include combined DEG and meta-analysis tables, a volcano plot, parameter TSV, and summary markdown under `results/meta/`.

## Scope

Only analysis code and configuration are included. `.gitignore` excludes local data and results. Numerical results may vary across dependency versions and hardware.
