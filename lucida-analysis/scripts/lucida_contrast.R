#!/usr/bin/env Rscript

# lucida_contrast.R
#
# R-only Lucida analysis for one preprocessed cohort H5AD. The preprocessing
# step (contrast subsetting, Scanpy QC metrics, and covariate recoding) lives in
# lucida_preproc.py and must be run before this script.

suppressPackageStartupMessages({
  library(lucida)
  library(GenomicDataStream)
  library(SingleCellExperiment)
  library(S4Vectors)
  library(ggplot2)
  library(dplyr)
})

DelayedArray::setAutoBlockSize(1e9)
options(stringsAsFactors = FALSE)

log_msg <- function(...) {
  cat(sprintf('[%s] ', format(Sys.time(), '%Y-%m-%d %H:%M:%S')), sprintf(...), '\n')
  flush.console()
}

usage <- function() {
  cat(
    'Usage:\n',
    '  Rscript lucida_contrast.R --file <cohort_preproc.h5ad> --cohort <name> --out_dir <dir>\n',
    '    --disease <disease> [--control CTRL] [--formula "~ dx + (1|Donor) + Age + Sex + PMI"]\n',
    '    [--null_formula "~ (1|Donor)"] [--coef dx] [--cell_type_col subclass]\n',
    '    [--random_dx false] [--random_seed auto|integer]\n',
    '    [--nthreads 4] [--chunk_size 1000] [--write_varpart true]\n',
    '    [--write_residualized false] [--plot_ncol 6] [--plot_dpi 220]\n',
    '    [--plot_width 12] [--plot_height 14]\n',
    sep = ''
  )
}

default_options <- list(
  control = 'CTRL',
  coef = 'dx',
  formula = '~ dx + (1|Donor) + Age + Sex + PMI + Chemistry + counts_per_gene + pct_counts_ribo',
  null_formula = '~ (1|Donor)',
  random_dx = 'false',
  random_seed = 'auto',
  nthreads = '8',
  chunk_size = '1000',
  write_varpart = 'true',
  write_residualized = 'false',
  varpart_genes = '',
  plot_ncol = '6',
  plot_dpi = '300',
  plot_width = '16',
  plot_height = '16',
  cell_type_col = 'subclass'
)

parse_cli_args <- function(args, defaults = list()) {
  if (!length(args) || any(args %in% c('-h', '--help'))) {
    usage()
    quit(save = 'no', status = 0)
  }
  if (length(args) %% 2 != 0) stop('Arguments must be provided as --key value pairs')
  keys <- args[c(TRUE, FALSE)]
  vals <- args[c(FALSE, TRUE)]
  if (any(!startsWith(keys, '--'))) stop('All argument names must start with --')
  out <- defaults
  provided <- as.list(vals)
  names(provided) <- sub('^--', '', keys)
  dashed <- names(provided)[grepl('-', names(provided), fixed = TRUE)]
  if (length(dashed)) {
    stop(sprintf('Use underscores in argument names, not hyphens: %s', paste(paste0('--', dashed), collapse = ', ')))
  }
  out[names(provided)] <- provided
  out
}

require_arg <- function(opts, key) {
  val <- opts[[key]]
  if (is.null(val) || identical(val, '')) stop(sprintf('Missing required argument --%s', key))
  val
}

parse_bool <- function(x, default = FALSE) {
  if (is.null(x) || identical(x, '')) return(default)
  val <- tolower(trimws(as.character(x)))
  if (val %in% c('1', 'true', 't', 'yes', 'y', 'on')) return(TRUE)
  if (val %in% c('0', 'false', 'f', 'no', 'n', 'off')) return(FALSE)
  stop(sprintf('Could not parse boolean value: %s', x))
}

normalize_formula <- function(x) {
  x <- trimws(x)
  if (!startsWith(x, '~')) x <- paste('~', x)
  as.formula(x)
}

stable_random_seed <- function(base_seed, ...) {
  labels <- unlist(list(...), recursive = TRUE, use.names = FALSE)
  text <- paste(c(as.character(base_seed), as.character(labels)), collapse = '|')
  bytes <- as.integer(charToRaw(enc2utf8(text)))
  seed <- as.numeric(base_seed) %% 2147483647
  for (b in bytes) {
    seed <- (seed * 131 + b) %% 2147483647
  }
  seed <- as.integer(seed)
  if (is.na(seed) || seed <= 0L) seed <- 1L
  seed
}

generate_time_seed <- function() {
  seed <- as.integer((as.numeric(format(Sys.time(), '%Y%m%d%H%M%S')) + as.integer(Sys.getpid())) %% 2147483646) + 1L
  if (is.na(seed) || seed <= 0L) seed <- 1L
  seed
}

read_h5ad_quiet <- function(file) {
  withCallingHandlers(
    GenomicDataStream::readH5AD(file, raw = TRUE, ondisk = TRUE, verbose = FALSE),
    warning = function(w) {
      if (grepl('could not find dimnames in this h5ad file', conditionMessage(w), fixed = TRUE)) {
        invokeRestart('muffleWarning')
      }
    }
  )
}

write_csv_gz <- function(df, path) {
  con <- gzfile(path, open = 'wt')
  on.exit(close(con), add = TRUE)
  utils::write.csv(df, con, row.names = FALSE)
}

write_table_auto <- function(df, base_path_no_ext) {
  if (requireNamespace('arrow', quietly = TRUE)) {
    out_path <- paste0(base_path_no_ext, '.parquet')
    arrow::write_parquet(df, out_path)
    return(out_path)
  }
  out_path <- paste0(base_path_no_ext, '.csv.gz')
  write_csv_gz(df, out_path)
  out_path
}

write_parameter_tsv <- function(path, kv) {
  df <- data.frame(parameter = names(kv), value = unname(as.character(kv)), stringsAsFactors = FALSE)
  write.table(df, path, sep = '\t', quote = FALSE, row.names = FALSE)
}

apply_random_dx <- function(sce, cohort, disease, control, random_seed) {
  cd <- as.data.frame(colData(sce))
  if (!'dx' %in% colnames(cd)) {
    stop('Preprocessed h5ad is missing required dx column for randomization')
  }
  if (!'Donor' %in% colnames(cd)) {
    stop('Preprocessed h5ad is missing required Donor column for donor-level randomization')
  }

  dx_observed <- suppressWarnings(as.numeric(as.character(cd$dx)))
  if (anyNA(dx_observed)) {
    dx_observed <- suppressWarnings(as.numeric(cd$dx))
  }
  if (anyNA(dx_observed)) {
    stop('Could not coerce dx to numeric before randomization')
  }
  if (!all(dx_observed %in% c(0, 1))) {
    stop('Expected dx to contain only 0/1 values before randomization')
  }

  donor <- as.character(cd$Donor)
  if (anyNA(donor) || any(!nzchar(donor))) {
    stop('Donor column contains missing or empty values before donor-level randomization')
  }

  dx_by_donor <- split(dx_observed, donor, drop = TRUE)
  mixed_donors <- names(dx_by_donor)[vapply(dx_by_donor, function(x) length(unique(x)) != 1L, logical(1))]
  if (length(mixed_donors)) {
    stop(sprintf(
      'Cannot run donor-level randomization because %d donor(s) have mixed dx labels before permutation; examples: %s',
      length(mixed_donors),
      paste(utils::head(mixed_donors, 5), collapse = ', ')
    ))
  }

  donor_ids <- sort(names(dx_by_donor))
  donor_dx_observed <- vapply(dx_by_donor[donor_ids], function(x) unique(x)[[1]], numeric(1))
  n_case_donors <- sum(donor_dx_observed == 1)
  n_control_donors <- length(donor_dx_observed) - n_case_donors
  if (n_case_donors == 0L || n_control_donors == 0L) {
    stop('Donor-level randomization requires at least one disease donor and one control donor')
  }

  draw_seed <- stable_random_seed(random_seed, cohort, disease, control, length(donor_ids), n_case_donors)
  old_seed <- if (exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE)) get('.Random.seed', envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE)) rm('.Random.seed', envir = .GlobalEnv)
    } else {
      assign('.Random.seed', old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(draw_seed)

  donor_dx_randomized <- sample(donor_dx_observed, length(donor_dx_observed), replace = FALSE)
  names(donor_dx_randomized) <- donor_ids
  randomized_dx <- unname(donor_dx_randomized[donor])
  n_random_cases <- sum(randomized_dx == 1)
  n_random_controls <- length(randomized_dx) - n_random_cases

  colData(sce)$dx_observed <- dx_observed
  colData(sce)$dx <- randomized_dx
  if ('Disease' %in% colnames(cd)) {
    colData(sce)$Disease_observed <- as.character(cd$Disease)
    colData(sce)$Disease <- ifelse(randomized_dx == 1, disease, control)
  }
  colData(sce)$dx_randomized <- TRUE
  colData(sce)$dx_random_seed <- draw_seed
  colData(sce)$dx_random_base_seed <- as.integer(random_seed)
  metadata(sce)$random_dx <- list(
    enabled = TRUE,
    mode = 'donor_label_permutation_in_memory_before_lucida_fit',
    unit = 'Donor',
    base_seed = as.integer(random_seed),
    draw_seed = draw_seed,
    cohort = cohort,
    disease = disease,
    control = control,
    n_cells = ncol(sce),
    n_cases = n_random_cases,
    n_controls = n_random_controls,
    n_donors = length(donor_ids),
    n_case_donors = n_case_donors,
    n_control_donors = n_control_donors
  )
  log_msg(
    'Randomized dx for %s at donor level before lucida fit: %s disease donors, %s control donors; resulting cells: %s cases, %s controls; draw_seed=%s',
    cohort,
    format(n_case_donors, big.mark = ','),
    format(n_control_donors, big.mark = ','),
    format(n_random_cases, big.mark = ','),
    format(n_random_controls, big.mark = ','),
    draw_seed
  )
  sce
}

prepare_sce <- function(file, disease, control, formula, null_formula = NULL, cell_type_col = 'subclass',
                        random_dx = FALSE, random_seed = 104729, cohort = '') {
  log_msg('Loading preprocessed h5ad %s', file)
  rhdf5::h5closeAll()
  sce <- read_h5ad_quiet(file)
  cd <- as.data.frame(colData(sce))

  if (!cell_type_col %in% colnames(cd)) {
    stop(sprintf('Preprocessed h5ad is missing required cell-type column: %s', cell_type_col))
  }
  formula_vars <- unique(c(all.vars(formula), if (!is.null(null_formula)) all.vars(null_formula) else character()))
  missing_formula_vars <- setdiff(formula_vars, colnames(cd))
  if (length(missing_formula_vars)) {
    stop(sprintf('Preprocessed h5ad is missing formula covariates: %s', paste(missing_formula_vars, collapse = ', ')))
  }

  if (random_dx) {
    sce <- apply_random_dx(sce, cohort = cohort, disease = disease, control = control, random_seed = random_seed)
  } else {
    metadata(sce)$random_dx <- list(enabled = FALSE)
  }

  attr(sce, 'cell_type_col') <- cell_type_col
  attr(sce, 'subject_id_col') <- if ('Donor' %in% colnames(cd)) 'Donor' else ''
  random_meta <- metadata(sce)$random_dx
  attr(sce, 'audit') <- list(
    n_cells_input = nrow(cd),
    n_cells_analyzed = ncol(sce),
    cell_type_col = cell_type_col,
    random_dx = isTRUE(random_meta$enabled),
    randomization_unit = if (is.null(random_meta$unit)) NA_character_ else random_meta$unit,
    random_base_seed = if (is.null(random_meta$base_seed)) NA_integer_ else random_meta$base_seed,
    random_draw_seed = if (is.null(random_meta$draw_seed)) NA_integer_ else random_meta$draw_seed,
    random_n_cases = if (is.null(random_meta$n_cases)) NA_integer_ else random_meta$n_cases,
    random_n_controls = if (is.null(random_meta$n_controls)) NA_integer_ else random_meta$n_controls,
    random_n_case_donors = if (is.null(random_meta$n_case_donors)) NA_integer_ else random_meta$n_case_donors,
    random_n_control_donors = if (is.null(random_meta$n_control_donors)) NA_integer_ else random_meta$n_control_donors
  )

  log_msg('%s ready for %s vs %s: %s genes x %s cells',
          basename(file), disease, control,
          format(nrow(sce), big.mark = ','), format(ncol(sce), big.mark = ','))
  sce
}

make_celltype_summary <- function(sce, cohort) {
  cell_type_col <- attr(sce, 'cell_type_col')
  subject_id_col <- attr(sce, 'subject_id_col')
  if (!nzchar(subject_id_col) || !'Disease' %in% colnames(colData(sce))) {
    stop('Preprocessed h5ad must contain obs/Donor and obs/Disease to write cell-type summary')
  }
  as.data.frame(colData(sce)) %>%
    dplyr::group_by(assay = .data[[cell_type_col]], Disease) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      n_donors = dplyr::n_distinct(.data[[subject_id_col]]),
      .groups = 'drop'
    ) %>%
    dplyr::mutate(cohort = cohort) %>%
    dplyr::select(cohort, assay, Disease, dplyr::everything()) %>%
    dplyr::arrange(assay, Disease)
}

render_diag_plots <- function(fit, out_dir, cohort, coef, plot_cfg) {
  ggsave(file.path(out_dir, paste0(cohort, '_dispEsts.png')), plotDispEsts(fit, ncol = plot_cfg$ncol), width = plot_cfg$width, height = plot_cfg$height, dpi = plot_cfg$dpi)
  ggsave(file.path(out_dir, paste0(cohort, '_MA.png')), plotMA(fit, coef = coef, ymax = 6, ncol = plot_cfg$ncol), width = plot_cfg$width, height = plot_cfg$height, dpi = plot_cfg$dpi)
  ggsave(file.path(out_dir, paste0(cohort, '_volcano.png')), plotVolcano(fit, coef = coef, ncol = plot_cfg$ncol), width = plot_cfg$width, height = plot_cfg$height, dpi = plot_cfg$dpi)
}

render_varpart_outputs <- function(sce, fit, out_dir, cohort, null_formula, nthreads, chunk_size, plot_cfg, varpart_genes) {
  cell_type_col <- attr(sce, 'cell_type_col')
  log_msg('Running null lucida fit for variance partition: %s', deparse(null_formula))
  fit_null <- lucida(sce, formula = null_formula, cluster_id = cell_type_col, verbose = TRUE, nthreads = nthreads, nReaders = 1, chunkSize = chunk_size)
  fit_null_path <- file.path(out_dir, paste0(cohort, '_lucida_fit_null.rds'))
  saveRDS(fit_null, fit_null_path)
  log_msg('Saved null fit for variance partition: %s', fit_null_path)

  vp <- fitVarPart(fit, fit_null)
  vp_path <- file.path(out_dir, paste0(cohort, '_varpart.rds'))
  saveRDS(vp, vp_path)
  log_msg('Saved variance-partition results: %s', vp_path)

  ggsave(file.path(out_dir, paste0(cohort, '_varpart_violin.png')), plotVarPart(sortCols(vp), ncol = plot_cfg$ncol), width = plot_cfg$width, height = plot_cfg$height, dpi = plot_cfg$dpi)
  log_msg('Saved variance-partition violin plot for %s', cohort)

  if (length(varpart_genes)) {
    vp_bar <- dplyr::filter(vp, ID %in% varpart_genes)
    if (nrow(vp_bar)) {
      ggsave(file.path(out_dir, paste0(cohort, '_varpart_bars.png')), plotPercentBars(sortCols(vp_bar), ncol = plot_cfg$ncol), width = plot_cfg$width, height = plot_cfg$height, dpi = plot_cfg$dpi)
      log_msg('Saved variance-partition bar plot for %s (%s genes)', cohort, nrow(vp_bar))
    } else {
      log_msg('Skipped variance-partition bar plot: requested genes not found in variance-partition table')
    }
  } else {
    log_msg('Skipping variance-partition bar plot: no --varpart_genes supplied')
  }
}

write_residualized_h5ad <- function(sce, input_file, out_dir, cohort, nthreads, chunk_size) {
  if (!requireNamespace('zellkonverter', quietly = TRUE)) {
    stop('Package zellkonverter is required for --write_residualized true')
  }

  log_msg('Running intercept-only all-cell lucida fit for residualized h5ad export')
  fit_all <- lucida(sce, formula = ~ 1, cluster_id = NULL, verbose = TRUE, nthreads = nthreads, nReaders = 1, chunkSize = chunk_size)
  res_all <- residuals(fit_all, sce)

  if (identical(dim(res_all), c(nrow(sce), ncol(sce)))) {
    residualized <- res_all
  } else if (identical(dim(res_all), c(ncol(sce), nrow(sce)))) {
    residualized <- t(res_all)
  } else {
    stop(sprintf(
      'Residualized matrix dimensions %s do not match SCE dimensions %s genes x %s cells',
      paste(dim(res_all), collapse = ' x '),
      format(nrow(sce), big.mark = ','),
      format(ncol(sce), big.mark = ',')
    ))
  }

  resid_path <- file.path(out_dir, paste0(cohort, '_resid.h5ad'))
  if (file.exists(resid_path)) unlink(resid_path)

  resid_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(residualized = residualized),
    colData = SummarizedExperiment::colData(sce),
    rowData = SummarizedExperiment::rowData(sce)
  )

  log_msg('Writing residualized-only h5ad to X with zellkonverter::writeH5AD (omitting original X and raw/X)')
  zellkonverter::writeH5AD(
    resid_sce,
    resid_path,
    X_name = 'residualized',
    compression = 'gzip',
    version = '0.10.9',
    reducedDims = FALSE,
    metadata = FALSE,
    colPairs = FALSE,
    rowPairs = FALSE,
    verbose = TRUE
  )
  log_msg('Saved residualized-only h5ad: %s', resid_path)
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE), defaults = default_options)

file <- require_arg(args, 'file')
out_dir <- require_arg(args, 'out_dir')
disease <- require_arg(args, 'disease')
cohort <- args[['cohort']]
if (is.null(cohort) || identical(cohort, '')) cohort <- sub('\\.h5ad$', '', basename(file))
control <- args[['control']]
coef <- args[['coef']]
formula <- normalize_formula(args[['formula']])
null_formula <- normalize_formula(args[['null_formula']])
nthreads <- as.integer(args[['nthreads']])
chunk_size <- as.integer(args[['chunk_size']])
write_varpart <- parse_bool(args[['write_varpart']])
write_residualized <- parse_bool(args[['write_residualized']])
random_dx <- parse_bool(args[['random_dx']])
random_seed_arg <- trimws(as.character(args[['random_seed']]))
if (tolower(random_seed_arg) == 'auto') {
  random_seed <- if (random_dx) generate_time_seed() else 104729L
} else {
  random_seed <- as.integer(random_seed_arg)
  if (is.na(random_seed)) stop(sprintf('Could not parse integer --random_seed: %s', args[['random_seed']]))
}
varpart_genes <- unique(Filter(nzchar, trimws(strsplit(args[['varpart_genes']], ',', fixed = TRUE)[[1]])))

plot_cfg <- list(
  ncol = as.integer(args[['plot_ncol']]),
  dpi = as.integer(args[['plot_dpi']]),
  width = as.numeric(args[['plot_width']]),
  height = as.numeric(args[['plot_height']])
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (write_varpart && is.null(null_formula)) stop('Missing required --null_formula when --write_varpart true')

sce <- prepare_sce(
  file = file,
  disease = disease,
  control = control,
  formula = formula,
  null_formula = null_formula,
  cell_type_col = args[['cell_type_col']],
  random_dx = random_dx,
  random_seed = random_seed,
  cohort = cohort
)
audit <- attr(sce, 'audit')

param_path <- file.path(out_dir, paste0(cohort, '_run_parameters.tsv'))
write_parameter_tsv(param_path, c(
  file = file,
  cohort = cohort,
  out_dir = out_dir,
  disease = disease,
  control = control,
  formula = deparse(formula),
  null_formula = deparse(null_formula),
  coef = coef,
  random_dx = random_dx,
  randomization_unit = audit$randomization_unit,
  random_seed = random_seed,
  random_draw_seed = audit$random_draw_seed,
  random_n_cases = audit$random_n_cases,
  random_n_controls = audit$random_n_controls,
  random_n_case_donors = audit$random_n_case_donors,
  random_n_control_donors = audit$random_n_control_donors,
  nthreads = nthreads,
  chunk_size = chunk_size,
  plot_ncol = plot_cfg$ncol,
  plot_dpi = plot_cfg$dpi,
  plot_width = plot_cfg$width,
  plot_height = plot_cfg$height,
  write_varpart = write_varpart,
  write_residualized = write_residualized,
  varpart_genes = if (length(varpart_genes)) paste(varpart_genes, collapse = ',') else '',
  cells_input = audit$n_cells_input,
  cells_analyzed = audit$n_cells_analyzed,
  cell_type_col = audit$cell_type_col,
  parquet_available = requireNamespace('arrow', quietly = TRUE)
))
log_msg('Saved run parameters: %s', param_path)

summary_path <- file.path(out_dir, paste0(cohort, '_celltype_summary.csv'))
write.csv(make_celltype_summary(sce, cohort), summary_path, row.names = FALSE)
log_msg('Saved cell-type summary: %s', summary_path)

cell_type_col <- attr(sce, 'cell_type_col')
log_msg('Skipping SCE RDS save because preprocessed h5ad is out-of-memory backed')

log_msg('Running lucida for %s', cohort)
fit <- lucida(sce, formula = formula, cluster_id = cell_type_col, verbose = TRUE, nthreads = nthreads, nReaders = 1, chunkSize = chunk_size)

fit_path <- file.path(out_dir, paste0(cohort, '_lucida_fit.rds'))
saveRDS(fit, fit_path)
log_msg('Saved fit: %s', fit_path)

res <- results(fit, coef = coef)
if ('cluster_id' %in% colnames(res)) res <- dplyr::rename(res, assay = cluster_id)
res$Dataset <- cohort

deg_path <- write_table_auto(as.data.frame(res), file.path(out_dir, paste0(cohort, '_deg')))
log_msg('Saved DEG table: %s (%s rows)', deg_path, format(nrow(res), big.mark = ','))
if (!requireNamespace('arrow', quietly = TRUE)) log_msg('arrow not available in lucida env; used csv.gz fallback instead of parquet')

render_diag_plots(fit, out_dir, cohort, coef, plot_cfg)
log_msg('Saved diagnostic PNGs for %s', cohort)

if (write_varpart) {
  render_varpart_outputs(sce, fit, out_dir, cohort, null_formula, nthreads, chunk_size, plot_cfg, varpart_genes)
} else {
  log_msg('Skipping variance-partition outputs (enable with --write_varpart true)')
}

if (write_residualized) {
  write_residualized_h5ad(sce, file, out_dir, cohort, nthreads, chunk_size)
} else {
  log_msg('Skipping residualized h5ad export (enable with --write_residualized true)')
}

rhdf5::h5closeAll()
