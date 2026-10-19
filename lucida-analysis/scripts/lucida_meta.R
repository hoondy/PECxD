#!/usr/bin/env Rscript

# lucida_meta.R
#
# Run disease-level fixed-effect meta-analysis across DEG tables that were
# already produced by lucida_contrast.R for individual cohorts, following
# https://gist.github.com/GabrielHoffman/d01e1b778fc0ca9f41a8232fab14478f :
# metafor::rma(yi = logFC, sei = se, method = 'FE') for each gene x cell type,
# then one FDR correction across all gene x cell-type rows.
#
# This script covers only cross-cohort aggregation. It does NOT run lucida itself.
#
# Example:
#   Rscript lucida_meta.R \
#     --out_dir out/SCZ_lucida_strict \
#     --disease SCZ \
#     --control CTRL \
#     --method FE \
#     --input out/SCZ_lucida_strict/CMC_finalmeta_deg.csv.gz \
#     --input out/SCZ_lucida_strict/SZBD_finalmeta_deg.csv.gz
#
# Notes:
# - The script prefers parquet output when the optional R package `arrow` is
#   available. Otherwise it falls back to compressed CSV (`.csv.gz`).
# - Standard errors are taken from `se` / `std.error` when present, or derived as
#   abs(logFC / t) if only the Lucida t-statistic is available.
# - adj.P.Val is BH (`p.adjust(..., 'fdr')`) across all gene x cell-type rows,
#   as in the gist. adj.P.Val.assay keeps the previous within-cell-type BH for
#   comparison.
# - A gene x cell type observed in only one cohort is still fitted (rma with
#   k = 1 returns that cohort's estimate with a Wald z-test); n.studies and
#   datasets record this so single-cohort rows can be filtered downstream.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(metafor)
})

options(stringsAsFactors = FALSE)

log_msg <- function(...) {
  cat(sprintf('[%s] ', format(Sys.time(), '%Y-%m-%d %H:%M:%S')), sprintf(...), '\n')
  flush.console()
}

usage <- function() {
  cat(
    'Usage:\n',
    '  Rscript lucida_meta.R --out_dir <dir> --disease <disease> [--control CTRL] [--method FE]\n',
    '    --input <cohort1_deg.{csv.gz|parquet}> --input <cohort2_deg.{csv.gz|parquet}> [--input ...]\n',
    '    [--nthreads 8] [--plot_ncol 6] [--plot_dpi 300] [--plot_width 16] [--plot_height 16]\n',
    sep = ''
  )
}

default_options <- list(
  control = 'CTRL',
  method = 'FE',
  nthreads = '8',
  plot_ncol = '6',
  plot_dpi = '300',
  plot_width = '16',
  plot_height = '16'
)

parse_cli_args <- function(args, defaults = list()) {
  if (!length(args) || any(args %in% c('-h', '--help'))) {
    usage()
    quit(save = 'no', status = 0)
  }

  out <- defaults
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, '--')) {
      stop('All argument names must start with --')
    }
    if (i == length(args)) {
      stop(sprintf('Missing value for %s', key))
    }
    val <- args[[i + 1L]]
    key <- gsub('-', '_', sub('^--', '', key), fixed = TRUE)
    if (identical(key, 'input')) {
      out[[key]] <- c(out[[key]], val)
    } else {
      out[[key]] <- val
    }
    i <- i + 2L
  }
  out
}

require_arg <- function(opts, key) {
  val <- opts[[key]]
  if (is.null(val) || identical(val, '') || (is.atomic(val) && !length(val))) {
    stop(sprintf('Missing required argument --%s', key))
  }
  val
}

write_csv_gz <- function(df, path) {
  con <- gzfile(path, open = 'wt')
  on.exit(close(con), add = TRUE)
  utils::write.csv(df, con, row.names = FALSE)
}

read_csv_gz <- function(path) {
  con <- gzfile(path, open = 'rt')
  on.exit(close(con), add = TRUE)
  utils::read.csv(con, check.names = FALSE)
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

read_table_auto <- function(path) {
  if (grepl('\\.parquet$', path, ignore.case = TRUE)) {
    if (!requireNamespace('arrow', quietly = TRUE)) {
      stop(sprintf('Cannot read parquet without arrow: %s', path))
    }
    return(as.data.frame(arrow::read_parquet(path)))
  }
  if (grepl('\\.csv\\.gz$', path, ignore.case = TRUE)) return(read_csv_gz(path))
  if (grepl('\\.csv$', path, ignore.case = TRUE)) return(utils::read.csv(path, check.names = FALSE))
  stop(sprintf('Unsupported input file type: %s', path))
}

write_parameter_tsv <- function(path, kv) {
  df <- data.frame(parameter = names(kv), value = unname(as.character(kv)), stringsAsFactors = FALSE)
  write.table(df, path, sep = '\t', quote = FALSE, row.names = FALSE)
}

plotMetaVolcano <- function(meta, nGenes = 5, size = 12, ymax = 300, xmax = 6, cutoff = 0.05, ncol = 3, assays = unique(meta$assay), ...) {
  ID <- assay <- estimate <- p.value <- adj.P.Val <- isSignif <- y <- NULL
  required_cols <- c('ID', 'assay', 'estimate', 'std.error', 'p.value', 'adj.P.Val')
  missing <- setdiff(required_cols, colnames(meta))
  if (length(missing)) {
    stop(sprintf('Meta-analysis table is missing required volcano columns: %s', paste(missing, collapse = ', ')))
  }
  if (!all(assays %in% meta$assay)) {
    txt <- paste(assays[!assays %in% meta$assay], collapse = ', ')
    stop('assays not valid: ', txt)
  }
  if (!requireNamespace('ggrepel', quietly = TRUE)) {
    stop('Package ggrepel is required to render meta-analysis volcano plots')
  }

  df <- meta %>%
    dplyr::filter(assay %in% assays) %>%
    dplyr::mutate(
      assay = factor(assay, assays),
      estimate = as.numeric(estimate),
      p.value = as.numeric(p.value),
      adj.P.Val = as.numeric(adj.P.Val)
    ) %>%
    dplyr::filter(
      is.finite(estimate),
      is.finite(p.value),
      p.value >= 0,
      is.finite(adj.P.Val)
    ) %>%
    dplyr::mutate(estimate = pmin(xmax, estimate)) %>%
    dplyr::mutate(estimate = pmax(-xmax, estimate)) %>%
    dplyr::mutate(y = pmin(-log10(pmax(p.value, .Machine$double.xmin)), ymax)) %>%
    dplyr::mutate(isSignif = factor(c('no', 'yes')[(adj.P.Val < cutoff) + 1], levels = c('no', 'yes'))) %>%
    dplyr::arrange(dplyr::desc(p.value))

  if (!nrow(df)) {
    stop('No finite rows are available for meta-analysis volcano plotting')
  }

  df2 <- df %>%
    dplyr::filter(!is.na(ID), nzchar(as.character(ID))) %>%
    dplyr::group_by(assay) %>%
    dplyr::slice_tail(n = nGenes)
  ymax <- max(df$y, na.rm = TRUE)

  fig <- ggplot(df, aes(estimate, y, color = isSignif)) +
    geom_point() +
    theme_classic(size) +
    theme(
      aspect.ratio = 1,
      legend.position = 'none',
      strip.background = element_rect('grey95'),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_color_manual(values = c('grey', 'darkred')) +
    ggrepel::geom_text_repel(
      data = df2,
      aes(estimate, y, label = ID),
      segment.size = 0.5,
      segment.color = 'black',
      color = 'black',
      force = 1,
      nudge_x = 0.005,
      nudge_y = 0.5
    ) +
    facet_wrap(~assay, ncol = ncol, ...) +
    xlab(bquote(meta*"-analysis" ~ effect ~ estimate)) +
    ylab(bquote(-log[10] ~ P))

  argsList <- match.call(expand.dots = FALSE)$...
  ags <- names(argsList)
  if (length(ags) == 0 || !('scales' %in% ags)) {
    fig <- fig + scale_y_continuous(limits = c(0, max(ymax * 1.02, ymax + 0.5)), expand = c(0, 1))
  }
  fig
}

derive_standard_error <- function(tab) {
  if ('se' %in% colnames(tab)) {
    se <- as.numeric(tab$se)
    se[!is.finite(se) | se <= 0] <- NA_real_
    return(se)
  }
  if ('std.error' %in% colnames(tab)) {
    se <- as.numeric(tab[['std.error']])
    se[!is.finite(se) | se <= 0] <- NA_real_
    return(se)
  }
  stat_col <- NULL
  if ('t' %in% colnames(tab)) stat_col <- 't'
  if (is.null(stat_col) && 'statistic' %in% colnames(tab)) stat_col <- 'statistic'
  if (is.null(stat_col) || !('logFC' %in% colnames(tab))) {
    stop('Cannot derive standard error: expected se/std.error or logFC plus t/statistic')
  }
  stat <- as.numeric(tab[[stat_col]])
  logfc <- as.numeric(tab[['logFC']])
  se <- abs(logfc / stat)
  se[!is.finite(se) | se <= 0] <- NA_real_
  se
}

# Fit metafor::rma(method = 'FE') to one gene x cell type. Q / I2 are the
# heterogeneity test and statistic across cohorts (Q = 0, I2 = 0 when k = 1).
rma_fe_one <- function(logFC, se) {
  fit <- metafor::rma(yi = logFC, sei = se, method = 'FE')
  c(
    estimate = as.numeric(fit$beta),
    std.error = fit$se,
    statistic = fit$zval,
    p.value = fit$pval,
    Q = fit$QE,
    Q.p.value = fit$QEp,
    I2 = fit$I2,
    n.studies = fit$k
  )
}

fixed_effect_meta <- function(tab, nthreads = 1L) {
  tab$std.error.meta <- derive_standard_error(tab)
  work <- tab %>%
    dplyr::filter(is.finite(logFC), is.finite(std.error.meta), std.error.meta > 0) %>%
    dplyr::arrange(assay, ID)

  # One rma() per gene x cell type, chunked across forked workers.
  key <- paste(work$assay, work$ID, sep = '\r')
  rows_by_group <- split(seq_len(nrow(work)), factor(key, levels = unique(key)))
  chunks <- split(seq_along(rows_by_group), cut(seq_along(rows_by_group), max(1L, nthreads), labels = FALSE))
  fits <- parallel::mclapply(chunks, function(idx) {
    t(vapply(rows_by_group[idx], function(r) rma_fe_one(work$logFC[r], work$std.error.meta[r]), numeric(8)))
  }, mc.cores = nthreads)
  stats <- as.data.frame(do.call(rbind, unname(fits)))

  first <- vapply(rows_by_group, `[`, integer(1), 1L)
  meta <- data.frame(ID = work$ID[first], assay = work$assay[first], stats, row.names = NULL)
  meta$n.studies <- as.integer(meta$n.studies)
  meta$method <- 'FE'
  meta$datasets <- vapply(rows_by_group, function(r) paste(sort(unique(work$Dataset[r])), collapse = ','), character(1))

  meta %>%
    dplyr::mutate(adj.P.Val = p.adjust(p.value, method = 'fdr')) %>%
    dplyr::group_by(assay) %>%
    dplyr::mutate(adj.P.Val.assay = p.adjust(p.value, method = 'BH')) %>%
    dplyr::ungroup() %>%
    dplyr::select(ID, assay, estimate, std.error, statistic, p.value, Q, Q.p.value, I2, n.studies, method, datasets, adj.P.Val, adj.P.Val.assay) %>%
    dplyr::arrange(assay, p.value, ID) %>%
    as.data.frame()
}

write_summary <- function(meta, out_file, disease, control, method, input_files) {
  sig_counts <- meta %>%
    dplyr::group_by(assay) %>%
    dplyr::summarise(
      n_genes = dplyr::n(),
      n_sig = sum(adj.P.Val < 0.05, na.rm = TRUE),
      best_fdr = min(adj.P.Val, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::arrange(dplyr::desc(n_sig), best_fdr)

  lines <- c(
    sprintf('# %s vs %s meta-analysis summary', disease, control),
    '',
    sprintf('Method: metafor::rma(method = "%s") by gene x cell type; FDR (BH) across all gene x cell-type rows', method),
    sprintf('Input DEG files (%d):', length(input_files)),
    paste0('- ', input_files),
    '',
    sprintf('Single-cohort rows (n.studies = 1): %d of %d, %d with FDR < 0.05',
      sum(meta$n.studies == 1), nrow(meta), sum(meta$n.studies == 1 & meta$adj.P.Val < 0.05, na.rm = TRUE)),
    '',
    '## Significant meta-analysis rows (FDR < 0.05) by cell type',
    ''
  )
  for (i in seq_len(nrow(sig_counts))) {
    row <- sig_counts[i, ]
    lines <- c(lines, sprintf('- %s: %d / %d significant (best FDR %.3e)', row$assay, row$n_sig, row$n_genes, row$best_fdr))
  }

  lines <- c(lines, '', '## Top 5 hits per cell type', '')
  for (assay_name in unique(meta$assay)) {
    top <- meta %>% dplyr::filter(assay == assay_name) %>% dplyr::slice_min(order_by = adj.P.Val, n = 5)
    lines <- c(lines, sprintf('### %s', assay_name), '', '| Gene | beta | z | p | FDR | studies |', '| --- | ---: | ---: | ---: | ---: | ---: |')
    for (i in seq_len(nrow(top))) {
      row <- top[i, ]
      lines <- c(lines, sprintf('| %s | %.4f | %.3f | %.3e | %.3e | %d |', row$ID, row$estimate, row$statistic, row$p.value, row$adj.P.Val, row$n.studies))
    }
    lines <- c(lines, '')
  }

  writeLines(lines, out_file)
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE), defaults = default_options)
out_dir <- require_arg(args, 'out_dir')
disease <- require_arg(args, 'disease')
input_files <- require_arg(args, 'input')
control <- args[['control']]
method <- args[['method']]
nthreads <- as.integer(args[['nthreads']])
plot_cfg <- list(
  ncol = as.integer(args[['plot_ncol']]),
  dpi = as.integer(args[['plot_dpi']]),
  width = as.numeric(args[['plot_width']]),
  height = as.numeric(args[['plot_height']])
)

if (!identical(method, 'FE')) {
  stop('This script currently implements metafor fixed-effect meta-analysis only (method FE)')
}
if (length(input_files) < 2) {
  stop('Meta-analysis requires at least two --input DEG files')
}
if (any(!file.exists(input_files))) {
  missing <- input_files[!file.exists(input_files)]
  stop(sprintf('Missing input files: %s', paste(missing, collapse = ', ')))
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
param_path <- file.path(out_dir, paste0(disease, '_meta_parameters.tsv'))
write_parameter_tsv(param_path, c(
  out_dir = out_dir,
  disease = disease,
  control = control,
  method = method,
  meta_engine = 'metafor_rma',
  metafor_version = as.character(utils::packageVersion('metafor')),
  fdr_scope = 'all_rows',
  nthreads = nthreads,
  n_inputs = length(input_files),
  input_files = paste(input_files, collapse = ';'),
  plot_ncol = plot_cfg$ncol,
  plot_dpi = plot_cfg$dpi,
  plot_width = plot_cfg$width,
  plot_height = plot_cfg$height,
  parquet_available = requireNamespace('arrow', quietly = TRUE)
))
log_msg('Saved meta parameters: %s', param_path)

# Load DEG tables produced independently by lucida_contrast.R.
tabs <- lapply(input_files, read_table_auto)
required_cols <- c('ID', 'assay', 'logFC', 'Dataset')
for (i in seq_along(tabs)) {
  missing <- setdiff(required_cols, colnames(tabs[[i]]))
  if (length(missing)) {
    stop(sprintf('Input %s is missing required columns: %s', input_files[[i]], paste(missing, collapse = ', ')))
  }
}

combined <- dplyr::bind_rows(tabs)
combined_path <- write_table_auto(combined, file.path(out_dir, paste0(disease, '_combined_deg')))
log_msg('Saved combined DEG table: %s (%s rows)', combined_path, format(nrow(combined), big.mark = ','))

meta <- fixed_effect_meta(combined, nthreads = nthreads)
meta_path <- write_table_auto(meta, file.path(out_dir, paste0(disease, '_meta')))
log_msg('Saved meta-analysis table: %s (%s rows)', meta_path, format(nrow(meta), big.mark = ','))
if (!requireNamespace('arrow', quietly = TRUE)) {
  log_msg('arrow not available in lucida env; used csv.gz fallback instead of parquet')
}

volcano_path <- file.path(out_dir, paste0(disease, '_meta_volcano.png'))
ggsave(
  volcano_path,
  plotMetaVolcano(meta, ncol = plot_cfg$ncol),
  width = plot_cfg$width,
  height = plot_cfg$height,
  dpi = plot_cfg$dpi
)
log_msg('Saved meta-analysis volcano plot: %s', volcano_path)

summary_path <- file.path(out_dir, paste0(disease, '_analysis_summary.md'))
write_summary(meta, summary_path, disease, control, method, input_files)
log_msg('Saved summary markdown: %s', summary_path)
