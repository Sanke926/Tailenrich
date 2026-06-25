#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0L) {
    stop("Cannot determine script path from commandArgs().")
  }
  normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
}

script_path <- get_script_path()
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

input_root <- normalizePath(Sys.getenv("TAILENRICH_INPUT_ROOT", file.path(repo_root, "sample_input")), winslash = "/", mustWork = TRUE)
output_root <- normalizePath(Sys.getenv("TAILENRICH_OUTPUT_ROOT", file.path(repo_root, "sample_output_rerun")), winslash = "/", mustWork = FALSE)
tailenrich_script <- normalizePath(
  Sys.getenv("TAILENRICH_SCRIPT", file.path(repo_root, "src", "tailEnrich.R")),
  winslash = "/",
  mustWork = TRUE
)

dataset_dirs_env <- trimws(Sys.getenv("TAILENRICH_DATASETS", unset = ""))
if (nzchar(dataset_dirs_env)) {
  dataset_dirs <- trimws(strsplit(dataset_dirs_env, ",", fixed = TRUE)[[1]])
  dataset_dirs <- dataset_dirs[nzchar(dataset_dirs)]
} else {
  dataset_dirs <- basename(list.dirs(input_root, full.names = TRUE, recursive = FALSE))
}

seed <- suppressWarnings(as.integer(Sys.getenv("TAILENRICH_SEED", "1")))
n_cores <- suppressWarnings(as.integer(Sys.getenv("TAILENRICH_N_CORES", "1")))
n_perm_env <- trimws(Sys.getenv("TAILENRICH_N_PERM", unset = ""))
n_perm <- if (nzchar(n_perm_env)) suppressWarnings(as.integer(n_perm_env)) else NA_integer_
fdr_cutoff <- suppressWarnings(as.numeric(Sys.getenv("TAILENRICH_FDR_CUTOFF", "0.05")))
log2fc_cutoff <- suppressWarnings(as.numeric(Sys.getenv("TAILENRICH_LOG2FC_CUTOFF", as.character(log2(1.5)))))

expr_filename <- "gene_TPM_by_salmon_covAdjusted.csv"
group_filename <- "used_samples_group.tsv"

stop_if_not_exists <- function(path, desc = "file") {
  if (!file.exists(path)) stop("Cannot find ", desc, ": ", path)
}

read_expr_matrix <- function(expr_file) {
  df <- read.csv(expr_file, check.names = FALSE)
  if (ncol(df) < 2L) stop("Expression file must have at least 2 columns: ", expr_file)

  gene_ids <- as.character(df[[1]])
  if (anyDuplicated(gene_ids)) stop("Duplicated gene IDs in expression matrix: ", expr_file)

  mat <- as.matrix(df[, -1, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "double")
  rownames(mat) <- gene_ids

  bad <- which(!is.finite(mat), arr.ind = TRUE)
  if (nrow(bad) > 0L) stop("Expression matrix contains NA/NaN/Inf values: ", expr_file)

  mat
}

read_group_table <- function(group_file) {
  df <- read.delim(group_file, sep = "\t", check.names = FALSE)
  need_cols <- c("sampleID", "group")
  miss <- setdiff(need_cols, colnames(df))
  if (length(miss) > 0L) {
    stop("Missing required columns in ", group_file, ": ", paste(miss, collapse = ", "))
  }

  df$sampleID <- as.character(df$sampleID)
  df$group <- suppressWarnings(as.numeric(df$group))
  if (anyNA(df$sampleID) || any(df$sampleID == "")) stop("Empty sampleID found in: ", group_file)
  if (anyNA(df$group)) stop("Non-numeric group values found in: ", group_file)
  if (!all(df$group %in% c(-1, 1))) stop("Group column must be coded as -1/1 in: ", group_file)
  if (anyDuplicated(df$sampleID)) stop("Duplicated sampleID in: ", group_file)

  df
}

align_expr_and_group <- function(expr_mat, group_df, dataset_name) {
  expr_samples <- colnames(expr_mat)
  idx <- match(expr_samples, group_df$sampleID)
  if (anyNA(idx)) {
    stop("Samples in expression matrix are missing from group table for dataset: ", dataset_name)
  }

  group_df2 <- group_df[idx, , drop = FALSE]
  if (!identical(as.character(group_df2$sampleID), as.character(expr_samples))) {
    stop("Sample order mismatch after alignment for dataset: ", dataset_name)
  }

  y <- as.numeric(group_df2$group)
  if (length(y) != ncol(expr_mat)) stop("length(group) != ncol(expr_mat) for dataset: ", dataset_name)

  list(expr = expr_mat, y = y)
}

resolve_ttailenrich_entry <- function(env, script_path) {
  preferred <- c("tailenrich", "tailEnrich", "tail_enrich")
  for (nm in preferred) {
    if (exists(nm, envir = env, inherits = FALSE) && is.function(get(nm, envir = env, inherits = FALSE))) {
      return(nm)
    }
  }

  hits <- ls(env, all.names = TRUE)
  hits_fun <- hits[vapply(hits, function(nm) is.function(get(nm, envir = env, inherits = FALSE)), logical(1))]
  hits_tail <- hits_fun[grepl("tail", hits_fun, ignore.case = TRUE)]
  if (length(hits_tail) >= 1L) return(hits_tail[1])

  stop("Cannot resolve TailEnrich entry function from script: ", script_path)
}

run_tailenrich_one <- function(
  counts,
  y,
  te_fun,
  seed = 1L,
  n_cores = 1L,
  n_perm = NA_integer_,
  fdr_cutoff = 0.05,
  log2fc_cutoff = log2(1.5)
) {
  args <- list(X = counts, y = y, seed = seed, n_cores = n_cores)
  if (!is.na(n_perm) && "n_perm" %in% names(formals(te_fun))) {
    args$n_perm <- n_perm
  }

  te_res <- do.call(te_fun, args)
  te_res <- as.data.frame(te_res, stringsAsFactors = FALSE, check.names = FALSE)

  genes <- rownames(counts)
  if (!is.null(rownames(te_res)) && all(genes %in% rownames(te_res))) {
    te_res <- te_res[genes, , drop = FALSE]
  } else if ("gene" %in% colnames(te_res) && all(genes %in% te_res$gene)) {
    te_res <- te_res[match(genes, te_res$gene), , drop = FALSE]
  } else if (nrow(te_res) != nrow(counts)) {
    stop("TailEnrich output cannot be aligned to input genes.")
  }

  cn0 <- colnames(te_res)
  cn <- tolower(cn0)
  find_col <- function(candidates) {
    cand <- tolower(candidates)
    hit <- match(cand, cn, nomatch = 0)
    if (all(hit == 0)) return(NULL)
    cn0[hit[hit > 0][1]]
  }

  p_raw_col <- find_col(c("p_pool", "pvalue", "p_value", "pval", "p.raw", "p_raw"))
  p_adj_col <- find_col(c("p_bh", "fdr", "padj", "p_adj", "adj.p.val", "qvalue", "q_value"))
  log2fc_col <- find_col(c("log2fc", "logfc", "log2foldchange"))
  if (is.null(p_raw_col) || is.null(p_adj_col)) {
    stop("TailEnrich output is missing p-value columns.")
  }

  p_raw <- suppressWarnings(as.numeric(te_res[[p_raw_col]]))
  p_adj <- suppressWarnings(as.numeric(te_res[[p_adj_col]]))
  p_raw[is.na(p_raw) | !is.finite(p_raw)] <- 1
  p_adj[is.na(p_adj) | !is.finite(p_adj)] <- 1

  log2fc <- if (is.null(log2fc_col)) rep(NA_real_, length(genes)) else suppressWarnings(as.numeric(te_res[[log2fc_col]]))
  score <- -log10(p_raw + 1e-15)

  te_extra <- te_res
  drop_cols <- unique(na.omit(c(p_raw_col, p_adj_col, log2fc_col, find_col(c("gene")))))
  keep_cols <- setdiff(colnames(te_extra), drop_cols)
  te_extra <- te_extra[, keep_cols, drop = FALSE]

  out_df <- data.frame(
    gene_id = genes,
    log2FC = log2fc,
    PValue = p_raw,
    FDR = p_adj,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (ncol(te_extra) > 0L) {
    out_df <- cbind(out_df, te_extra, stringsAsFactors = FALSE)
  }
  out_df$score <- score

  sig_df <- out_df[
    is.finite(out_df$FDR) &
      out_df$FDR < fdr_cutoff &
      is.finite(out_df$log2FC) &
      abs(out_df$log2FC) >= log2fc_cutoff,
    ,
    drop = FALSE
  ]
  list(full = out_df, sig = sig_df)
}

te_env <- new.env(parent = globalenv())
sys.source(tailenrich_script, envir = te_env)
te_entry <- resolve_ttailenrich_entry(te_env, tailenrich_script)
te_fun <- get(te_entry, envir = te_env)

if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

message("Using TailEnrich function: ", te_entry)
message("Input root: ", input_root)
message("Output root: ", output_root)
message("Datasets: ", paste(dataset_dirs, collapse = ", "))

for (dataset_name in dataset_dirs) {
  dataset_dir <- file.path(input_root, dataset_name)
  stop_if_not_exists(dataset_dir, "dataset directory")

  expr_file <- file.path(dataset_dir, expr_filename)
  group_file <- file.path(dataset_dir, group_filename)
  stop_if_not_exists(expr_file, expr_filename)
  stop_if_not_exists(group_file, group_filename)

  expr_mat <- read_expr_matrix(expr_file)
  group_df <- read_group_table(group_file)
  aligned <- align_expr_and_group(expr_mat, group_df, dataset_name)

  out_dir <- file.path(output_root, dataset_name, "tailEnrich_output")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  t0 <- proc.time()[3]
    res <- run_tailenrich_one(
      counts = aligned$expr,
      y = aligned$y,
      te_fun = te_fun,
      seed = seed,
      n_cores = n_cores,
      n_perm = n_perm,
      fdr_cutoff = fdr_cutoff,
      log2fc_cutoff = log2fc_cutoff
    )
  t1 <- proc.time()[3]

  write.table(res$full, file = file.path(out_dir, "tailEnrich.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(res$sig, file = file.path(out_dir, "tailEnrich_sig.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  message(sprintf(
    "[DONE] %s | genes=%d | samples=%d | sig=%d | %.3f sec",
    dataset_name,
    nrow(aligned$expr),
    ncol(aligned$expr),
    nrow(res$sig),
    as.numeric(t1 - t0)
  ))
}
