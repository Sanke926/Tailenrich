#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    if (i == length(args)) {
      stop("Missing value for argument: ", key)
    }
    out[[substring(key, 3)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

stop_if_missing <- function(args, key) {
  if (is.null(args[[key]]) || !nzchar(args[[key]])) {
    stop("Missing required argument --", key)
  }
}

parse_sample_txt <- function(sample_file) {
  lines <- readLines(sample_file, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0L) {
    stop("Sample file is empty: ", sample_file)
  }

  group_code_map <- c("exp" = 1, "ctrl" = -1)
  sample_list <- vector("list", length(lines))

  for (i in seq_along(lines)) {
    tokens <- trimws(strsplit(lines[[i]], ",", fixed = TRUE)[[1]])
    tokens <- tokens[nzchar(tokens)]
    group_label <- tolower(tokens[[1]])
    if (!group_label %in% names(group_code_map)) {
      stop("Unsupported group label in sample file: ", tokens[[1]])
    }

    sample_list[[i]] <- data.frame(
      sampleID = tokens[-1],
      group = rep(group_code_map[[group_label]], length(tokens) - 1L),
      group_label = rep(group_label, length(tokens) - 1L),
      stringsAsFactors = FALSE
    )
  }

  sample_df <- do.call(rbind, sample_list)
  if (anyDuplicated(sample_df$sampleID)) {
    stop("Duplicate sample IDs found in sample file: ", sample_file)
  }
  sample_df
}

read_expression_matrix <- function(expr_file) {
  expr_df <- read.csv(expr_file, check.names = FALSE)
  expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
  suppressWarnings(storage.mode(expr_mat) <- "numeric")
  if (anyNA(expr_mat)) {
    stop("Expression matrix contains non-numeric values: ", expr_file)
  }
  rownames(expr_mat) <- expr_df[[1]]
  colnames(expr_mat) <- colnames(expr_df)[-1]
  list(feature_col = colnames(expr_df)[1], mat = expr_mat)
}

prepare_selected_covariates <- function(cov_df, selected_covariates) {
  if (length(selected_covariates) == 0L) {
    empty <- matrix(nrow = nrow(cov_df), ncol = 0L)
    attr(empty, "complete_idx") <- rep(TRUE, nrow(cov_df))
    return(empty)
  }

  selected_df <- cov_df[, selected_covariates, drop = FALSE]
  complete_idx <- complete.cases(selected_df)
  selected_df2 <- selected_df[complete_idx, , drop = FALSE]

  for (nm in colnames(selected_df2)) {
    if (is.character(selected_df2[[nm]]) || is.logical(selected_df2[[nm]])) {
      selected_df2[[nm]] <- as.factor(selected_df2[[nm]])
    }
  }

  mm <- model.matrix(~ ., data = selected_df2)
  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  }
  attr(mm, "complete_idx") <- complete_idx
  mm
}

adjust_counts <- function(counts, covmodel) {
  samples <- colnames(counts)
  log2counts <- t(log2(counts + 1))

  covmodel <- as.matrix(cbind(1, covmodel))
  B <- solve(t(covmodel) %*% covmodel) %*% t(covmodel) %*% log2counts

  if (ncol(covmodel) >= 3L) {
    yadjust <- log2counts - covmodel[, 3:ncol(covmodel), drop = FALSE] %*% B[3:ncol(covmodel), , drop = FALSE]
  } else {
    yadjust <- log2counts
  }

  counts_corrected <- t(2^(yadjust) - 1)
  counts_corrected[counts_corrected < 0] <- 0
  rownames(counts_corrected) <- rownames(counts)
  colnames(counts_corrected) <- samples
  counts_corrected
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  for (key in c("expr", "covariates", "sample-file", "used-covariates", "out-expr", "out-groups")) {
    stop_if_missing(args, key)
  }

  expr_info <- read_expression_matrix(args[["expr"]])
  expr_mat <- expr_info$mat
  sample_df <- parse_sample_txt(args[["sample-file"]])
  cov_df <- read.csv(args[["covariates"]], check.names = FALSE)
  selected_covariates <- readLines(args[["used-covariates"]], warn = FALSE)
  selected_covariates <- trimws(selected_covariates)
  selected_covariates <- selected_covariates[nzchar(selected_covariates)]

  idx_sample <- match(colnames(expr_mat), sample_df$sampleID)
  idx_cov <- match(colnames(expr_mat), cov_df$sampleID)
  if (anyNA(idx_sample)) {
    stop("Some expression samples are missing from the sample file.")
  }
  if (anyNA(idx_cov)) {
    stop("Some expression samples are missing from the covariate file.")
  }

  sample_df <- sample_df[idx_sample, , drop = FALSE]
  cov_df <- cov_df[idx_cov, , drop = FALSE]

  covariate_matrix <- prepare_selected_covariates(cov_df, selected_covariates)
  complete_idx <- attr(covariate_matrix, "complete_idx")

  expr_mat <- expr_mat[, complete_idx, drop = FALSE]
  sample_df <- sample_df[complete_idx, , drop = FALSE]
  group_vec <- as.numeric(sample_df$group)
  cov_model <- cbind(group = group_vec, covariate_matrix)
  adjusted <- adjust_counts(expr_mat, cov_model)

  out_expr <- data.frame(
    setNames(list(rownames(adjusted)), expr_info$feature_col),
    as.data.frame(adjusted, check.names = FALSE),
    check.names = FALSE
  )
  write.csv(out_expr, args[["out-expr"]], row.names = FALSE, quote = TRUE)
  write.table(sample_df, args[["out-groups"]], sep = "\t", row.names = FALSE, quote = FALSE)
}

main()
