############################################################
# TailEnrich standalone
#
# Independent script:
#   - does NOT require source("tailEnrich.R")
#   - directly provides TailEnrich()
#
# Statistic:
#   - L2H statistic := h * (1 - x)
#   - H2L statistic := h * (1 - x)
#   - best-stat     := max(PE_L2H, PE_H2L)
#
# Logic:
#   - U-shape removed
#   - log2FC computed only for best direction
#   - p-values are based only on pooled best-stat permutation
#   - report a permutation p value and BH-adjusted FDR
#
# Tie-break:
#   - uses a group-independent reproducible tie_key
#   - avoids tie clustering caused by original sample order
############################################################

suppressPackageStartupMessages({
  library(parallel)
})

# ============================================================
# Tie-break configuration
# ============================================================
.TE_USE_TIEBREAK  <- TRUE
.TE_TIEBREAK_SEED <- 20260205L

.make_tie_key <- function(n, seed = .TE_TIEBREAK_SEED) {
  if (!isTRUE(.TE_USE_TIEBREAK)) return(NULL)

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }

  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(as.integer(seed[1]))
  sample.int(n)
}

# ============================================================
# Core PE score
# ============================================================
.norm_pos_01 <- function(max_pos, n) {
  if (n <= 1L) return(0)
  (max_pos - 1) / (n - 1)
}

compute_tail_score <- function(gene_expr, y, tie_key = NULL) {
  nd <- sum(y == 1)
  n  <- length(y)

  sorted_idx <- if (is.null(tie_key)) {
    order(gene_expr)
  } else {
    order(gene_expr, tie_key)
  }

  sorted_gene <- gene_expr[sorted_idx]
  y_sorted    <- y[sorted_idx]

  y_center <- y_sorted - mean(y_sorted)
  cum_vals <- cumsum(y_center)
  denom    <- sum(y_center[y_center > 0])

  if (denom == 0) {
    ratios <- rep(0, n)
  } else {
    ratios <- cum_vals / denom
  }

  max_pos  <- which.max(ratios)
  peak_h   <- ratios[max_pos]
  peak_x   <- .norm_pos_01(max_pos, n)
  tail_score <- peak_h * (1 - peak_x)
  critical <- sorted_gene[max_pos]

  if (denom == 0 || nd < 1L) {
    area <- 0
  } else {
    k <- min(nd, n)
    area <- sum(ratios[seq_len(k)]) / n
  }

  c(tail_score, critical, area, peak_h, peak_x)
}

tailenrich_core <- function(X, y, tie_key = NULL) {
  X <- as.matrix(X)
  y <- as.numeric(y)

  p <- nrow(X)
  n <- ncol(X)

  if (length(y) != n) stop("length(y) must equal ncol(X).")

  PE <- matrix(0, p, 10)
  rownames(PE) <- rownames(X)

  for (j in seq_len(p)) {
    gene_expr <- X[j, ]

    rt_low  <- compute_tail_score(gene_expr,  y, tie_key = tie_key)
    rt_high <- compute_tail_score(-gene_expr, y, tie_key = if (is.null(tie_key)) NULL else -tie_key)

    PE[j, ] <- c(
      rt_low[1], rt_low[2], rt_low[3], rt_low[4], rt_low[5],
      rt_high[1], -rt_high[2], rt_high[3], rt_high[4], rt_high[5]
    )
  }

  PE <- as.data.frame(PE, stringsAsFactors = FALSE)
  colnames(PE) <- c(
    "PE_L2H", "Critical_L2H", "Area_L2H", "PeakH_L2H", "PeakX_L2H",
    "PE_H2L", "Critical_H2L", "Area_H2L", "PeakH_H2L", "PeakX_H2L"
  )
  PE
}

# ============================================================
# Precompute orders
# ============================================================
precompute_order_L2H <- function(X, tie_key = NULL) {
  X <- as.matrix(X)
  p <- nrow(X)
  n <- ncol(X)

  ordL <- matrix(0L, p, n)
  for (j in seq_len(p)) {
    ordL[j, ] <- if (is.null(tie_key)) {
      order(X[j, ])
    } else {
      order(X[j, ], tie_key)
    }
  }
  ordL
}

precompute_order_H2L <- function(X, tie_key = NULL) {
  X <- as.matrix(X)
  p <- nrow(X)
  n <- ncol(X)

  ordH <- matrix(0L, p, n)
  for (j in seq_len(p)) {
    ordH[j, ] <- if (is.null(tie_key)) {
      order(-X[j, ])
    } else {
      order(-X[j, ], -tie_key)
    }
  }
  ordH
}

# ============================================================
# hmax-defined log2FC for the best direction
# Denominator controls are selected from the enrichment-tail start,
# without requiring them to lie outside 1:k. The target count equals
# n_exp_inside. If fewer controls are available, use all available
# controls; if none are available, set log2FC to zero.
# ============================================================
.summarize_hmax_region <- function(x, y_sorted, ord, k) {
  inside_positions <- seq_len(k)
  inside_exp_positions <- inside_positions[y_sorted[inside_positions] == 1]
  tail_ctrl_positions <- which(y_sorted == -1)

  n_exp_inside <- length(inside_exp_positions)
  n_ctrl_target <- n_exp_inside
  n_ctrl_available <- length(tail_ctrl_positions)
  n_ctrl_selected <- min(n_ctrl_target, n_ctrl_available)
  selected_ctrl_positions <- if (n_ctrl_selected > 0L) {
    tail_ctrl_positions[seq_len(n_ctrl_selected)]
  } else {
    integer(0)
  }
  selected_ctrl_first_position <- if (n_ctrl_selected > 0L) selected_ctrl_positions[1] else NA_integer_
  selected_ctrl_last_position <- if (n_ctrl_selected > 0L) selected_ctrl_positions[n_ctrl_selected] else NA_integer_

  list(
    n_exp_inside = n_exp_inside,
    n_ctrl_target = n_ctrl_target,
    n_ctrl_available = n_ctrl_available,
    n_ctrl_selected = n_ctrl_selected,
    ctrl_selection_complete = n_ctrl_selected == n_ctrl_target,
    selected_ctrl_first_position = selected_ctrl_first_position,
    selected_ctrl_last_position = selected_ctrl_last_position,
    selected_ctrl_positions = paste(selected_ctrl_positions, collapse = ","),
    enriched_mean_exp = if (n_exp_inside > 0L) mean(x[ord[inside_exp_positions]]) else NA_real_,
    denominator_mean_ctrl = if (n_ctrl_selected > 0L) mean(x[ord[selected_ctrl_positions]]) else NA_real_,
    zero_ctrl_log2fc_fallback = n_ctrl_selected == 0L
  )
}

.compute_hmax_log2fc_one <- function(x, y, ord, peak_x_expected,
                                    direction, peak_tol = 1e-12) {
  n <- length(ord)
  if (length(x) != n || length(y) != n) {
    stop("x, y, and ord must have the same length.")
  }
  if (!identical(sort(as.integer(ord)), seq_len(n))) {
    stop("ord must be a permutation of seq_len(length(x)).")
  }

  y_sorted <- y[ord]
  y_centered <- y_sorted - mean(y_sorted)
  denom_h <- sum(y_centered[y_centered > 0])
  if (!is.finite(denom_h) || denom_h <= 0) {
    stop("Cannot determine hmax: centered-label denominator is not positive.")
  }

  running_score <- cumsum(y_centered) / denom_h
  k <- which.max(running_score)
  peak_x_observed <- .norm_pos_01(k, n)
  peak_matches <- is.finite(peak_x_expected) &&
    abs(peak_x_observed - peak_x_expected) < peak_tol

  region <- .summarize_hmax_region(x = x, y_sorted = y_sorted, ord = ord, k = k)
  n_exp_inside <- region$n_exp_inside
  M <- region$enriched_mean_exp
  C <- region$denominator_mean_ctrl

  valid_k <- length(k) == 1L && is.finite(k) && k >= 1L && k <= n
  valid_numerator <- is.finite(M) && M >= 0
  valid_denominator <- region$n_ctrl_selected == 0L || (is.finite(C) && C >= 0)
  diagnostic_ok <- valid_k && n_exp_inside > 0L && valid_numerator &&
    valid_denominator && peak_matches
  log2fc <- if (!diagnostic_ok) {
    NA_real_
  } else if (region$zero_ctrl_log2fc_fallback) {
    0
  } else {
    log2((M + 1) / (C + 1))
  }

  list(
    log2FC = log2fc,
    peak_position = as.integer(k),
    peak_x_observed = peak_x_observed,
    peak_x_expected = peak_x_expected,
    peak_matches = peak_matches,
    n_exp_inside = n_exp_inside,
    n_ctrl_target = region$n_ctrl_target,
    n_ctrl_available = region$n_ctrl_available,
    n_ctrl_selected = region$n_ctrl_selected,
    ctrl_selection_complete = region$ctrl_selection_complete,
    selected_ctrl_first_position = region$selected_ctrl_first_position,
    selected_ctrl_last_position = region$selected_ctrl_last_position,
    selected_ctrl_positions = region$selected_ctrl_positions,
    enriched_mean_exp = M,
    denominator_mean_ctrl = C,
    zero_ctrl_log2fc_fallback = region$zero_ctrl_log2fc_fallback,
    diagnostic_ok = diagnostic_ok,
    direction = direction,
    sorted_group_summary = paste(y_sorted, collapse = ",")
  )
}

.compute_hmax_log2fc_all <- function(X, y, ordL, ordH, PE_obs,
                                     direction_best, peak_tol = 1e-12) {
  p <- nrow(X)
  out <- vector("list", p)

  for (j in seq_len(p)) {
    direction <- direction_best[j]
    if (!direction %in% c("L2H", "H2L")) {
      stop("Unknown direction_best at row ", j, ": ", direction)
    }
    ord <- if (direction == "L2H") ordL[j, ] else ordH[j, ]
    peak_x_expected <- if (direction == "L2H") PE_obs$PeakX_L2H[j] else PE_obs$PeakX_H2L[j]
    out[[j]] <- .compute_hmax_log2fc_one(
      x = X[j, ],
      y = y,
      ord = ord,
      peak_x_expected = peak_x_expected,
      direction = direction,
      peak_tol = peak_tol
    )
  }

  data.frame(
    log2FC = vapply(out, `[[`, numeric(1), "log2FC"),
    peak_position = vapply(out, `[[`, integer(1), "peak_position"),
    peak_x_observed = vapply(out, `[[`, numeric(1), "peak_x_observed"),
    peak_x_expected = vapply(out, `[[`, numeric(1), "peak_x_expected"),
    peak_matches = vapply(out, `[[`, logical(1), "peak_matches"),
    n_exp_inside = vapply(out, `[[`, integer(1), "n_exp_inside"),
    n_ctrl_target = vapply(out, `[[`, integer(1), "n_ctrl_target"),
    n_ctrl_available = vapply(out, `[[`, integer(1), "n_ctrl_available"),
    n_ctrl_selected = vapply(out, `[[`, integer(1), "n_ctrl_selected"),
    ctrl_selection_complete = vapply(out, `[[`, logical(1), "ctrl_selection_complete"),
    selected_ctrl_first_position = vapply(out, `[[`, integer(1), "selected_ctrl_first_position"),
    selected_ctrl_last_position = vapply(out, `[[`, integer(1), "selected_ctrl_last_position"),
    selected_ctrl_positions = vapply(out, `[[`, character(1), "selected_ctrl_positions"),
    enriched_mean_exp = vapply(out, `[[`, numeric(1), "enriched_mean_exp"),
    denominator_mean_ctrl = vapply(out, `[[`, numeric(1), "denominator_mean_ctrl"),
    zero_ctrl_log2fc_fallback = vapply(out, `[[`, logical(1), "zero_ctrl_log2fc_fallback"),
    diagnostic_ok = vapply(out, `[[`, logical(1), "diagnostic_ok"),
    sorted_group_summary = vapply(out, `[[`, character(1), "sorted_group_summary"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

# ============================================================
# Scheduler-aware default cores
# ============================================================
.get_alloc_cores <- function() {
  env <- Sys.getenv(
    c("SLURM_CPUS_PER_TASK", "NSLOTS", "PBS_NP", "LSB_DJOB_NUMPROC", "OMP_NUM_THREADS"),
    unset = NA
  )
  env <- suppressWarnings(as.integer(env))
  env <- env[!is.na(env) & env > 0]

  if (length(env) > 0) return(env[1])

  dc <- parallel::detectCores(logical = TRUE)
  if (is.na(dc) || dc < 1) return(1L)
  as.integer(dc)
}

# ============================================================
# Fast tail score under precomputed orders
# ============================================================
compute_tail_by_ord_fast <- function(y, ord, mean_y, denom) {
  y_sorted <- y[ord]
  cum_vals <- cumsum(y_sorted - mean_y)
  max_pos <- which.max(cum_vals)
  peak_h <- max(cum_vals) / denom
  peak_x <- .norm_pos_01(max_pos, length(ord))
  peak_h * (1 - peak_x)
}

tail_all_genes_preord_fast <- function(y, ord_mat, mean_y, denom) {
  p <- nrow(ord_mat)
  out <- numeric(p)

  for (j in seq_len(p)) {
    out[j] <- compute_tail_by_ord_fast(y, ord_mat[j, ], mean_y, denom)
  }

  out
}

# ============================================================
# Main function
# ============================================================
TailEnrich <- function(expression_matrix, sample_groups,
                       n_permutations = NULL, seed = NULL, n_cores = 1L,
                       fdr_cutoff = 0.05, log2fc_cutoff = log2(1.5)) {
  X <- as.matrix(expression_matrix)
  y <- as.numeric(sample_groups)

  if (!is.numeric(X) || !nrow(X) || !ncol(X) || any(!is.finite(X)) || any(X < 0)) {
    stop("expression_matrix must contain finite, nonnegative numeric values.")
  }
  if (is.null(rownames(X)) || anyNA(rownames(X)) || any(!nzchar(rownames(X))) || anyDuplicated(rownames(X))) {
    stop("expression_matrix must have unique, nonempty gene IDs as row names.")
  }

  if (!all(y %in% c(-1, 1))) {
    stop("sample_groups must use 1 for disease samples and -1 for control samples.")
  }

  p <- nrow(X)
  n <- ncol(X)

  if (length(y) != n) {
    stop("length(sample_groups) must equal ncol(expression_matrix).")
  }

  if (length(unique(y)) != 2L) stop("sample_groups must contain both disease and control samples.")
  if (!is.numeric(fdr_cutoff) || length(fdr_cutoff) != 1L || !is.finite(fdr_cutoff) ||
      fdr_cutoff <= 0 || fdr_cutoff > 1) {
    stop("fdr_cutoff must be a number in (0, 1].")
  }
  if (!is.numeric(log2fc_cutoff) || length(log2fc_cutoff) != 1L ||
      !is.finite(log2fc_cutoff) || log2fc_cutoff < 0) {
    stop("log2fc_cutoff must be a finite, nonnegative number.")
  }

  if (is.null(n_permutations)) {
    n_perm <- max(50L, as.integer(ceiling(2.5e6 / p)))
  } else {
    if (!is.numeric(n_permutations) || length(n_permutations) != 1L ||
        !is.finite(n_permutations) || n_permutations < 1 ||
        n_permutations != floor(n_permutations) || n_permutations > .Machine$integer.max) {
      stop("n_permutations must be a positive integer.")
    }
    n_perm <- as.integer(n_permutations)
  }

  if (is.null(n_cores)) {
    alloc <- .get_alloc_cores()
    n_cores <- max(1L, alloc - 1L)
  } else {
    if (!is.numeric(n_cores) || length(n_cores) != 1L || !is.finite(n_cores) ||
        n_cores < 1 || n_cores != floor(n_cores) || n_cores > .Machine$integer.max) {
      stop("n_cores must be a positive integer.")
    }
    n_cores <- as.integer(n_cores)
  }
  if (.Platform$OS.type == "windows" && n_cores > 1L) stop("Use n_cores = 1 on Windows.")
  n_jobs <- min(n_cores, n_perm)

  base <- n_perm %/% n_jobs
  rem  <- n_perm %% n_jobs
  chunk_sizes <- rep(base, n_jobs)
  if (rem > 0L) chunk_sizes[seq_len(rem)] <- chunk_sizes[seq_len(rem)] + 1L

  if (!is.null(seed)) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(as.integer(seed[1]))
  }

  tie_key <- .make_tie_key(n)

  ordL <- precompute_order_L2H(X, tie_key = tie_key)
  ordH <- precompute_order_H2L(X, tie_key = tie_key)

  PE_obs <- tailenrich_core(X, y, tie_key = tie_key)

  T_L_obs <- PE_obs$PE_L2H
  T_R_obs <- PE_obs$PE_H2L

  T_best_obs     <- pmax(T_L_obs, T_R_obs)
  direction_best <- ifelse(T_L_obs >= T_R_obs, "L2H", "H2L")

  # ----------------------------------------------------------
  # TailEnrich log2FC for the best direction, with hmax defining 1:k.
  # Numerator: mean expression of experimental samples in 1:k.
  # Denominator: mean of the first n_exp_inside controls encountered from
  # the same enrichment-tail start, without restricting controls to after k.
  # If fewer controls are available, use all available controls; if none
  # are available, log2FC is zero.
  # The hmax sample is included in the enriched region.
  # ----------------------------------------------------------
  fc_diag <- .compute_hmax_log2fc_all(
    X = X,
    y = y,
    ordL = ordL,
    ordH = ordH,
    PE_obs = PE_obs,
    direction_best = direction_best
  )

  if (any(!fc_diag$diagnostic_ok)) {
    bad <- which(!fc_diag$diagnostic_ok)
    bad_gene <- rownames(X)[bad]
    bad_table <- data.frame(
      gene_id = bad_gene,
      direction = direction_best[bad],
      peak_position = fc_diag$peak_position[bad],
      n_exp_inside = fc_diag$n_exp_inside[bad],
      n_ctrl_target = fc_diag$n_ctrl_target[bad],
      n_ctrl_available = fc_diag$n_ctrl_available[bad],
      n_ctrl_selected = fc_diag$n_ctrl_selected[bad],
      ctrl_selection_complete = fc_diag$ctrl_selection_complete[bad],
      peak_matches = fc_diag$peak_matches[bad],
      enriched_mean_exp = fc_diag$enriched_mean_exp[bad],
      denominator_mean_ctrl = fc_diag$denominator_mean_ctrl[bad],
      zero_ctrl_log2fc_fallback = fc_diag$zero_ctrl_log2fc_fallback[bad],
      sorted_group_summary = fc_diag$sorted_group_summary[bad],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    stop(
      "TailEnrich hmax log2FC diagnostics failed for ", length(bad), " gene(s).\n",
      paste(capture.output(print(utils::head(bad_table, 20), row.names = FALSE)), collapse = "\n")
    )
  }

  log2FC <- fc_diag$log2FC

  # ----------------------------------------------------------
  # constants for fast permutation
  # ----------------------------------------------------------
  nd <- sum(y == 1)
  if (nd == 0L || nd == n) stop("Invalid y: need both classes present.")

  mean_y <- mean(y)
  denom  <- nd * (1 - mean_y)
  if (denom == 0) stop("Invalid y: denom==0")

  perm_worker <- function(B) {
    exceed_best <- numeric(p)

    for (b in seq_len(B)) {
      y_perm <- sample(y)

      peL_perm <- tail_all_genes_preord_fast(y_perm, ordL, mean_y, denom)
      peR_perm <- tail_all_genes_preord_fast(y_perm, ordH, mean_y, denom)

      T_best_perm <- pmax(peL_perm, peR_perm)

      sp_best <- sort(T_best_perm)
      n_lt    <- findInterval(T_best_obs, sp_best, left.open = TRUE)
      exceed_best <- exceed_best + (p - n_lt)
    }

    list(best = exceed_best)
  }

  exceed_list <- parallel::mclapply(
    X           = chunk_sizes,
    FUN         = perm_worker,
    mc.cores    = n_jobs,
    mc.set.seed = TRUE
  )

  if (any(vapply(exceed_list, inherits, logical(1), what = "try-error"))) {
    stop("At least one mclapply worker failed. Run warnings(); then re-run with n_cores=1 to see the exact error.")
  }

  ex_best <- Reduce(`+`, lapply(exceed_list, `[[`, "best"))

  denom_pool <- n_perm * p + 1
  p_pool <- (ex_best + 1) / denom_pool
  p_bh   <- p.adjust(p_pool, method = "BH")

  res <- data.frame(
    gene_id = rownames(X),
    direction = ifelse(direction_best == "L2H", "left", "right"),
    enrichment_score = T_best_obs,
    log2fc = log2FC,
    p_value = p_pool,
    fdr = p_bh,
    is_significant = is.finite(p_bh) & p_bh < fdr_cutoff &
      is.finite(log2FC) & abs(log2FC) >= log2fc_cutoff,
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )

  rownames(res) <- rownames(X)
  as.data.frame(res)
}

