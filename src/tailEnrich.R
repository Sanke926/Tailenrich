############################################################
# TailEnrich standalone (clean version)
#
# Independent script:
#   - does NOT require source("tailEnrich.R")
#   - directly provides tailenrich()
#
# Statistic:
#   - L2H statistic := PE_L2H = h * (1 - x)
#   - H2L statistic := PE_H2L = h * (1 - x)
#   - best-stat     := max(PE_L2H, PE_H2L)
#
# Logic:
#   - U-shape removed
#   - log2FC computed only for best direction
#   - p-values are based only on pooled best-stat permutation
#   - keep only final p_pool and p_bh
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
  tail_score  <- peak_h * (1 - peak_x)
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
compute_tail_score_by_order_fast <- function(y, ord, mean_y, denom) {
  y_sorted <- y[ord]
  cum_vals <- cumsum(y_sorted - mean_y)
  max_pos <- which.max(cum_vals)
  peak_h <- max(cum_vals) / denom
  peak_x <- .norm_pos_01(max_pos, length(ord))
  peak_h * (1 - peak_x)
}

tail_scores_all_genes_by_order_fast <- function(y, ord_mat, mean_y, denom) {
  p <- nrow(ord_mat)
  out <- numeric(p)

  for (j in seq_len(p)) {
    out[j] <- compute_tail_score_by_order_fast(y, ord_mat[j, ], mean_y, denom)
  }

  out
}

# ============================================================
# Main function
# ============================================================
tailenrich <- function(X, y,
                          n_perm    = NULL,
                          seed      = NULL,
                          n_cores   = NULL,
                          verbose   = FALSE,
                          fc_thresh = 1.5) {
  X <- as.matrix(X)
  y <- as.numeric(y)

  if (!all(y %in% c(-1, 1))) {
    stop("TailEnrich requires y in {-1, +1}. Please recode labels.")
  }

  p <- nrow(X)
  n <- ncol(X)

  if (length(y) != n) {
    stop("length(y) must equal ncol(X).")
  }

  if (!is.numeric(fc_thresh) || length(fc_thresh) != 1L || is.na(fc_thresh) || fc_thresh <= 1) {
    stop("fc_thresh must be a single numeric value > 1 (e.g., 1.5).")
  }
  fc_log2_thresh <- log2(fc_thresh)

  if (is.null(n_perm)) {
    n_perm <- max(50L, as.integer(ceiling(2.5e6 / p)))
  } else {
    n_perm <- as.integer(n_perm)
  }
  if (n_perm < 1L) stop("n_perm must be a positive integer.")

  if (is.null(n_cores)) {
    alloc <- .get_alloc_cores()
    n_cores <- max(1L, alloc - 1L)
  } else {
    n_cores <- max(1L, as.integer(n_cores))
  }
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

  # kept for backward compatibility with previous output schema
  TE_score_L2H <- abs(PE_obs$Area_L2H)

  # ----------------------------------------------------------
  # log2FC only for best direction
  #   - L2H: use exp samples lower than ctrl-min
  #   - H2L: use exp samples higher than ctrl-max
  #   - offset fixed to 1
  #   - empty tail => 0
  # ----------------------------------------------------------
  eps <- .Machine$double.eps
  safe_log2_ratio_offset1 <- function(M, C) {
    A <- max(M + 1, eps)
    B <- max(C + 1, eps)
    log2(A / B)
  }

  idx_exp  <- which(y ==  1)
  idx_ctrl <- which(y == -1)

  log2FC <- numeric(p)

  for (j in seq_len(p)) {
    x <- X[j, ]
    ctrl <- x[idx_ctrl]
    expv <- x[idx_exp]

    C_L <- suppressWarnings(min(ctrl, na.rm = TRUE))
    C_R <- suppressWarnings(max(ctrl, na.rm = TRUE))

    if (!is.finite(C_L) || !is.finite(C_R)) {
      log2FC[j] <- 0
      next
    }

    if (direction_best[j] == "L2H") {
      S <- expv[is.finite(expv) & (expv < C_L)]
      if (length(S) == 0L) {
        log2FC[j] <- 0
      } else {
        M <- mean(S, na.rm = TRUE)
        log2FC[j] <- if (is.finite(M)) safe_log2_ratio_offset1(M, C_L) else 0
      }
    } else {
      S <- expv[is.finite(expv) & (expv > C_R)]
      if (length(S) == 0L) {
        log2FC[j] <- 0
      } else {
        M <- mean(S, na.rm = TRUE)
        log2FC[j] <- if (is.finite(M)) safe_log2_ratio_offset1(M, C_R) else 0
      }
    }
  }

  pass_fc <- abs(log2FC) >= fc_log2_thresh

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

      peL_perm <- tail_scores_all_genes_by_order_fast(y_perm, ordL, mean_y, denom)
      peR_perm <- tail_scores_all_genes_by_order_fast(y_perm, ordH, mean_y, denom)

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
    PE_obs,
    direction_best = direction_best,
    TE_score_L2H   = TE_score_L2H,
    TE_score_any   = T_best_obs,
    log2FC         = log2FC,
    fc_thresh      = fc_thresh,
    pass_fc        = pass_fc,
    p_pool         = p_pool,
    p_bh           = p_bh,
    stat_name      = "PE",
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )

  rownames(res) <- rownames(X)
  as.data.frame(res)
}

