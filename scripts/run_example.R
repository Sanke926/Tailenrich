#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 4L) {
  stop("Usage: Rscript scripts/run_example.R expression.csv groups.tsv [output_dir] [cores]")
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "..", "src", "TailEnrich.R"))
expr <- read.csv(args[1], row.names = 1, check.names = FALSE)
if (!all(vapply(expr, is.numeric, logical(1)))) stop("Expression columns must be numeric")
expression_matrix <- as.matrix(expr)
group <- read.delim(args[2], check.names = FALSE)
if (!all(c("sampleID", "group") %in% names(group))) stop("Need sampleID and group columns")
if (anyDuplicated(group$sampleID) || anyDuplicated(colnames(expression_matrix))) stop("Duplicate sample IDs")
if (anyDuplicated(rownames(expression_matrix)) || !nrow(expression_matrix)) stop("Need nonempty unique gene IDs")
idx <- match(colnames(expression_matrix), group$sampleID)
if (anyNA(idx)) stop("Expression samples missing from groups.tsv")
sample_groups <- as.numeric(group$group[idx])
if (any(!is.finite(expression_matrix)) || any(expression_matrix < 0)) stop("Expression must be finite and nonnegative")
cores <- if (length(args) >= 4L) as.integer(args[4]) else 1L
if (is.na(cores) || cores < 1L) stop("cores must be a positive integer")
if (.Platform$OS.type == "windows" && cores > 1L) stop("Windows requires cores=1")
outdir <- if (length(args) >= 3L) args[3] else "results/example"
tailenrich_output <- TailEnrich(expression_matrix, sample_groups, seed = 1, n_cores = cores)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(value, name) write.table(value, file.path(outdir, name),
  sep = "\t", quote = FALSE, row.names = FALSE)
write_tsv(tailenrich_output, "tailEnrich_results.tsv")
write_tsv(tailenrich_output[tailenrich_output$is_significant, , drop = FALSE], "tailEnrich_sig.tsv")
capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))
message(sprintf("Tested: %d; significant: %d",
                nrow(tailenrich_output), sum(tailenrich_output$is_significant)))
