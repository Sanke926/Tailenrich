#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Usage: export_display_data.R expression.csv groups.tsv results.tsv edger_reference.tsv output_dir")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "..", "src", "TailEnrich.R"))
expression_matrix <- as.matrix(read.csv(args[1], row.names = 1, check.names = FALSE))
group_table <- read.delim(args[2], check.names = FALSE)
idx <- match(colnames(expression_matrix), group_table$sampleID)
if (anyNA(idx)) stop("Unmatched samples")
sample_groups <- as.numeric(group_table$group[idx])
results <- read.delim(args[3], check.names = FALSE)
reference <- read.delim(args[4], check.names = FALSE)
idx <- match(reference$gene_id, results$gene_id)
if (anyNA(idx) || anyDuplicated(reference$gene_id)) stop("Invalid display gene reference")
selected <- results[idx, , drop = FALSE]
selected$edger_fdr <- reference$edger_fdr
symbols <- read.delim(file.path(args[5], "gene_symbols.tsv"), check.names = FALSE)
if (anyDuplicated(symbols$gene_id)) stop("Duplicate gene IDs in gene_symbols.tsv")
selected$gene_symbol <- symbols$gene_symbol[match(selected$gene_id, symbols$gene_id)]
if (anyNA(selected$gene_symbol) || any(!nzchar(selected$gene_symbol))) stop("Missing gene symbols")
if (!all(selected$is_significant & is.finite(selected$edger_fdr) & selected$edger_fdr >= 0.05)) {
  stop("Display genes must be TailEnrich significant and edgeR nonsignificant")
}
selected$enrichment <- paste0(selected$direction, "-tail enriched gene")
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
tie_key <- .make_tie_key(ncol(expression_matrix))
parts <- list()
for (i in seq_len(nrow(selected))) {
  gene <- selected$gene_id[i]
  direction <- selected$direction[i]
  values <- as.numeric(expression_matrix[gene, ])
  observed <- tailenrich_core(matrix(values, nrow = 1), sample_groups, tie_key)
  best_direction <- if (observed$PE_L2H >= observed$PE_H2L) "left" else "right"
  score <- max(observed$PE_L2H, observed$PE_H2L)
  if (direction != best_direction || !isTRUE(all.equal(score, selected$enrichment_score[i], tolerance = 1e-10))) {
    stop("Display curve differs from result: ", gene)
  }
  ord <- if (direction == "left") order(values, tie_key) else order(-values, -tie_key)
  groups <- sample_groups[ord]
  centered <- groups - mean(groups)
  running <- cumsum(centered) / sum(centered[centered > 0])
  peak <- which.max(running)
  disease_positions <- which(seq_along(ord) <= peak & groups == 1)
  control_positions <- head(which(groups == -1), length(disease_positions))
  reconstructed_fc <- log2((mean(values[ord[disease_positions]]) + 1) /
                           (mean(values[ord[control_positions]]) + 1))
  if (!isTRUE(all.equal(reconstructed_fc, selected$log2fc[i], tolerance = 1e-10))) {
    stop("Display log2fc differs from result: ", gene)
  }
  parts[[i]] <- data.frame(gene_id = gene, direction = direction,
    position = seq_along(ord), scan_x = (seq_along(ord) - 1) / (length(ord) - 1),
    running_score = running, group = groups)
}
dir.create(args[5], showWarnings = FALSE, recursive = TRUE)
write.table(selected, file.path(args[5], "display_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, parts), file.path(args[5], "display_curves.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("Exported selected-direction curves and verified enrichment scores and log2fc")
