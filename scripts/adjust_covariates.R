#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: adjust_covariates.R gene_tpm.tsv groups.tsv covariates.csv output.csv")
raw <- read.delim(args[1], check.names = FALSE)
X <- as.matrix(raw[, -1, drop = FALSE])
rownames(X) <- raw[[1]]
if (!is.numeric(X) || any(!is.finite(X)) || any(X < 0)) stop("Invalid TPM matrix")
groups <- read.delim(args[2], check.names = FALSE)
cov <- read.csv(args[3], check.names = FALSE)
selected <- c("sex", "AgeAtDeath", "ApoE", "RIN", "Source")
if (!all(c("sampleID", "group") %in% names(groups))) stop("Invalid groups file")
if (!all(c("sampleID", selected) %in% names(cov))) stop("Missing covariate columns")
if (anyDuplicated(groups$sampleID) || anyDuplicated(cov$sampleID)) stop("Duplicate sample IDs")
gi <- match(colnames(X), groups$sampleID)
ci <- match(colnames(X), cov$sampleID)
if (anyNA(gi) || anyNA(ci)) stop("Unmatched sample IDs")
y <- as.numeric(groups$group[gi])
if (!all(y %in% c(-1, 1)) || length(unique(y)) != 2L) stop("Need both groups, coded -1 and 1")
cov <- cov[ci, selected, drop = FALSE]
if (any(!complete.cases(cov))) stop("Missing covariate values")
for (nm in names(cov)) {
  if (is.character(cov[[nm]]) || is.logical(cov[[nm]])) cov[[nm]] <- factor(cov[[nm]])
}
covmat <- model.matrix(~ ., data = cov)
covmat <- covmat[, colnames(covmat) != "(Intercept)", drop = FALSE]
design <- cbind(`(Intercept)` = 1, group = y, covmat)
if (qr(design)$rank != ncol(design)) stop("Covariate design is rank deficient")
log_tpm <- t(log2(X + 1))
coefficients <- solve(crossprod(design), crossprod(design, log_tpm))
adjusted_log <- log_tpm - design[, 3:ncol(design), drop = FALSE] %*%
  coefficients[3:ncol(design), , drop = FALSE]
adjusted <- t(2^adjusted_log - 1)
adjusted[adjusted < 0] <- 0
output <- data.frame(geneID = rownames(X), adjusted, check.names = FALSE)
dir.create(dirname(args[4]), recursive = TRUE, showWarnings = FALSE)
write.csv(output, args[4], row.names = FALSE, quote = FALSE)
