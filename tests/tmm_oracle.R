#!/usr/bin/env Rscript
# edgeR calcNormFactors(method="TMM") reference. Emits sample<TAB>norm.factor,
# factors formatted to 10 decimals, matching rsomics-tmm-norm's output layout.
# Usage: Rscript tmm_oracle.R counts.tsv
suppressMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
counts_path <- args[1]

x <- read.delim(counts_path, row.names = 1, check.names = FALSE)
x <- as.matrix(x)

f <- calcNormFactors(x, method = "TMM")

cat("sample\tnorm.factor\n")
for (i in seq_along(f)) {
  cat(sprintf("%s\t%.10f\n", colnames(x)[i], f[i]))
}
