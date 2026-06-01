# rsomics-tmm-norm

TMM (trimmed mean of M-values) per-sample normalization factors for a gene x
sample count matrix — a Rust port of edgeR's `calcNormFactors(method = "TMM")`.

```
rsomics-tmm-norm counts.tsv -o factors.tsv
```

Input is a TSV: a header row of sample IDs (the first column header is the gene
column, conventionally `gene`), then one row per gene with the gene ID followed
by integer-ish counts. Output is `sample<TAB>norm.factor`, one row per sample,
scaled so the factors' geometric mean is 1.

## Method

All-zero genes are dropped. Library sizes are column sums over the remaining
genes. The reference sample is the one whose 0.75 quantile of counts divided by
its library size (`f75`) is closest to the mean `f75`; when the median `f75` is
effectively zero the reference falls back to the column of largest summed
square-root mass. For each sample the factor is the variance-weighted trimmed
mean of the per-gene log2 ratios M = log2((obs/Nobs)/(ref/Nref)) against the
reference, trimming `logratioTrim = 0.3` of M and `sumTrim = 0.05` of the mean
abundance A, using asymptotic-variance weights `1/v`. Genes with a zero count in
either the sample or the reference fall out via the finite-value filter on M.

## Origin

This crate is an independent Rust reimplementation of edgeR's TMM normalization
based on:
- The published method: Robinson MD, Oshlack A. "A scaling normalization method
  for differential expression analysis of RNA-seq data." Genome Biology 2010,
  11:R25. doi:10.1186/gb-2010-11-3-r25 — which gives the M/A log-ratio
  definitions, the asymptotic-variance weights, the double trim, and the
  reference-by-upper-quartile selection.
- The public conventions of R's `rank` (average ties) and `quantile` type-7,
  and the documented default parameters (`logratioTrim = 0.3`,
  `sumTrim = 0.05`).
- Black-box behaviour testing against the upstream binary.

No source code from the GPL upstream (edgeR) was used as reference during
implementation. Test fixtures are independently generated (a fixed-seed
negative-binomial count matrix), and the reference values are captured by
running the upstream binary as a black box (`tests/tmm_oracle.R` via
`calcNormFactors`), with a committed golden so CI validates without R.

License: MIT OR Apache-2.0.
Upstream credit: edgeR (https://bioconductor.org/packages/edgeR/), GPL (>=2).
