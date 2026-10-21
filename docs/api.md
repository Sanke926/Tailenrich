# Function Arguments

```r
TailEnrich(expression_matrix, sample_groups,
           n_permutations = NULL, seed = NULL, n_cores = 1L,
           fdr_cutoff = 0.05, log2fc_cutoff = log2(1.5))
```

| Argument | Meaning | Default |
| --- | --- | --- |
| `expression_matrix` | Numeric matrix with genes in rows and samples in columns; row names must be unique gene IDs | Required |
| `sample_groups` | Labels in matrix column order, where `1` represents disease samples and `-1` represents control samples | Required |
| `n_permutations` | Number of label permutations | `max(50, ceiling(2.5e6 / nrow(expression_matrix)))` |
| `seed` | Random seed for label permutations | `NULL`; the example uses `1` |
| `n_cores` | Number of workers | `1` |
| `fdr_cutoff` | Upper FDR cutoff; a gene must have `fdr < fdr_cutoff` | `0.05` |
| `log2fc_cutoff` | Lower absolute log2FC cutoff; a gene must have `abs(log2fc) >= log2fc_cutoff` | `log2(1.5)` |

The function returns every tested gene, using the seven columns described in [Output Columns](output_fields.md). The two cutoffs determine `is_significant`; they do not filter the returned rows or change the statistics or p values.

The example was run on Windows with R 4.5.2, 174 permutations, seed 1, and one worker. Keep the input order, seed, R environment, and worker count fixed to reproduce the results exactly. Multiple workers use separate permutation streams and are supported on Unix-like systems; use one worker on Windows. Expression ties are resolved using a fixed ordering key generated independently of group labels with seed 20260205.
