# Output Columns

The R return value and exported TSV files share the following seven columns. The full table retains all tested genes; the significant table contains only the rows where `is_significant` is `TRUE`.

| Column | Meaning |
| --- | --- |
| `gene_id` | Input gene identifier |
| `direction` | `left` indicates enrichment at the low-expression tail; `right` indicates enrichment at the high-expression tail |
| `enrichment_score` | Larger of the left- and right-tail statistics |
| `log2fc` | Log2 ratio of the disease-sample mean to the selected-control mean, with 1 added to each mean |
| `p_value` | Pooled permutation p value for the larger of the two tail statistics |
| `fdr` | Benjamini-Hochberg-adjusted p value across all tested genes |
| `is_significant` | Whether both `fdr < fdr_cutoff` and `abs(log2fc) >= log2fc_cutoff` hold |

Each tail statistic is the maximum normalized running-score height multiplied by one minus the normalized position where that maximum is first reached. Scan positions range from 0 to 1. The larger statistic determines the direction; ties are assigned to the left. For each label permutation, the larger tail statistic is also calculated for every gene, then pooled across genes and permutations to form the null distribution. The function reports one gene-level p value, not separate directional p values.

For `log2fc`, the disease-sample mean is calculated over the region from the selected tail to the first running-score maximum, including that position. If this region contains n disease samples, the control-sample mean uses the first n controls encountered from the same tail, or all controls if fewer are available. These controls need not lie outside the enriched region. A pseudocount of 1 is added to both means before taking their log2 ratio. This effect size differs from the overall mean difference between disease and control samples, and its sign does not define enrichment direction.

Peak coordinates, areas, control-selection positions, and other intermediate diagnostics are not included in the main result table. Data used to plot the example curves are stored separately under `sample_output/Mayo_TCX/display/`.
