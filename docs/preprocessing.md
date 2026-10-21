# Mayo TCX Preprocessing

The preprocessing workflow starts from quantified transcript TPM and Salmon NumReads matrices. Transcripts are retained if they have TPM >= 1 in at least half the selected samples (75 of 149). Transcript version suffixes are removed before mapping, and values for retained transcripts are summed by gene. This gives 37,576 transcripts and 14,412 genes. Three retained transcripts missing from the NumReads matrix were recovered from the original per-sample quantification files. NumReads are used for the edgeR comparison; TailEnrich uses the TPM matrix.

The covariate-adjustment script takes the filtered gene TPM TSV file, the sample group TSV file, and a covariate CSV file containing `sampleID`, `sex`, `AgeAtDeath`, `ApoE`, `RIN`, and `Source`. Character and logical covariates are converted to factors; numeric covariates retain their input types. Samples are matched by identifier. Missing covariate values or a rank-deficient design cause an error.

The files needed to repeat this adjustment are included under `sample_input/Mayo_TCX/`:

| File | Contents |
| --- | --- |
| `gene_tpm_er_filtered.tsv` | Gene TPM values before covariate adjustment |
| `used_samples_group.tsv` | Sample identifiers and disease/control groups |
| `used_covariates.csv` | Sample identifiers and covariate values |
| `gene_TPM_covAdjusted.csv` | Adjusted matrix used by the example analysis |

Run the following command from the repository root. The regenerated matrix is written to `results/preprocessing/` to preserve the supplied example input.

```bash
Rscript scripts/adjust_covariates.R sample_input/Mayo_TCX/gene_tpm_er_filtered.tsv sample_input/Mayo_TCX/used_samples_group.tsv sample_input/Mayo_TCX/used_covariates.csv results/preprocessing/gene_TPM_covAdjusted.csv
```

The linear regression model includes an intercept, disease/control group, and covariates. Only the fitted covariate contribution is subtracted from log2(TPM + 1); the intercept and group contribution are retained. Values are transformed back to the expression scale, and negative values are set to zero. This adjustment is performed before calling `TailEnrich()`.

```text
Z = \log_2(X + 1)
Z_{adjusted} = Z - C\widehat{\beta}_C
X_{adjusted} = \max(2^{Z_{adjusted}} - 1, 0)
```

This workflow starts from quantified expression matrices and does not include read alignment or transcript quantification.
