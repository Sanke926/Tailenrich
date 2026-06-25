# TailEnrich

This package is a GitHub-ready public bundle for the `TailEnrich` method. It includes:

- the core R implementation,
- a portable shell command for rerunning the method,
- one representative real-data sample input,
- the covariate-adjustment code and raw files used to prepare that sample input,
- the corresponding sample output,
- a concise README with method description, dependencies, figures, and result interpretation.

## Repository Contents

- `src/tailEnrich.R`: core TailEnrich implementation.
- `scripts/run_tailEnrich.R`: repository-relative R runner.
- `scripts/run_sample.sh`: one-command example for the bundled sample dataset.
- `scripts/prepare_covariate_adjusted_input.R`: sample-level covariate-adjustment script used before TailEnrich.
- `scripts/plot_sample_volcano.py`: volcano plot for the bundled sample dataset.
- `scripts/plot_selected_pe_curves.py`: annotated PE-curve panel for selected genes.
- `sample_input/MSBB-BM36/`: representative input dataset.
- `sample_output/MSBB-BM36/`: TailEnrich output plus one `edgeR` reference table.
- `results/selected_genes.tsv`: genes used for the PE-curve figure.
- `figures/`: README figures.

## Core Dependencies

### Required to run TailEnrich

- `bash` or `sh`
- `R >= 4.2`
- base R `parallel`

The core TailEnrich run in this package does not require extra CRAN or Bioconductor packages.

### Optional for rebuilding figures

- `Python >= 3.10`
- `matplotlib`

## Input Format

Each dataset directory should contain:

1. `gene_TPM_by_salmon_covAdjusted.csv`
   - first column: gene ID
   - remaining columns: expression values
   - header row: sample IDs
2. `used_samples_group.tsv`
   - required columns: `sampleID`, `group`
   - `group` must be coded as `-1` and `1`

For the bundled public sample, the repository also includes the files used to generate the covariate-adjusted matrix:

- `gene_TPM_by_salmon.csv`: raw TPM matrix before covariate removal
- `MSBB36_control_isoform_covarites.csv`: sample-level covariate table
- `isoform_sample.txt`: original case/control sample definition
- `used_covariates.txt`: covariates selected for removal in this dataset

The bundled sample dataset is:

- `sample_input/MSBB-BM36/gene_TPM_by_salmon.csv`
- `sample_input/MSBB-BM36/gene_TPM_by_salmon_covAdjusted.csv`
- `sample_input/MSBB-BM36/MSBB36_control_isoform_covarites.csv`
- `sample_input/MSBB-BM36/isoform_sample.txt`
- `sample_input/MSBB-BM36/used_covariates.txt`
- `sample_input/MSBB-BM36/used_samples_group.tsv`

## Quick Start

Run the bundled sample:

```bash
bash scripts/run_sample.sh
```

If `Rscript` is not on `PATH`:

```bash
RSCRIPT_BIN=/path/to/Rscript bash scripts/run_sample.sh
```

Main environment variables:

- `TAILENRICH_INPUT_ROOT`: dataset root, default `sample_input`
- `TAILENRICH_OUTPUT_ROOT`: output root, default `sample_output_rerun`
- `TAILENRICH_DATASETS`: comma-separated dataset names, default `MSBB-BM36`
- `TAILENRICH_N_CORES`: number of cores, default `1`
- `TAILENRICH_N_PERM`: permutation count, default `1000`
- `TAILENRICH_SEED`: random seed, default `1`
- `TAILENRICH_FDR_CUTOFF`: significance cutoff, default `0.05`

## Preprocessing Workflow

The bundled sample follows the same analysis sequence used in the full real-data workflow:

1. Start from the raw TPM matrix `gene_TPM_by_salmon.csv`.
2. Read the sample-level covariate table `MSBB36_control_isoform_covarites.csv`.
3. Keep the covariates listed in `used_covariates.txt`.
4. Use `isoform_sample.txt` to define `exp = 1` and `ctrl = -1`, then write `used_samples_group.tsv`.
5. Remove covariate effects from the TPM matrix and write `gene_TPM_by_salmon_covAdjusted.csv`.
6. Run TailEnrich on the covariate-adjusted TPM matrix.

The bundled preprocessing script reproduces this sample-level step:

```bash
Rscript scripts/prepare_covariate_adjusted_input.R \
  --expr sample_input/MSBB-BM36/gene_TPM_by_salmon.csv \
  --covariates sample_input/MSBB-BM36/MSBB36_control_isoform_covarites.csv \
  --sample-file sample_input/MSBB-BM36/isoform_sample.txt \
  --used-covariates sample_input/MSBB-BM36/used_covariates.txt \
  --out-expr sample_input/MSBB-BM36/gene_TPM_by_salmon_covAdjusted.csv \
  --out-groups sample_input/MSBB-BM36/used_samples_group.tsv
```

In the broader paper workflow, TailEnrich is run on the covariate-adjusted TPM matrix, while the comparison `edgeR` analysis is run on raw count data with its own covariate-aware design.

## Bundled Example

The bundled public example uses `MSBB-BM36`.

Bundled input:

- `sample_input/MSBB-BM36/gene_TPM_by_salmon.csv`
- `sample_input/MSBB-BM36/gene_TPM_by_salmon_covAdjusted.csv`
- `sample_input/MSBB-BM36/MSBB36_control_isoform_covarites.csv`
- `sample_input/MSBB-BM36/isoform_sample.txt`
- `sample_input/MSBB-BM36/used_covariates.txt`
- `sample_input/MSBB-BM36/used_samples_group.tsv`

Bundled output:

- `sample_output/MSBB-BM36/tailEnrich.tsv`
- `sample_output/MSBB-BM36/tailEnrich_sig.tsv`
- `sample_output/MSBB-BM36/edgeR.tsv`

In this sample dataset:

- TailEnrich significant genes: `253`
- `edgeR` significant genes: `1445`
- selected PE-curve examples: genes that are significant in TailEnrich but not significant in `edgeR`

## Output Columns

The main TailEnrich result table includes:

- `gene_id`
- `log2FC`
- `PValue`
- `FDR`
- `direction_best`
- `TE_score_any`
- `score`

Interpretation:

- for figure display, `right` denotes the high-expression tail.
- for figure display, `left` denotes the low-expression tail.
- `log2FC` is the direction-specific tail fold change.
- `PValue` and `FDR` summarize permutation-based significance.

## Figures

### Volcano Plot

Red points satisfy both `TailEnrich FDR < 0.05` and `|log2FC| >= log2(1.5)`.
Non-significant genes are shown in gray.

![MSBB-BM36 volcano plot](figures/MSBB-BM36_volcano.png)

### Annotated PE Curves

These genes were chosen from the same sample dataset under the rule:

- significant in TailEnrich,
- non-significant in `edgeR`,
- not chosen simply by top significance rank.

The exact list is stored in `results/selected_genes.tsv`.

![Selected PE curves](figures/MSBB-BM36_selected_pe_curves.png)

## Results and Interpretation

For the representative sample dataset `MSBB-BM36`, TailEnrich detects a substantial set of significant genes (`253`) while still leaving many genes outside the conventional `edgeR` hit list.

The point of the bundled PE-curve examples is not to show the strongest possible genes, but to show genes where TailEnrich identifies tail-structured signal while `edgeR` remains non-significant. That is the main behavior this public package is intended to illustrate.

In practical terms, the example suggests:

- TailEnrich can recover genes with subgroup-tail structure that do not necessarily appear as standard `edgeR` hits.
- The method should be interpreted as complementary to conventional DE analysis, not merely as another ranking of the same global-shift signal.

## Rebuild Figures

```bash
python scripts/plot_sample_volcano.py
python scripts/plot_selected_pe_curves.py
```

## Directory Layout

```text
TailEnrich_github_repo/
|- README.md
|- figures/
|- results/
|- sample_input/
|- sample_output/
|- scripts/
`- src/
```
