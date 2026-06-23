# PairMarker

`PairMarker` is an R package for pairwise gene marker discovery using binary
phenotype labels. It screens gene pairs using within-sample relative expression
comparisons, chi-square testing, Fisher's exact test, and multiple-testing
correction.

The package was originally developed for immunotherapy response analysis, but
the core workflow is applicable to any expression matrix with matched binary or
case-control annotations, such as responder versus non-responder, sensitive
versus resistant, disease versus normal, mutant versus wildtype, or high-risk
versus low-risk groups.

![PairMarker workflow overview](man/figures/immunePair-workflow.png)

## Installation

You can install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("twl-00/PairMarker", upgrade = "never")
```

## Quick Start

```r
library(PairMarker)

expr_file <- system.file(
  "extdata",
  "example_exp.txt",
  package = "PairMarker"
)

clinical_file <- system.file(
  "extdata",
  "example_clinical.txt",
  package = "PairMarker"
)

result <- run_pair_marker_analysis(
  expr_file = expr_file,
  clinical_file = clinical_file,
  out_dir = tempdir(),
  phenotype_col = "response",
  positive_label = "response",
  negative_label = "non_response",
  main_delta = 0.25,
  delta_list = c(0, 0.25, 0.5),
  dataset_name = "example"
)

result$sensitivity
head(result$pairs)
head(result$sig_pairs)
```

Basic downstream plots can be generated directly from the result object:

```r
plot_pair_roc(result, gene1 = "GENE_A", gene2 = "GENE_B")
```

If the clinical table contains survival time and event columns, a Kaplan-Meier
curve can also be drawn for a selected gene pair:

```r
plot_pair_survival(
  result,
  gene1 = "GENE_A",
  gene2 = "GENE_B",
  time_col = "OS_time",
  event_col = "OS_status"
)
```

## Response and Non-response Definitions

In `PairMarker`, `response` and `non_response` are example labels, not fixed
package-specific categories. The package only requires a binary phenotype
definition:

- `phenotype_col` is the column in the clinical table that contains the group
  label.
- `positive_label` is the label coded as 1. In an immunotherapy study, this is
  often the response group.
- `negative_label` is the label coded as 0. In an immunotherapy study, this is
  often the non-response group.

For immunotherapy datasets, the response group is usually defined by the
clinical response criteria used in the original study. For example, patients
with complete response or partial response may be coded as `response`, while
patients with progressive disease may be coded as `non_response`. Stable
disease can be handled according to the study design: it may be grouped with
response, grouped with non-response, or excluded by setting `negative_label`
explicitly and leaving stable-disease samples outside the two selected labels.

The label names do not need to be literally `response` and `non_response`.
For example, the same workflow can compare `sensitive` versus `resistant`,
`tumor` versus `normal`, `mutant` versus `wildtype`, or `high_risk` versus
`low_risk`. The important point is that the two groups should represent a
biologically or clinically meaningful contrast.

If `negative_label = NULL`, all samples that are not `positive_label` are coded
as 0. This is useful for one-versus-rest comparisons. If `negative_label` is
supplied, only samples matching `positive_label` or `negative_label` are used,
and samples with other labels are excluded from pair screening.

## Generic Binary Phenotype Examples

Common use cases:

| Application | `phenotype_col` | `positive_label` | `negative_label` | Interpretation |
| --- | --- | --- | --- | --- |
| Immunotherapy response | `"response"` | `"response"` | `"non_response"` | Gene pairs associated with treatment response status. |
| Drug sensitivity | `"drug_status"` | `"sensitive"` | `"resistant"` | Gene pairs associated with sensitivity or resistance to a drug. |
| Disease comparison | `"disease_status"` | `"tumor"` | `"normal"` | Gene pairs associated with disease versus control status. |
| Mutation status | `"TP53_status"` | `"mutant"` | `"wildtype"` | Gene pairs associated with a mutation-defined subgroup. |
| Risk group | `"risk_group"` | `"high"` | `"low"` | Gene pairs associated with high-risk versus low-risk labels. |
| Molecular subtype | `"subtype"` | `"classical"` | `NULL` | Gene pairs associated with one subtype versus all other subtypes. |

Immunotherapy response example:

```r
run_pair_marker_analysis(
  expr_file = expr_file,
  clinical_file = clinical_file,
  phenotype_col = "response",
  positive_label = "response",
  negative_label = "non_response"
)
```

Drug sensitivity:

```r
run_pair_marker_analysis(
  expr_file = expr_file,
  clinical_file = clinical_file,
  phenotype_col = "drug_status",
  positive_label = "sensitive",
  negative_label = "resistant"
)
```

Mutation status:

```r
run_pair_marker_analysis(
  expr_file = expr_file,
  clinical_file = clinical_file,
  phenotype_col = "TP53_status",
  positive_label = "mutant",
  negative_label = "wildtype"
)
```

## Method Overview

`PairMarker` uses the following workflow:

1. Read an expression matrix and a clinical annotation table.
2. Remove duplicated or missing gene names.
3. Filter genes by the proportion of samples with non-zero expression.
4. Match samples between the expression matrix and clinical table.
5. Convert the selected phenotype labels to a binary vector, where
   `positive_label` is coded as 1 and `negative_label` is coded as 0.
6. Apply `log2(expression + 1)` transformation internally.
7. For each gene pair, compare the two genes within each sample. If the
   expression difference is larger than `delta`, the pair is assigned to one of
   two relative-expression states.
8. Use a chi-square test as a fast screening step.
9. Apply Fisher's exact test to screened pairs.
10. Adjust Fisher p-values using Benjamini-Hochberg correction.

In this context, a significant gene pair means that the relative expression
pattern between two genes differs between the positive and negative phenotype
groups.

## Input Files

The expression file should be a tab-delimited text file. The first column should
contain gene names, and the remaining columns should contain samples.

```text
gene    S1    S2    S3
GENE_A  12    11    10
GENE_B  2     1     2
```

The clinical file should be a tab-delimited text file. The first column should
contain sample IDs, and one column should contain the phenotype label.

```text
sample  response      OS_time  OS_status  PFS_time  PFS_status
S1      response      36       0          18        0
S2      non_response  12       1          5         1
```

Sample IDs in the clinical file should match sample names in the expression
matrix.

Important input requirements:

- Rows of the expression file represent genes and columns represent samples.
- Expression values should be numeric raw or normalized expression values. The
  package applies `log2(expression + 1)` internally before pairwise comparison.
- Duplicated gene names are removed, keeping the first occurrence.
- Genes with too many zero-expression samples are removed according to
  `min_nonzero_prop`.
- Genes with missing expression values after filtering are removed.
- The first column of the clinical file is used as sample IDs.
- The value supplied to `positive_label` must appear in the `phenotype_col`
  column.
- If supplied, the value of `negative_label` must also appear in
  `phenotype_col`.
- Survival columns are optional. If provided, time columns such as `OS_time` or
  `PFS_time` should be numeric, and event columns such as `OS_status` or
  `PFS_status` should use 1 for event and 0 for censored.

## Parameters

Common parameters in `run_pair_marker_analysis()`:

| Parameter | Description | Default |
| --- | --- | --- |
| `expr_file` | Path to the tab-delimited expression matrix. | Required |
| `clinical_file` | Path to the tab-delimited clinical annotation file. | Required |
| `out_dir` | Output directory. If `NULL`, result files are not written. | `NULL` |
| `gene_col` | Gene-name column in the expression file. If `NULL`, the package detects common names such as `gene`, `Gene`, or uses the first column. | `NULL` |
| `phenotype_col` | Clinical column containing phenotype labels. | `"response"` |
| `positive_label` | Label coded as the positive group, 1. | `"response"` |
| `negative_label` | Optional label coded as the negative group, 0. Other labels are excluded when this is supplied. | `NULL` |
| `response_col` | Deprecated alias for `phenotype_col`. | `NULL` |
| `response_label` | Deprecated alias for `positive_label`. | `NULL` |
| `main_delta` | Delta cutoff used in the main pairwise analysis. Larger values require a stronger expression difference between two genes. | `0.25` |
| `delta_list` | Delta cutoffs used for sensitivity analysis. | `c(0, 0.25, 0.5)` |
| `dataset_name` | Prefix used for written output file names. | Expression file name |
| `sig_cutoff` | Adjusted p-value cutoff for significant gene pairs. | `0.05` |
| `min_nonzero_prop` | Minimum proportion of samples with expression greater than zero for keeping a gene. | `0.5` |
| `min_prop` | Minimum proportion allowed for one relative-expression state in a gene pair. | `0.05` |
| `max_prop` | Maximum proportion allowed for one relative-expression state in a gene pair. | `0.95` |
| `chisq_cutoff` | Chi-square p-value cutoff used for screening. | `0.01` |
| `min_valid_prop` | Minimum proportion of samples with a valid pairwise comparison for a gene pair. | `0.5` |

## Main Output

`run_pair_marker_analysis()` returns a list containing:

- `expr`: filtered and sample-aligned expression matrix
- `clinical`: sample-aligned clinical annotation table
- `phenotype`: binary phenotype vector used in the analysis
- `response`: deprecated alias of `phenotype`, kept for backward compatibility
- `sensitivity`: number of gene pairs passing chi-square screening under each
  delta cutoff
- `pairs`: all screened gene pairs with chi-square p-values, Fisher p-values,
  odds ratios, and adjusted p-values
- `sig_pairs`: significant gene pairs after adjusted p-value filtering
- `output_paths`: paths of written result files

The `pairs` and `sig_pairs` tables contain these main columns:

| Column | Description |
| --- | --- |
| `gene1`, `gene2` | Gene pair tested. |
| `nr0`, `nr1` | Counts of negative-class samples in relative-expression state 0 or 1. |
| `r0`, `r1` | Counts of positive-class samples in relative-expression state 0 or 1. |
| `chisq_p` | Chi-square p-value used in the screening step. |
| `OR` | Odds ratio from Fisher's exact test. |
| `fisher_p` | Fisher's exact test p-value. |
| `adjusted_p` | Benjamini-Hochberg adjusted Fisher p-value. |

When `out_dir` is provided, the workflow writes:

- `<dataset_name>_delta_sensitivity_summary.txt`
- `<dataset_name>_chisq_fisher_bh.txt`
- `<dataset_name>_sig_pairs_adjP<sig_cutoff>.txt`, if significant pairs are
  found

## Multi-dataset Integration

After running pair-marker screening in several independent datasets,
`PairMarker` can optionally integrate the result tables to identify recurrent
and directionally consistent gene pairs. This step is useful when a single
dataset produces many significant pairs and you want to prioritize pairs that
are repeatedly observed across cohorts.

The integration step uses:

- how many datasets contain the pair as significant
- which datasets support the pair
- the minimum adjusted p-value across datasets
- the mean and median log2 odds ratio
- OR-direction consistency across datasets
- an evidence score based on `-log10(adjusted_p)`

Example using result files written by `run_pair_marker_analysis()`:

```r
pair_files <- c(
  GSE91061 = "GSE91061_chisq_fisher_bh.txt",
  GSE78220 = "GSE78220_chisq_fisher_bh.txt",
  PRJEB23709 = "PRJEB23709_chisq_fisher_bh.txt"
)

integrated <- integrate_pair_results(
  pair_results = pair_files,
  fdr_cutoff = 0.05,
  min_datasets = 2,
  min_direction_consistency = 0.67
)

head(integrated$summary)
```

The same function also accepts in-memory result tables:

```r
integrated <- integrate_pair_results(
  list(
    cohort1 = result1$pairs,
    cohort2 = result2$pairs,
    cohort3 = result3$pairs
  ),
  min_datasets = 2
)
```

To write the integrated table:

```r
write_integrated_pair_results(
  integrated,
  file = "stable_gene_pairs.tsv",
  include_all = TRUE
)
```

The main output is `integrated$summary`, where higher-priority pairs are ranked
by more supporting datasets, stronger evidence scores, smaller adjusted
p-values, and more consistent OR directions. When `canonicalize_pairs = TRUE`,
the default, pairs such as `GENE_A|GENE_B` and `GENE_B|GENE_A` are treated as
the same pair and the OR direction is harmonized.

## Downstream Plots

`PairMarker` provides base R plotting functions for quick downstream inspection:

- `plot_pair_roc()`: ROC curve for one or more gene pairs, using the continuous
  expression difference `gene1 - gene2` by default
- `plot_pair_survival()`: Kaplan-Meier curve when survival time and event
  columns are available

Both plotting functions draw to the active graphics device by default. To save
a plot directly, provide a file path:

```r
plot_pair_roc(
  result,
  gene1 = "GENE_A",
  gene2 = "GENE_B",
  file = "pair_roc.pdf"
)

plot_pair_survival(
  result,
  gene1 = "GENE_A",
  gene2 = "GENE_B",
  time_col = "OS_time",
  event_col = "OS_status",
  file = "pair_os_survival.pdf"
)
```

## Result Interpretation

Each gene pair is converted into a relative-expression comparison within each
sample. For example, state 1 means that `gene1` is higher than `gene2` by more
than the selected delta cutoff, while state 0 means that `gene1` is lower than
`gene2` by more than the cutoff. Samples with smaller differences are excluded
from that pairwise comparison.

A significant adjusted p-value indicates that the distribution of the two
relative-expression states is different between positive and negative phenotype
groups. Such pairs may be candidate pairwise markers, but they should be
validated in independent datasets before being treated as robust biomarkers.

## Citation

If you use `PairMarker` in your research, please cite the GitHub repository:

```text
https://github.com/twl-00/PairMarker
```

## License

This package is licensed under the MIT License.
