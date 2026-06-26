# PairMarker

Pairwise gene marker discovery for binary phenotypes in R.

`PairMarker` screens gene pairs using within-sample relative expression
comparisons, chi-square screening, Fisher's exact test, and multiple-testing
correction. It was developed for immunotherapy response analysis, but the same
workflow can be used for any matched expression matrix and binary annotation,
such as responder versus non-responder, sensitive versus resistant, tumor
versus normal, or mutant versus wildtype.

![PairMarker workflow overview](man/figures/PairMarker-workflow.png)

## Installation

```r
install.packages("remotes")
remotes::install_github("twl-00/PairMarker", upgrade = "never")
```

## Quick Start

```r
library(PairMarker)

expr_file <- system.file("extdata", "example_exp.txt", package = "PairMarker")
clinical_file <- system.file("extdata", "example_clinical.txt", package = "PairMarker")

result <- run_pair_marker_analysis(
  expr_file = expr_file,
  clinical_file = clinical_file,
  out_dir = tempdir(),
  phenotype_col = "response",
  positive_label = "response",
  negative_label = "non_response",
  dataset_name = "example"
)

head(result$sig_pairs)
```

## Documentation

- [User guide](docs/user-guide.md): phenotype definitions, input formats,
  parameters, outputs, multi-dataset integration, plotting, and interpretation.
- [Workflow PDF](man/figures/PairMarker-workflow-editable.pdf): editable
  workflow overview.
- [Function reference](man/): package help pages generated from R documentation.

## Main Functions

- `run_pair_marker_analysis()`: run pairwise marker screening for one dataset.
- `integrate_pair_results()`: prioritize recurrent pairs across datasets.
- `plot_pair_roc()` and `plot_pair_survival()`: inspect selected gene pairs.
- `plot_integrated_pair_roc()` and `plot_integrated_pair_survival()`: inspect
  prioritized integrated pairs.

## Citation

If you use `PairMarker` in your research, please cite this repository:

```text
https://github.com/twl-00/PairMarker
```

## License

MIT.
