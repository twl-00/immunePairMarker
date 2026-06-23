# PairMarker 0.2.0

- Renamed the package from `immunePair` to `PairMarker`.
- Generalized the workflow from immune response labels to binary phenotype labels.
- Added `phenotype_col`, `positive_label`, and `negative_label` arguments while
  keeping `response_col` and `response_label` as deprecated aliases.
- Added optional multi-dataset integration functions for recurrent and
  directionally consistent gene-pair prioritization:
  `integrate_pair_results()`, `read_pair_result()`, and
  `write_integrated_pair_results()`.

# immunePair 0.1.0

Initial public release of `immunePair`.

## Main features

- Added a complete pairwise gene marker analysis workflow.
- Added Rcpp-based pairwise chi-square screening for gene pairs.
- Added Fisher's exact test and adjusted p-value calculation.
- Added delta sensitivity analysis.
- Added example expression and clinical datasets.
- Added runnable examples in the README and function documentation.
- Added basic package tests.
- Passed `R CMD check --no-manual` with `Status: OK`.
