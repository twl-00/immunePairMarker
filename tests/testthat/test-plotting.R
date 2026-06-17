test_that("pair state helper encodes pairwise expression states", {
  expr <- matrix(
    c(10, 10, 1, 1, 1, 1, 10, 10),
    nrow = 2,
    byrow = TRUE
  )
  rownames(expr) <- c("A", "B")
  colnames(expr) <- paste0("S", 1:4)

  pair_state_and_score <- getFromNamespace(".pair_state_and_score", "immunePair")
  pair <- pair_state_and_score(expr, "A", "B", delta = 0.25)

  expect_equal(unname(pair$state), c(1L, 1L, 0L, 0L))
  expect_equal(names(pair$state), colnames(expr))
})

test_that("ROC helper returns high AUC for ordered scores", {
  response <- c(1L, 1L, 0L, 0L)
  score <- c(0.9, 0.8, 0.2, 0.1)

  compute_roc <- getFromNamespace(".compute_roc", "immunePair")
  roc <- compute_roc(response, score)

  expect_equal(roc$auc, 1)
  expect_false(roc$reversed)
})

test_that("HR labels do not print rounded zero", {
  format_hr <- getFromNamespace(".format_hr", "immunePair")

  expect_equal(format_hr(0.0001), "< 0.01")
  expect_equal(format_hr(0.8), "= 0.80")
})

test_that("plot pair selector prefers significant pairs", {
  result <- list(
    pairs = data.frame(
      gene1 = c("A", "C"),
      gene2 = c("B", "D"),
      adjusted_p = c(0.02, 0.03)
    ),
    sig_pairs = data.frame(
      gene1 = "E",
      gene2 = "F",
      adjusted_p = 0.001
    )
  )

  select_pair_table <- getFromNamespace(".select_pair_table", "immunePair")
  selected <- select_pair_table(result, top_n = 1)

  expect_equal(selected$gene1, "E")
  expect_equal(selected$gene2, "F")
})
