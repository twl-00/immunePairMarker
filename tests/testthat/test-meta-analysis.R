test_that("integrate_pair_results summarizes recurrent significant pairs", {
  ds1 <- data.frame(
    gene1 = c("A", "C"),
    gene2 = c("B", "D"),
    OR = c(4, 2),
    fisher_p = c(0.001, 0.02),
    adjusted_p = c(0.01, 0.2)
  )
  ds2 <- data.frame(
    gene1 = c("A", "B"),
    gene2 = c("B", "C"),
    OR = c(8, 3),
    fisher_p = c(0.002, 0.03),
    adjusted_p = c(0.02, 0.04)
  )

  integrated <- integrate_pair_results(
    list(ds1 = ds1, ds2 = ds2),
    fdr_cutoff = 0.05,
    min_datasets = 2
  )

  expect_equal(nrow(integrated$summary), 1)
  expect_equal(integrated$summary$pair_id, "A|B")
  expect_equal(integrated$summary$n_dataset, 2)
  expect_equal(integrated$summary$direction_consistency, 1)
})

test_that("integrate_pair_results supports legacy single pair column", {
  ds1 <- data.frame(
    gene = c("A|B", "C|D"),
    OR = c(4, 2),
    fisher_p = c(0.001, 0.02),
    adjusted_p = c(0.01, 0.2)
  )
  ds2 <- data.frame(
    gene = c("A|B", "B|C"),
    OR = c(8, 3),
    fisher_p = c(0.002, 0.03),
    adjusted_p = c(0.02, 0.04)
  )

  integrated <- integrate_pair_results(
    list(ds1 = ds1, ds2 = ds2),
    fdr_cutoff = 0.05,
    min_datasets = 2
  )

  expect_equal(integrated$summary$pair_id, "A|B")
  expect_equal(integrated$summary$gene1, "A")
  expect_equal(integrated$summary$gene2, "B")
})

test_that("canonical pair integration flips reversed pair direction", {
  ds1 <- data.frame(
    gene1 = "A",
    gene2 = "B",
    OR = 4,
    fisher_p = 0.001,
    adjusted_p = 0.01
  )
  ds2 <- data.frame(
    gene1 = "B",
    gene2 = "A",
    OR = 0.25,
    fisher_p = 0.002,
    adjusted_p = 0.02
  )

  integrated <- integrate_pair_results(
    list(ds1 = ds1, ds2 = ds2),
    fdr_cutoff = 0.05,
    min_datasets = 2,
    min_direction_consistency = 1
  )

  expect_equal(integrated$summary$pair_id, "A|B")
  expect_equal(integrated$summary$direction_consistency, 1)
  expect_equal(integrated$summary$mean_logOR, 2)
})
