#' Read one pair-marker result table for multi-dataset integration
#'
#' @param path Path to a tab-delimited pair result table.
#' @param dataset Dataset name used in the integrated result.
#' @param fdr_cutoff Adjusted p-value cutoff used to mark significant pairs.
#' @param pair_col Optional column containing a pair identifier. If `NULL`, the
#'   function uses `gene1` and `gene2` when available.
#' @param gene1_col,gene2_col Columns containing the two genes in each pair.
#' @param or_col Column containing the odds ratio.
#' @param p_col Column containing Fisher's exact test p-value.
#' @param adjusted_p_col Column containing adjusted p-values. If absent, values
#'   are calculated from `p_col`.
#' @return A data frame with one row per gene pair in one dataset.
#' @details
#' A minimal input table should contain `gene1`, `gene2`, `OR`, `fisher_p`, and
#' `adjusted_p`. Older tables with a single pair column such as
#' `gene = "CD8A|PDCD1"` are also supported.
#' @export
read_pair_result <- function(
    path,
    dataset = tools::file_path_sans_ext(basename(path)),
    fdr_cutoff = 0.05,
    pair_col = NULL,
    gene1_col = "gene1",
    gene2_col = "gene2",
    or_col = "OR",
    p_col = "fisher_p",
    adjusted_p_col = "adjusted_p") {

  x <- utils::read.table(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  .standardize_pair_result(
    x,
    dataset = dataset,
    fdr_cutoff = fdr_cutoff,
    pair_col = pair_col,
    gene1_col = gene1_col,
    gene2_col = gene2_col,
    or_col = or_col,
    p_col = p_col,
    adjusted_p_col = adjusted_p_col
  )
}

#' Integrate significant gene pairs across datasets
#'
#' @param pair_results A named list of data frames, or a named character vector
#'   of file paths. Names are used as dataset names.
#' @param fdr_cutoff Adjusted p-value cutoff used to define significant pairs.
#' @param min_datasets Minimum number of datasets in which a pair must be
#'   significant.
#' @param min_direction_consistency Minimum OR-direction consistency across
#'   significant occurrences. A value of 1 requires all significant occurrences
#'   to have the same direction.
#' @param pair_col Optional column containing a pair identifier. If `NULL`, the
#'   function uses `gene1` and `gene2` when available.
#' @param gene1_col,gene2_col Columns containing the two genes in each pair.
#' @param or_col Column containing the odds ratio.
#' @param p_col Column containing Fisher's exact test p-value.
#' @param adjusted_p_col Column containing adjusted p-values. If absent, values
#'   are calculated from `p_col`.
#' @param canonicalize_pairs If `TRUE`, gene-pair identifiers are sorted so
#'   `A|B` and `B|A` are treated as the same pair.
#' @return A list with `summary`, `all_pairs`, and `significant_pairs`.
#' @examples
#' cohort1_pairs <- data.frame(
#'   gene1 = c("CD8A", "CXCL9", "MKI67"),
#'   gene2 = c("PDCD1", "LAG3", "TOP2A"),
#'   OR = c(5.2, 3.1, 0.8),
#'   fisher_p = c(0.0004, 0.006, 0.20),
#'   adjusted_p = c(0.004, 0.03, 0.50)
#' )
#' cohort2_pairs <- data.frame(
#'   gene1 = c("CD8A", "CXCL9", "GZMB"),
#'   gene2 = c("PDCD1", "LAG3", "IFNG"),
#'   OR = c(4.4, 2.6, 0.7),
#'   fisher_p = c(0.001, 0.008, 0.10),
#'   adjusted_p = c(0.006, 0.04, 0.30)
#' )
#' cohort3_pairs <- data.frame(
#'   gene1 = c("CD8A", "CXCL10", "MKI67"),
#'   gene2 = c("PDCD1", "TIGIT", "TOP2A"),
#'   OR = c(6.0, 2.2, 1.3),
#'   fisher_p = c(0.0002, 0.02, 0.09),
#'   adjusted_p = c(0.003, 0.08, 0.25)
#' )
#'
#' integrated <- integrate_pair_results(
#'   list(cohort1 = cohort1_pairs, cohort2 = cohort2_pairs, cohort3 = cohort3_pairs),
#'   min_datasets = 2
#' )
#' integrated$summary
#' @export
integrate_pair_results <- function(
    pair_results,
    fdr_cutoff = 0.05,
    min_datasets = 2,
    min_direction_consistency = 0.67,
    pair_col = NULL,
    gene1_col = "gene1",
    gene2_col = "gene2",
    or_col = "OR",
    p_col = "fisher_p",
    adjusted_p_col = "adjusted_p",
    canonicalize_pairs = TRUE) {

  if (is.null(pair_results) || length(pair_results) == 0) {
    stop("pair_results must contain at least one dataset.", call. = FALSE)
  }
  if (min_datasets < 1) {
    stop("min_datasets must be at least 1.", call. = FALSE)
  }
  if (min_direction_consistency < 0 || min_direction_consistency > 1) {
    stop("min_direction_consistency must be between 0 and 1.", call. = FALSE)
  }

  dataset_names <- names(pair_results)
  if (is.null(dataset_names)) {
    dataset_names <- rep("", length(pair_results))
  }
  missing_names <- is.na(dataset_names) | dataset_names == ""
  dataset_names[missing_names] <- paste0("dataset", which(missing_names))

  rows <- vector("list", length(pair_results))
  for (i in seq_along(pair_results)) {
    item <- pair_results[[i]]
    dataset <- dataset_names[i]

    if (is.character(item) && length(item) == 1) {
      rows[[i]] <- read_pair_result(
        path = item,
        dataset = dataset,
        fdr_cutoff = fdr_cutoff,
        pair_col = pair_col,
        gene1_col = gene1_col,
        gene2_col = gene2_col,
        or_col = or_col,
        p_col = p_col,
        adjusted_p_col = adjusted_p_col
      )
    } else if (is.data.frame(item)) {
      rows[[i]] <- .standardize_pair_result(
        item,
        dataset = dataset,
        fdr_cutoff = fdr_cutoff,
        pair_col = pair_col,
        gene1_col = gene1_col,
        gene2_col = gene2_col,
        or_col = or_col,
        p_col = p_col,
        adjusted_p_col = adjusted_p_col
      )
    } else {
      stop("Each pair_results element must be a file path or data frame.", call. = FALSE)
    }
  }

  all_pairs <- do.call(rbind, rows)
  rownames(all_pairs) <- NULL

  if (canonicalize_pairs) {
    all_pairs <- .canonicalize_pair_result_rows(all_pairs)
  }

  significant_pairs <- all_pairs[all_pairs$sig, , drop = FALSE]
  if (nrow(significant_pairs) == 0) {
    summary <- .empty_integrated_summary()
  } else {
    significant_pairs <- significant_pairs[order(
      significant_pairs$pair_id,
      significant_pairs$dataset,
      significant_pairs$adjusted_p
    ), , drop = FALSE]
    significant_pairs <- significant_pairs[!duplicated(
      paste(significant_pairs$dataset, significant_pairs$pair_id, sep = "\r")
    ), , drop = FALSE]

    split_pairs <- split(significant_pairs, significant_pairs$pair_id)
    summary <- do.call(rbind, lapply(names(split_pairs), function(pair_id) {
      .summarize_integrated_pair(pair_id, split_pairs[[pair_id]])
    }))
    rownames(summary) <- NULL

    keep <- summary$n_dataset >= min_datasets &
      summary$direction_consistency >= min_direction_consistency
    keep[is.na(keep)] <- FALSE
    summary <- summary[keep, , drop = FALSE]
    if (nrow(summary) > 0) {
      summary <- summary[order(
        -summary$n_dataset,
        -summary$evidence_score,
        summary$min_adjusted_p,
        summary$pair_id
      ), , drop = FALSE]
      rownames(summary) <- NULL
    }
  }

  list(
    summary = summary,
    all_pairs = all_pairs,
    significant_pairs = significant_pairs
  )
}

#' Write integrated pair-marker results
#'
#' @param integrated Result returned by [integrate_pair_results()].
#' @param file Output path for the summary table.
#' @param include_all If `TRUE`, also write all standardized pair results and
#'   significant pair occurrences using file-name suffixes.
#' @return Invisibly returns `integrated`.
#' @export
write_integrated_pair_results <- function(integrated, file, include_all = FALSE) {
  if (!is.list(integrated) || is.null(integrated$summary)) {
    stop("integrated must be a result returned by integrate_pair_results().", call. = FALSE)
  }

  utils::write.table(
    integrated$summary,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  if (include_all) {
    prefix <- tools::file_path_sans_ext(file)
    utils::write.table(
      integrated$all_pairs,
      file = paste0(prefix, "_all_pairs.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    utils::write.table(
      integrated$significant_pairs,
      file = paste0(prefix, "_significant_occurrences.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  invisible(integrated)
}

.standardize_pair_result <- function(
    x,
    dataset,
    fdr_cutoff,
    pair_col,
    gene1_col,
    gene2_col,
    or_col,
    p_col,
    adjusted_p_col) {

  if (!is.data.frame(x)) {
    stop("Pair result must be a data frame.", call. = FALSE)
  }
  if (!or_col %in% names(x)) {
    stop("Pair result is missing OR column: ", or_col, call. = FALSE)
  }
  if (!p_col %in% names(x)) {
    stop("Pair result is missing p-value column: ", p_col, call. = FALSE)
  }

  if (!is.null(pair_col) && pair_col %in% names(x)) {
    pair_id <- as.character(x[[pair_col]])
    genes <- .split_pair_id(pair_id)
    gene1 <- genes$gene1
    gene2 <- genes$gene2
  } else if (
    "gene" %in% names(x) &&
      !(gene1_col %in% names(x)) &&
      !(gene2_col %in% names(x))
  ) {
    pair_id <- as.character(x[["gene"]])
    genes <- .split_pair_id(pair_id)
    gene1 <- genes$gene1
    gene2 <- genes$gene2
  } else {
    if (!gene1_col %in% names(x) || !gene2_col %in% names(x)) {
      stop(
        "Pair result must contain gene1/gene2 columns or a pair identifier column.",
        call. = FALSE
      )
    }
    gene1 <- as.character(x[[gene1_col]])
    gene2 <- as.character(x[[gene2_col]])
    pair_id <- paste(gene1, gene2, sep = "|")
  }

  or_value <- suppressWarnings(as.numeric(x[[or_col]]))
  fisher_p <- suppressWarnings(as.numeric(x[[p_col]]))
  if (adjusted_p_col %in% names(x)) {
    adjusted_p <- suppressWarnings(as.numeric(x[[adjusted_p_col]]))
  } else {
    adjusted_p <- stats::p.adjust(fisher_p, method = "BH")
  }

  or_cap <- pmin(pmax(or_value, 1e-6), 1e6)
  log_or <- log2(or_cap)
  direction <- ifelse(log_or >= 0, 1L, -1L)
  sig <- !is.na(adjusted_p) & adjusted_p < fdr_cutoff

  data.frame(
    dataset = dataset,
    pair_id = pair_id,
    gene1 = gene1,
    gene2 = gene2,
    OR = or_value,
    logOR = log_or,
    direction = direction,
    fisher_p = fisher_p,
    adjusted_p = adjusted_p,
    sig = sig,
    stringsAsFactors = FALSE
  )
}

.split_pair_id <- function(pair_id) {
  pair_id <- as.character(pair_id)
  parts <- strsplit(pair_id, "\\s*(?:_vs_|__|--|[|/:,;~_])\\s*", perl = TRUE)
  gene1 <- vapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))
  gene2 <- vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
  list(gene1 = gene1, gene2 = gene2)
}

.canonicalize_pair_result_rows <- function(x) {
  missing_genes <- is.na(x$gene1) | x$gene1 == "" | is.na(x$gene2) | x$gene2 == ""
  swap <- !missing_genes & x$gene1 > x$gene2
  if (any(swap)) {
    old_gene1 <- x$gene1[swap]
    x$gene1[swap] <- x$gene2[swap]
    x$gene2[swap] <- old_gene1

    swap_idx <- which(swap)
    x$OR[swap_idx] <- 1 / x$OR[swap_idx]
    x$logOR[swap] <- -x$logOR[swap]
    x$direction[swap] <- -x$direction[swap]
  }

  x$pair_id <- ifelse(
    missing_genes,
    x$pair_id,
    paste(x$gene1, x$gene2, sep = "|")
  )
  x
}

.summarize_integrated_pair <- function(pair_id, x) {
  valid_direction <- !is.na(x$direction)
  direction_sum <- sum(x$direction[valid_direction], na.rm = TRUE)
  direction_n <- sum(valid_direction)
  direction_consistency <- if (direction_n == 0) NA_real_ else abs(direction_sum) / direction_n
  evidence <- sum(-log10(pmax(x$adjusted_p, 1e-300)), na.rm = TRUE)
  data.frame(
    pair_id = pair_id,
    gene1 = x$gene1[1],
    gene2 = x$gene2[1],
    n_dataset = length(unique(x$dataset)),
    datasets = paste(sort(unique(x$dataset)), collapse = ","),
    min_adjusted_p = min(x$adjusted_p, na.rm = TRUE),
    mean_logOR = mean(x$logOR, na.rm = TRUE),
    median_logOR = stats::median(x$logOR, na.rm = TRUE),
    direction_consistency = direction_consistency,
    evidence_score = evidence,
    stringsAsFactors = FALSE
  )
}

.empty_integrated_summary <- function() {
  data.frame(
    pair_id = character(),
    gene1 = character(),
    gene2 = character(),
    n_dataset = integer(),
    datasets = character(),
    min_adjusted_p = numeric(),
    mean_logOR = numeric(),
    median_logOR = numeric(),
    direction_consistency = numeric(),
    evidence_score = numeric(),
    stringsAsFactors = FALSE
  )
}
