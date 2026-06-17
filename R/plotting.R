#' Plot ROC curves for gene pairs
#'
#' @param result Result list returned by [run_pair_marker_analysis()].
#' @param gene1,gene2 Optional gene names. If omitted, top pairs are selected
#'   from the result table.
#' @param top_n Number of top pairs to plot when `gene1` and `gene2` are not
#'   supplied.
#' @param delta Pairwise expression difference cutoff.
#' @param score Deprecated. ROC analysis always uses the continuous expression
#'   difference `gene1 - gene2`.
#' @param response_col Clinical column containing response labels.
#' @param response_label Label treated as response.
#' @param use_sig_pairs If `TRUE`, use `result$sig_pairs` when available.
#' @param direction If `auto`, scores are reversed when the initial AUC is below
#'   0.5.
#' @param log_transform If `TRUE`, use `log2(expression + 1)` before comparison.
#' @param file Optional output path ending in `.pdf`, `.png`, `.jpg`, or `.jpeg`.
#' @param width Plot width in inches when `file` is provided.
#' @param height Plot height in inches when `file` is provided.
#' @param main Plot title.
#' @return A data frame with AUC values, invisibly.
#' @export
plot_pair_roc <- function(
    result,
    gene1 = NULL,
    gene2 = NULL,
    top_n = 1,
    delta = 0.25,
    score = "difference",
    response_col = "response",
    response_label = "response",
    use_sig_pairs = TRUE,
    direction = c("auto", "as_is"),
    log_transform = TRUE,
    file = NULL,
    width = 6,
    height = 6,
    main = NULL) {

  if (!identical(score, "difference")) {
    warning(
      "score is deprecated; plot_pair_roc() now uses the continuous expression difference.",
      call. = FALSE
    )
  }
  direction <- match.arg(direction)
  expr <- .get_result_expr(result)
  resp <- .get_result_response(
    result,
    response_col = response_col,
    response_label = response_label
  )
  pairs <- .resolve_plot_pairs(
    result,
    gene1 = gene1,
    gene2 = gene2,
    top_n = top_n,
    use_sig_pairs = use_sig_pairs
  )

  close_device <- .open_plot_device(file, width = width, height = height)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    close_device()
  }, add = TRUE)

  single_pair <- nrow(pairs) == 1
  if (is.null(main)) {
    main <- "Pair ROC"
  }

  graphics::plot(
    c(0, 1),
    c(0, 1),
    type = "n",
    xlim = c(-0.02, 1.02),
    ylim = c(-0.02, 1.02),
    xlab = "1 - Specificity",
    ylab = "Sensitivity",
    main = if (single_pair) "" else main,
    xaxs = "i",
    yaxs = "i",
    axes = FALSE
  )
  graphics::segments(0, 0, 1, 1, lty = 2, col = "grey70")
  roc_ticks <- seq(0, 1, by = 0.2)
  graphics::axis(1, at = roc_ticks, labels = sprintf("%.1f", roc_ticks))
  graphics::axis(2, at = roc_ticks, labels = sprintf("%.1f", roc_ticks), las = 1)
  graphics::box()

  colors <- .plot_colors(nrow(pairs))
  auc_rows <- vector("list", nrow(pairs))

  for (i in seq_len(nrow(pairs))) {
    pair_scores <- .pair_state_and_score(
      expr = expr,
      gene1 = pairs$gene1[i],
      gene2 = pairs$gene2[i],
      delta = delta,
      log_transform = log_transform
    )
    marker <- pair_scores$score
    common <- intersect(names(marker), names(resp))
    roc <- .compute_roc(
      response = resp[common],
      score = marker[common],
      direction = direction
    )

    graphics::lines(
      roc$fpr,
      roc$tpr,
      col = colors[i],
      lwd = 3,
      type = "s"
    )
    auc_rows[[i]] <- data.frame(
      gene1 = pairs$gene1[i],
      gene2 = pairs$gene2[i],
      auc = roc$auc,
      n_valid = roc$n_valid,
      reversed = roc$reversed,
      stringsAsFactors = FALSE
    )
  }

  auc_df <- do.call(rbind, auc_rows)
  pair_labels <- paste(auc_df$gene1, auc_df$gene2, sep = " / ")
  legend_labels <- paste0(
    pair_labels,
    " (AUC=",
    sprintf("%.3f", auc_df$auc),
    ")"
  )

  if (single_pair) {
    graphics::title(
      main = paste0(main, "\n", legend_labels[1]),
      line = 1
    )
  }
  graphics::legend(
    "bottomright",
    legend = legend_labels,
    col = colors,
    lwd = 3,
    bty = "n",
    cex = 0.8
  )

  message(
    "ROC AUC: ",
    paste(legend_labels, collapse = "; ")
  )

  invisible(auc_df)
}

#' Plot a Kaplan-Meier curve for a gene pair
#'
#' @param result Result list returned by [run_pair_marker_analysis()].
#' @param gene1,gene2 Optional gene names. If omitted, the top pair is used.
#' @param delta Pairwise expression difference cutoff.
#' @param time_col Clinical column containing survival time.
#' @param event_col Clinical column containing event status, where 1 indicates
#'   event and 0 indicates censored.
#' @param use_sig_pairs If `TRUE`, use `result$sig_pairs` when available.
#' @param log_transform If `TRUE`, use `log2(expression + 1)` before comparison.
#' @param file Optional output path ending in `.pdf`, `.png`, `.jpg`, or `.jpeg`.
#' @param width Plot width in inches when `file` is provided.
#' @param height Plot height in inches when `file` is provided.
#' @param main Plot title.
#' @return A list with the plotted data, survival fit, log-rank p-value, and
#'   Cox hazard ratio, invisibly.
#' @export
plot_pair_survival <- function(
    result,
    gene1 = NULL,
    gene2 = NULL,
    delta = 0.25,
    time_col,
    event_col,
    use_sig_pairs = TRUE,
    log_transform = TRUE,
    file = NULL,
    width = 6,
    height = 6,
    main = NULL) {

  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for survival plots.", call. = FALSE)
  }

  expr <- .get_result_expr(result)
  clinical <- .get_result_clinical(result)
  pairs <- .resolve_plot_pairs(
    result,
    gene1 = gene1,
    gene2 = gene2,
    top_n = 1,
    use_sig_pairs = use_sig_pairs
  )

  missing_cols <- setdiff(c(time_col, event_col), names(clinical))
  if (length(missing_cols) > 0) {
    stop("clinical data is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  pair_scores <- .pair_state_and_score(
    expr = expr,
    gene1 = pairs$gene1[1],
    gene2 = pairs$gene2[1],
    delta = delta,
    log_transform = log_transform
  )

  common <- intersect(names(pair_scores$state), rownames(clinical))
  surv_df <- data.frame(
    time = as.numeric(clinical[common, time_col]),
    event = as.integer(clinical[common, event_col]),
    state = pair_scores$state[common],
    stringsAsFactors = FALSE
  )
  surv_df <- surv_df[stats::complete.cases(surv_df), , drop = FALSE]
  if (nrow(surv_df) == 0 || length(unique(surv_df$state)) < 2) {
    stop("At least two non-missing pair states are required for survival plotting.", call. = FALSE)
  }

  state_labels <- c(
    paste(pairs$gene1[1], "<", pairs$gene2[1]),
    paste(pairs$gene1[1], ">", pairs$gene2[1])
  )
  surv_df$group <- factor(
    ifelse(surv_df$state == 1, state_labels[2], state_labels[1]),
    levels = state_labels
  )

  surv_obj <- survival::Surv(surv_df$time, surv_df$event)
  fit <- survival::survfit(surv_obj ~ group, data = surv_df)
  logrank <- survival::survdiff(surv_obj ~ group, data = surv_df)
  logrank_p <- stats::pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)
  cox_hr <- NA_real_
  cox_fit <- suppressWarnings(try(survival::coxph(surv_obj ~ group, data = surv_df), silent = TRUE))
  if (!inherits(cox_fit, "try-error")) {
    cox_hr <- unname(exp(stats::coef(cox_fit))[1])
  }

  close_device <- .open_plot_device(file, width = width, height = height)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    close_device()
  }, add = TRUE)

  if (is.null(main)) {
    main <- paste0(
      pairs$gene1[1],
      " / ",
      pairs$gene2[1],
      "\nlog-rank p = ",
      format.pval(logrank_p, digits = 3),
      if (is.finite(cox_hr)) paste0(", HR ", .format_hr(cox_hr)) else ""
    )
  }

  graphics::par(mar = .compact_margins(c(4.5, 4.5, 3, 1)))
  graphics::plot(
    fit,
    col = c("#0072B2", "#D55E00"),
    lwd = 2,
    xlab = time_col,
    ylab = "Survival probability",
    main = main,
    mark.time = TRUE
  )
  graphics::legend(
    "bottomleft",
    legend = levels(surv_df$group),
    col = c("#0072B2", "#D55E00"),
    lwd = 2,
    bty = "n"
  )

  invisible(list(
    data = surv_df,
    fit = fit,
    logrank_p = logrank_p,
    cox_hr = cox_hr
  ))
}

.select_pair_table <- function(result, top_n = 10, use_sig_pairs = TRUE) {
  if (is.data.frame(result)) {
    pair_df <- result
  } else if (is.list(result)) {
    if (use_sig_pairs && !is.null(result$sig_pairs) && nrow(result$sig_pairs) > 0) {
      pair_df <- result$sig_pairs
    } else if (!is.null(result$pairs)) {
      pair_df <- result$pairs
    } else {
      stop("result must contain a pairs table.", call. = FALSE)
    }
  } else {
    stop("result must be a result list or a pair result data frame.", call. = FALSE)
  }

  required_cols <- c("gene1", "gene2")
  missing_cols <- setdiff(required_cols, names(pair_df))
  if (length(missing_cols) > 0) {
    stop("pair table is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (nrow(pair_df) == 0) {
    stop("pair table is empty.", call. = FALSE)
  }

  sort_col <- intersect(c("adjusted_p", "fisher_p", "chisq_p"), names(pair_df))
  if (length(sort_col) > 0) {
    pair_df <- pair_df[order(pair_df[[sort_col[1]]]), , drop = FALSE]
  }
  pair_df[seq_len(min(top_n, nrow(pair_df))), , drop = FALSE]
}

.resolve_plot_pairs <- function(result, gene1 = NULL, gene2 = NULL, top_n = 1, use_sig_pairs = TRUE) {
  if (!is.null(gene1) || !is.null(gene2)) {
    if (is.null(gene1) || is.null(gene2)) {
      stop("gene1 and gene2 must be supplied together.", call. = FALSE)
    }
    if (length(gene1) != length(gene2)) {
      stop("gene1 and gene2 must have the same length.", call. = FALSE)
    }
    return(data.frame(gene1 = gene1, gene2 = gene2, stringsAsFactors = FALSE))
  }
  .select_pair_table(result, top_n = top_n, use_sig_pairs = use_sig_pairs)
}

.get_result_expr <- function(result) {
  if (!is.list(result) || is.null(result$expr)) {
    stop("result must be a list containing an expr matrix.", call. = FALSE)
  }
  result$expr
}

.get_result_clinical <- function(result) {
  if (!is.list(result) || is.null(result$clinical)) {
    stop("result must be a list containing clinical data.", call. = FALSE)
  }
  result$clinical
}

.get_result_response <- function(result, response_col = "response", response_label = "response") {
  clinical <- if (is.list(result) && !is.null(result$clinical)) result$clinical else NULL
  if (is.list(result) && !is.null(result$response)) {
    resp <- as.integer(result$response)
    if (is.null(names(resp))) {
      if (!is.null(clinical)) {
        names(resp) <- rownames(clinical)
      } else if (!is.null(result$expr)) {
        names(resp) <- colnames(result$expr)
      }
    }
    return(resp)
  }
  if (is.null(clinical)) {
    stop("result must contain response or clinical data.", call. = FALSE)
  }
  resp <- make_response_vector(
    clinical,
    response_col = response_col,
    response_label = response_label
  )
  names(resp) <- rownames(clinical)
  resp
}

.pair_state_and_score <- function(expr, gene1, gene2, delta = 0.25, log_transform = TRUE) {
  mat <- as.matrix(expr)
  suppressWarnings(storage.mode(mat) <- "numeric")
  missing_genes <- setdiff(c(gene1, gene2), rownames(mat))
  if (length(missing_genes) > 0) {
    stop("gene not found in expression matrix: ", paste(missing_genes, collapse = ", "), call. = FALSE)
  }
  if (log_transform) {
    mat <- log2(mat + 1)
  }

  diff_score <- mat[gene1, ] - mat[gene2, ]
  state <- rep(NA_integer_, length(diff_score))
  state[diff_score > delta] <- 1L
  state[diff_score < -delta] <- 0L
  names(state) <- colnames(mat)
  names(diff_score) <- colnames(mat)
  list(state = state, score = diff_score)
}

.compute_roc <- function(response, score, direction = "auto") {
  keep <- !is.na(response) & !is.na(score)
  response <- as.integer(response[keep])
  score <- as.numeric(score[keep])
  if (length(response) == 0 || length(unique(response)) < 2) {
    stop("ROC analysis requires both response and non-response samples.", call. = FALSE)
  }
  if (length(unique(score)) < 2) {
    stop("ROC analysis requires at least two unique score values.", call. = FALSE)
  }

  roc <- .roc_points(response, score)
  reversed <- FALSE
  if (direction == "auto" && roc$auc < 0.5) {
    roc <- .roc_points(response, -score)
    reversed <- TRUE
  }
  roc$reversed <- reversed
  roc$n_valid <- length(response)
  roc
}

.roc_points <- function(response, score) {
  pos_n <- sum(response == 1)
  neg_n <- sum(response == 0)

  ord <- order(score, decreasing = TRUE)
  response <- response[ord]
  score <- score[ord]

  score_groups <- rle(score)
  threshold_end <- cumsum(score_groups$lengths)
  tp <- cumsum(response == 1)[threshold_end]
  fp <- cumsum(response == 0)[threshold_end]

  fpr <- c(0, fp / neg_n)
  tpr <- c(0, tp / pos_n)
  auc <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)

  list(fpr = fpr, tpr = tpr, auc = auc)
}

.open_plot_device <- function(file = NULL, width = 7, height = 5, res = 150) {
  if (is.null(file)) {
    return(function() invisible(NULL))
  }

  out_dir <- dirname(file)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  ext <- tolower(tools::file_ext(file))
  if (ext == "pdf") {
    grDevices::pdf(file, width = width, height = height)
  } else if (ext == "png") {
    grDevices::png(file, width = width, height = height, units = "in", res = res)
  } else if (ext %in% c("jpg", "jpeg")) {
    grDevices::jpeg(file, width = width, height = height, units = "in", res = res)
  } else {
    stop("Unsupported plot file extension: ", ext, call. = FALSE)
  }
  function() grDevices::dev.off()
}

.plot_colors <- function(n) {
  palette <- c(
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#E69F00",
    "#56B4E9",
    "#000000",
    "#999999"
  )
  rep(palette, length.out = n)
}

.format_hr <- function(hr) {
  if (!is.finite(hr)) {
    return(NA_character_)
  }
  if (hr < 0.01) {
    return("< 0.01")
  }
  paste0("= ", sprintf("%.2f", hr))
}

.compact_margins <- function(mar) {
  din <- tryCatch(graphics::par("din"), error = function(e) c(7, 5))
  line_height <- tryCatch(
    graphics::par("csi") * graphics::par("mex"),
    error = function(e) 0.2
  )
  if (!is.finite(line_height) || line_height <= 0) {
    line_height <- 0.2
  }

  max_width_lines <- max(2.5, din[1] / line_height - 1)
  max_height_lines <- max(2.5, din[2] / line_height - 1)
  min_mar <- c(1.8, 2.2, 1.5, 0.5)

  mar <- pmax(mar, min_mar)
  if (sum(mar[c(2, 4)]) >= max_width_lines) {
    scale <- (max_width_lines * 0.85) / sum(mar[c(2, 4)])
    mar[c(2, 4)] <- pmax(mar[c(2, 4)] * scale, min_mar[c(2, 4)])
  }
  if (sum(mar[c(1, 3)]) >= max_height_lines) {
    scale <- (max_height_lines * 0.85) / sum(mar[c(1, 3)])
    mar[c(1, 3)] <- pmax(mar[c(1, 3)] * scale, min_mar[c(1, 3)])
  }
  mar
}
