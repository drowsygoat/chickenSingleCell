#' Plot Gene Expression per Cluster across Samples using Patchwork
#'
#' Creates boxplots (with quasirandom jitter) of gene expression per cluster across samples.
#' Uses patchwork to display all plots in a single layout.
#'
#' @param expr_list Output from AverageExpressionPerClusterPerSample().
#' @param long_df Optional precomputed long-format data frame (overrides expr_list).
#' @param genes Character vector of gene names to include. If NULL, uses all.
#' @param base_width Numeric, base width per sample (inches). Default: 1.5
#' @param base_height Numeric, base height per cluster (inches). Default: 3
#' @param save_as_pdf Logical, if TRUE saves the plots as PDF. Default: FALSE
#' @param subset Numeric between 0 and 1, to randomly downsample data (optional).
#' @param threads Number of parallel threads to use for plotting. Default: 1
#'
#' @return A named list containing:
#'   - linearplots: list of ggplot objects (linear scale)
#'   - logplots: list of ggplot objects (log scale)
#'   - detection_plot_clusters: ggplot object
#'   - detection_plot_samples: ggplot object
#' @export
PlotBulkExpressionViolins <- function(expr_list = NULL,
                                      long_df = NULL,
                                      genes = NULL,
                                      base_width = 1.5,
                                      base_height = 3,
                                      save_as_pdf = FALSE,
                                      subset = NULL,
                                      threads = 1,
                                      y_max = NULL) {

  required_pkgs <- c("ggplot2", "tidyr", "dplyr", "purrr")

  lapply(required_pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
    }
  })

  if (is.null(long_df)) {
    if (length(expr_list) == 0) stop("Either expr_list or long_df must be provided.")

    # Combine into long format
    long_df <- purrr::map_dfr(names(expr_list), function(sample_name) {
      mat <- expr_list[[sample_name]]
      if (is.null(mat)) return(NULL)
      mat <- as.data.frame(mat)
      mat$gene <- rownames(mat)
      long <- tidyr::pivot_longer(mat, -gene, names_to = "cluster", values_to = "expression") |>
        dplyr::filter(expression > 0)
      long$sample <- gsub("_.*$", "", sample_name)
      long$assay <- sub("^.*_", "", sample_name)
      return(long)
    })
  }

  if (!is.null(genes)) {
    long_df <- dplyr::filter(long_df, gene %in% genes)
  }

  # Clean sample names
  long_df$sample <- gsub("^sample_|_obj[0-9]+$", "", long_df$sample)

  # Optional subset
  if (!is.null(subset)) {
    if (!is.numeric(subset) || subset <= 0 || subset > 1) {
      stop("subset must be a number between 0 and 1.")
    }
    set.seed(42)
    long_df <- dplyr::slice_sample(long_df, prop = subset)
  }

  # Split by cluster
  split_data <- long_df |>
    split(~cluster) |>
    purrr::imap(function(df, cluster_name) {
      before <- dplyr::n_distinct(df$gene)
      df <- dplyr::filter(df, expression > 0)
      after <- dplyr::n_distinct(df$gene)
      dropped <- before - after
      message(sprintf(
        "Cluster %s: dropped %d gene(s) with all-zero expression. Total before: %d, after: %d",
        cluster_name, dropped, before, after
      ))
      return(df)
    })

  # Generate plots and save if save_as_pdf = TRUE
  linearplots <- save_cluster_expression_plots(split_data, save_as_pdf = TRUE, log_scale = FALSE, file_prefix = "cluster_expression", threads = threads, y_max = y_max)

  logplots <- save_cluster_expression_plots(split_data, save_as_pdf = TRUE, log_scale = TRUE, file_prefix = "cluster_expression", threads = threads)

  detection_plot_clusters <- PlotCP10KDetectionStats(group_by = "cluster", save_as_pdf = TRUE)

  detection_plot_samples  <- PlotCP10KDetectionStats(group_by = "sample", save_as_pdf = TRUE)

  return(list(
    linearplots = linearplots,
    logplots = logplots,
    detection_plot_clusters = detection_plot_clusters,
    detection_plot_samples = detection_plot_samples
  ))
}
