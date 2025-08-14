#' Plot Gene Expression per Cluster across Samples using Patchwork
#'
#' Creates boxplots (with quasirandom jitter) of gene expression per cluster across samples.
#' Uses patchwork to display all plots in a single layout.
#'
#' @param expr_list Output from AverageExpressionPerClusterPerSample().
#' @param long_df Optional precomputed long-format data frame (overrides expr_list).
#' @param genes Character vector of gene names to include. If NULL, uses all.
#' @param save_as_pdf Logical, if TRUE saves the plots as PDF. Default: FALSE
#' @param subset Numeric between 0 and 1, to randomly downsample data (optional).
#' @param threads Number of parallel threads to use for plotting. Default: 1
#' @param y_max Optional numeric to set the y-axis maximum.
#' @param plot_dir Directory to save plots. Default: "plots"
#'
#' @return A named list of plots
#' @export
PlotBulkExpressionViolins2 <- function(expr_list = NULL,
                                      long_df = NULL,
                                      genes = NULL,
                                      save_as_pdf = FALSE,
                                      subset = NULL,
                                      threads = 1,
                                      y_max = NULL,
                                      plot_dir = "plots") {

  required_pkgs <- c("ggplot2", "tidyr", "dplyr", "purrr")
  invisible(lapply(required_pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
    }
  }))

  # Helper to reshape expression list into long format
  generate_long_df <- function(expr_list, subset_genes = NULL, min_expression = 0) {
    purrr::map_dfr(names(expr_list), function(sample_name) {
      mat <- expr_list[[sample_name]]
      if (is.null(mat)) return(NULL)

      # Subset genes before reshaping
      if (!is.null(subset_genes)) {
        mat <- mat[rownames(mat) %in% subset_genes, , drop = FALSE]
      }

      mat <- as.data.frame(mat)
      mat$gene <- rownames(mat)

      long <- tidyr::pivot_longer(
        mat,
        cols = -gene,
        names_to = "cluster",
        values_to = "expression"
      ) |>
        dplyr::filter(expression > min_expression)

      long$sample <- gsub("_.*$", "", sample_name)
      long$assay <- sub("^.*_", "", sample_name)
      return(long)
    })
  }

  # Prepare long_df
  if (is.null(long_df)) {
    if (length(expr_list) == 0) stop("Either expr_list or long_df must be provided.")
    long_df <- generate_long_df(expr_list, subset_genes = genes)  }

  # Filter genes
  if (!is.null(genes)) {
    long_df <- dplyr::filter(long_df, gene %in% genes)
  }

  # Clean sample names
  long_df$sample <- gsub("^sample_|_obj[0-9]+$", "", long_df$sample)

  # Optional downsampling
  if (!is.null(subset)) {
    if (!is.numeric(subset) || subset <= 0 || subset > 1) {
      stop("subset must be a numeric value between 0 and 1.")
    }
    set.seed(42)
    long_df <- dplyr::slice_sample(long_df, prop = subset)
  }

  # Split by cluster and filter
  split_data <- long_df |>
    split(~cluster) |>
    purrr::imap(function(df, cluster_name) {
      before <- dplyr::n_distinct(df$gene)
      df <- dplyr::filter(df, expression > 0)
      after <- dplyr::n_distinct(df$gene)
      if (before != after) {
        message(sprintf("Cluster %s: dropped %d gene(s) with all-zero expression", cluster_name, before - after))
      }
      df
    })

  # Generate plots

  # comb_facets <- save_combined_cluster_expression_plot(
  #   split_data, plot_dir = plot_dir, log_scale = TRUE,
  #   y_max = y_max, file_prefix = "cluster_expression_comb_facets", save_as_pdf = save_as_pdf, facet_clusters = FALSE, color_by_cluster = FALSE,combine_clusters = TRUE
  # )

  # comb_color <- save_combined_cluster_expression_plot(
  #   split_data, plot_dir = plot_dir, log_scale = TRUE,
  #   y_max = y_max, file_prefix = "cluster_expression_comb_color", save_as_pdf = save_as_pdf, facet_clusters = FALSE, color_by_cluster = TRUE,combine_clusters = TRUE
  # )

  # separated_facets <- save_combined_cluster_expression_plot(
  #   split_data, plot_dir = plot_dir, log_scale = TRUE,
  #   y_max = y_max, file_prefix = "cluster_expression_separated_facets", save_as_pdf = save_as_pdf,  facet_clusters = TRUE, color_by_cluster = FALSE, combine_clusters = FALSE
  # )

  linearplot_comb <- save_cluster_expression_plots(
    split_data, save_as_pdf = save_as_pdf, log_scale = TRUE,
    file_prefix = "cluster_expression_no_dodge", threads = threads,
    y_max = y_max, plot_dir = plot_dir, combine_clusters = TRUE, fill_by_clusters = FALSE
  )

  logplot_comb <- save_cluster_expression_plots(
    split_data, save_as_pdf = save_as_pdf, log_scale = TRUE,
    file_prefix = "cluster_expression_dodge", threads = threads,
    y_max = y_max, plot_dir = plot_dir, combine_clusters = TRUE, fill_by_clusters = TRUE, point = 0.01
  )

  linearplots <- save_cluster_expression_plots(
    split_data, save_as_pdf = save_as_pdf, log_scale = FALSE,
    file_prefix = "cluster_expression", threads = threads,
    y_max = y_max, plot_dir = plot_dir
  )

  logplots <- save_cluster_expression_plots(
    split_data, save_as_pdf = save_as_pdf, log_scale = TRUE,
    file_prefix = "cluster_expression", threads = threads,
    y_max = y_max, plot_dir = plot_dir
  )

  detection_plot_clusters <- PlotCP10KDetectionStats(
    group_by = "cluster", save_as_pdf = save_as_pdf, plot_dir = plot_dir
  )

  detection_plot_samples <- PlotCP10KDetectionStats(
    group_by = "sample", save_as_pdf = save_as_pdf, plot_dir = plot_dir
  )

  return(list(
    linearplots = linearplots,
    logplots = logplots,
    logplot_comb = logplot_comb,
    linearplot_comb = linearplot_comb,
    detection_plot_clusters = detection_plot_clusters,
    detection_plot_samples = detection_plot_samples
  ))
}
