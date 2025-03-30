#' Plot Gene Expression per Cluster across Samples using Patchwork
#'
#' Creates boxplots (with quasirandom jitter) of gene expression per cluster across samples.
#' Uses patchwork to display all plots in a single layout.
#'
#' @param expr_list Output from AverageExpressionPerClusterPerSample().
#' @param genes Character vector of gene names to include. If NULL, uses all.
#' @param base_width Numeric, base width per sample (inches). Default: 1.5
#' @param base_height Numeric, base height per cluster (inches). Default: 3
#' @param plot Logical, if TRUE saves the combined plot as PDF. Default: FALSE
#'
#' @return A patchwork object of all cluster plots. Width and height attributes are attached for saving.
#' @export
#'
#' @examples
#' plots <- PlotBulkExpressionViolins(expr_list, genes = c("Gene1", "Gene2"))
#' ggsave("cluster_plots.pdf", plots, width = attr(plots, "width"), height = attr(plots, "height"))
PlotBulkExpressionViolins <- function(expr_list,
                                      genes = NULL,
                                      base_width = 1.5,
                                      base_height = 3,
                                      plot = FALSE,
                                      subset = NULL,
                                      file_prefix) {
  # Required packages
  required_pkgs <- c("ggplot2", "tidyr", "dplyr", "patchwork", "purrr", "ggbeeswarm")
  lapply(required_pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
    }
  })

  if (length(expr_list) == 0) stop("Empty input list")

  # Combine matrices into long format
  long_df <- purrr::map_dfr(names(expr_list), function(sample_name) {
    mat <- expr_list[[sample_name]]
    if (is.null(mat)) return(NULL)
    mat <- as.data.frame(mat)
    mat$gene <- rownames(mat)
    long <- tidyr::pivot_longer(mat, -gene, names_to = "cluster", values_to = "expression") |> dplyr::filter(expression > 0)

    long$sample <- gsub("_.*$", "", sample_name)

    # Extract assay from sample name (assumes format: anything_<assay>)
    assay <- sub("^.*_", "", sample_name)
    long$assay <- assay
    return(long)
  })

  if (!is.null(genes)) {
    long_df <- dplyr::filter(long_df, gene %in% genes)
  }

  # Clean sample names (remove prefixes/suffixes)
  long_df$sample <- gsub("^sample_|_obj[0-9]+$", "", long_df$sample)

  # Subset the data if a fraction is given
  if (!is.null(subset)) {
    if (!is.numeric(subset) || subset <= 0 || subset > 1) {
      stop("subset must be a number between 0 and 1.")
    }
    set.seed(42)  # for reproducibility
    long_df <- dplyr::slice_sample(long_df, prop = subset)
  }

  # Save the long_df for future reference
  if (!dir.exists("data")) dir.create("data")
  saveRDS(long_df, file = file.path("data", "long_df_cluster_expression.rds"))

  split_data <- long_df |>
  split(~cluster) |>
  purrr::imap(function(df, cluster_name) {
    before <- dplyr::n_distinct(df$gene)
    df <- df |>
      # dplyr::group_by(gene) |>
      # dplyr::filter(sum(expression > 0) > 0) |>
      dplyr::filter(expression > 0) # |>
      # dplyr::ungroup()
    after <- dplyr::n_distinct(df$gene)
    dropped <- before - after
    # if (dropped > 0) {
      message(sprintf(
        "Cluster %s: dropped %d gene(s) with all-zero expression. Total before: %d, after: %d",
        cluster_name, dropped, before, after
      ))
    # }
    return(df)
  })

  # Optional saving
  if (plot) {
    save_cluster_expression_plots(split_data, plot = TRUE, log_scale = FALSE, file_prefix = "cluster_expression")
    
    save_cluster_expression_plots(split_data, plot = TRUE, log_scale = TRUE,  file_prefix = "cluster_expression")

    # Optionally generate CP10K detection stats plots
    detection_plot_clusters <- PlotCP10KDetectionStats()
    detection_plot_samples <- PlotCP10KDetectionStatsPerSample()

    detection_file_clusters <- file.path("plots", "CP10K_detection_stats_clusters.pdf")

    detection_file_samples  <- file.path("plots", "CP10K_detection_stats_samples.pdf")

    ggplot2::ggsave(detection_file_clusters, detection_plot_clusters, width = 10, height = 6)

    message(sprintf("CP10K detection plot per cluster saved to: %s", detection_file_clusters))

    ggplot2::ggsave(detection_file_samples, detection_plot_samples, width = 10, height = 6)

    message(sprintf("CP10K detection plot per sample saved to: %s", detection_file_samples))
  }

  return(combined)
}