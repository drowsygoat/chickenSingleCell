#' Plot Gene Detection Thresholds per Cluster
#'
#' Plots the number of genes detected per cluster across CP10K thresholds.
#' Assumes expression values are already normalized to CP10K in the input RDS.
#'
#' @param long_df_path Path to the saved long-format expression RDS file. Default: "data/long_df_cluster_expression.rds".
#' @param save_plot Logical. If TRUE, saves the plot to "plots/CP10K_detection_stats.pdf".
#'
#' @return A ggplot2 object.
#' @export
PlotCP10KDetectionStats <- function(long_df_path = "data/long_df_cluster_expression.rds", save_plot = FALSE) {
  required_pkgs <- c("dplyr", "tidyr", "ggplot2")
  invisible(lapply(required_pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("Missing required package: %s", pkg))
  }))

  long_df <- readRDS(long_df_path)

  thresholds <- c(
    "≥1 cp10k in ≥1 sample" = function(x) any(x >= 1),
    "≥1 cp10k in all samples" = function(x) all(x >= 1),
    ">5 cp10k in all samples" = function(x) all(x > 5),
    ">10 cp10k in all samples" = function(x) all(x > 10),
    ">50 cp10k in all samples" = function(x) all(x > 50)
  )

  # Precompute CP10K detection flags for each cluster/gene/sample
  summary_df <- long_df |>
    dplyr::group_by(cluster, gene, sample) |>
    dplyr::summarise(expr = dplyr::first(expression), .groups = "drop")

  # Spread into nested list: one data frame per cluster-gene
  detection_df <- summary_df |>
    dplyr::group_by(cluster, gene) |>
    dplyr::summarise(
      across_expr = list(expr),  # list of expression values per gene across samples
      .groups = "drop"
    ) |>
    tidyr::expand_grid(threshold = names(thresholds)) |>
    dplyr::rowwise() |>
    dplyr::mutate(detected = thresholds[[threshold]](across_expr)) |>
    dplyr::ungroup() |>
    dplyr::filter(detected) |>
    dplyr::count(cluster, threshold)

  detection_df$threshold <- factor(detection_df$threshold, levels = names(thresholds))

  p <- ggplot2::ggplot(detection_df, ggplot2::aes(x = cluster, y = n, fill = threshold)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::labs(title = "Gene Detection per Cluster (CP10K thresholds)",
                  x = "Cluster", y = "Number of genes", fill = "Detection Threshold") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    )

  if (save_plot) {
    if (!dir.exists("plots")) dir.create("plots")
    ggplot2::ggsave("plots/CP10K_detection_stats.pdf", p, width = 10, height = 6)
    message("Saved plot to plots/CP10K_detection_stats.pdf")
  }

  return(p)
}