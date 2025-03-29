#' Plot Gene Expression per Cluster across Samples using Patchwork
#'
#' Creates boxplots of gene expression per cluster across samples,
#' using patchwork to display all plots in a single layout.
#'
#' @param expr_list Output from AverageExpressionPerClusterPerSample().
#' @param genes Character vector of gene names to include. If NULL, uses all.
#' @param base_width Numeric, base width per sample (inches). Default: 1.5
#' @param base_height Numeric, base height per cluster (inches). Default: 3
#'
#' @return A patchwork object of all cluster plots.
#' @export
PlotBulkExpressionViolins <- function(expr_list,
                                         genes = NULL,
                                         base_width = 1.5,
                                         base_height = 3) {
  if (length(expr_list) == 0) stop("Empty input list")

  # Combine matrices into long format
  long_df <- purrr::map_dfr(names(expr_list), function(sample_name) {
    mat <- expr_list[[sample_name]]
    if (is.null(mat)) return(NULL)
    mat <- as.data.frame(mat)
    mat$gene <- rownames(mat)
    long <- tidyr::pivot_longer(mat, -gene, names_to = "cluster", values_to = "expression")
    long$sample <- sample_name
    return(long)
  })

  if (!is.null(genes)) {
    long_df <- long_df %>% filter(gene %in% genes)
  }

  # Clean sample names
  long_df$sample <- gsub("^sample_|_obj[0-9]+$", "", long_df$sample)

  # Plot per cluster
  split_data <- split(long_df, long_df$cluster)

  plots <- lapply(names(split_data), function(cluster) {
    df <- split_data[[cluster]]
    ggplot(df, aes(x = sample, y = expression, fill = gene)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.8, position = position_dodge(width = 0.75)) +
      labs(title = paste("Cluster", cluster), x = "Sample", y = "Expression") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })

  n_clusters <- length(plots)
  n_samples <- length(unique(long_df$sample))

  # Dynamic layout
  ncol <- ceiling(sqrt(n_clusters))
  nrow <- ceiling(n_clusters / ncol)

  plot_width <- max(base_width * n_samples, 5)
  plot_height <- base_height * nrow

  combined <- wrap_plots(plots, ncol = ncol) +
    plot_annotation(title = "Expression per Gene per Cluster")

  attr(combined, "width") <- plot_width
  attr(combined, "height") <- plot_height

  return(combined)
}
