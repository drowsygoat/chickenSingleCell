#' Plot Combined Gene Expression per Cluster across Samples
#'
#' Creates one big boxplot (with quasirandom jitter) of gene expression across samples,
#' with clusters optionally shown as fill or facet.
#'
#' @param expr_lists Named list of expr_list objects (from AverageExpressionPerClusterPerSample),
#'                   each representing a normalization method.
#' @param genes Character vector of gene names to include. If NULL, uses all.
#' @param save_as_pdf Logical, if TRUE saves the plot as PDF. Default: TRUE
#' @param subset Numeric between 0 and 1, to randomly downsample data (optional).
#' @param y_max Optional numeric to set the y-axis maximum.
#' @param plot_dir Directory to save plots. Default: "plots"
#' @param linear_scale Logical, if TRUE uses linear y-axis instead of log10. Default: FALSE
#' @param facet_clusters Logical, whether to facet by cluster. Default: FALSE
#' @param color_by_cluster Logical, whether to use cluster as fill. Default: FALSE
#' @param combine_clusters Logical, required TRUE for this function to run. Default: TRUE
#' @param file_prefix Optional file name prefix for saving plots.
#'
#' @return A single ggplot object (or NULL if no data)
#' @export
save_combined_cluster_expression_plot <- function(expr_lists,
                                                  genes = NULL,
                                                  save_as_pdf = TRUE,
                                                  subset = NULL,
                                                  y_max = NULL,
                                                  plot_dir = "plots",
                                                  log_scale = TRUE,
                                                  facet_clusters = FALSE,
                                                  color_by_cluster = FALSE,
                                                  combine_clusters = TRUE,
                                                  file_prefix = "combined_expression") {
  required_pkgs <- c("ggplot2", "tidyr", "dplyr", "purrr", "ggbeeswarm", "rasterpdf")
  invisible(lapply(required_pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
    }
  }))

  if (!combine_clusters) {
    stop("This version of save_combined_cluster_expression_plot only supports combine_clusters = TRUE")
  }

  # Combine all long-format data
  long_df <- purrr::imap_dfr(expr_lists, function(expr_list, method_name) {
    purrr::map_dfr(names(expr_list), function(sample_name) {
      mat <- expr_list[[sample_name]]
      if (is.null(mat)) return(NULL)
      if (!is.null(genes)) mat <- mat[rownames(mat) %in% genes, , drop = FALSE]

      mat <- as.data.frame(mat)
      mat$gene <- rownames(mat)

      tidyr::pivot_longer(mat, cols = -gene, names_to = "cluster", values_to = "expression") %>%
        dplyr::filter(expression > 0) %>%
        dplyr::mutate(
          sample = gsub("_.*$", "", sample_name),
          method = method_name
        )
    })
  })

  if (nrow(long_df) == 0) {
    warning("No expression data available after filtering.")
    return(NULL)
  }

  # Optional downsampling
  if (!is.null(subset)) {
    if (!is.numeric(subset) || subset <= 0 || subset > 1) {
      stop("subset must be a numeric value between 0 and 1.")
    }
    set.seed(42)
    long_df <- dplyr::slice_sample(long_df, prop = subset)
  }

  # Plot
  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = sample, y = expression)) +
    ggbeeswarm::geom_quasirandom(
      ggplot2::aes_string(fill = if (color_by_cluster) "cluster" else "method"),
      color = "black", alpha = 0.5, groupOnX = TRUE
    ) +
    ggplot2::geom_boxplot(
      ggplot2::aes_string(fill = if (color_by_cluster) "cluster" else "method"),
      outlier.shape = NA, alpha = 0.3, width = 0.5,
      color = "black", position = ggplot2::position_dodge(width = 0.7)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Combined Cluster Expression", y = "Expression", x = "Sample") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (facet_clusters) {
    p <- p + ggplot2::facet_wrap(~cluster)
  }

  if (log_scale) {
    if (!is.null(y_max)) {
      p <- p + ggplot2::coord_cartesian(ylim = c(0, y_max))
    }
  } else {
    p <- p + ggplot2::scale_y_log10()
    if (!is.null(y_max)) {
      p <- p + ggplot2::coord_cartesian(ylim = c(NA, y_max))
    }
  }

  # Save
  if (save_as_pdf) {
    if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
    rasterpdf::raster_pdf(file = file.path(plot_dir, paste0(file_prefix, ".pdf")),
                          width = if (facet_clusters) 14 else 10,
                          height = if (facet_clusters) 10 else 6,
                          res = 300)
    print(p)
    grDevices::dev.off()
  }

  return(p)
}
