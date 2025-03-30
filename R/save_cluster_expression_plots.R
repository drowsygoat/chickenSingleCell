#' Save Cluster Expression Plots to Rasterized PDF or Return Plot List
#'
#' This function creates one expression plot per cluster from a list of data frames,
#' and either saves the plots to a rasterized multi-page PDF or returns them as a list
#' of ggplot objects. Each plot shows expression values across samples, with boxplots
#' and overlaid jittered points. Users can choose between linear and log10-transformed
#' y-axes.
#'
#' @param split_data A named list of data frames, one per cluster. Each data frame
#'   must contain the columns: `sample`, `expression`, and `assay`.
#' @param log_scale Logical. If TRUE, a log10-transformed y-axis is used for all plots.
#'   Default: FALSE (linear y-axis).
#' @param file_prefix Character. Prefix for the output PDF filename. The plot will be saved as
#'   `plots/<file_prefix>_linear.pdf` or `plots/<file_prefix>_log10.pdf` depending on the y-axis scale.
#' @param plot Logical. If TRUE (default), saves the plots to PDF. If FALSE, returns a list of ggplot objects instead.
#'
#' @return Either the filename (if `plot = TRUE`) or a named list of ggplot objects (if `plot = FALSE`).
#' @export
save_cluster_expression_plots <- function(split_data, log_scale = FALSE, file_prefix = "cluster_expression", plot = TRUE) {
  # Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("The 'ggplot2' package is required but not installed.")
  if (!requireNamespace("ggbeeswarm", quietly = TRUE)) stop("The 'ggbeeswarm' package is required but not installed.")
  if (plot && !requireNamespace("ragg", quietly = TRUE)) stop("The 'ragg' package is required for rasterized PDF output.")

  # Use pbapply if available, fallback to lapply
  if (requireNamespace("pbapply", quietly = TRUE)) {
    message("Generating plots for each cluster...")
    cluster_plotter <- pbapply::pblapply
  } else {
    message("Package 'pbapply' not installed. Using standard lapply...")
    cluster_plotter <- lapply
  }

  # Generate a ggplot for each cluster
  plotlist <- cluster_plotter(names(split_data), function(cluster_name) {
    cluster_df <- split_data[[cluster_name]]

    p <- ggplot2::ggplot(cluster_df, ggplot2::aes(x = sample, y = expression)) +
      ggplot2::geom_boxplot(
        ggplot2::aes(fill = assay),
        outlier.shape = NA,
        alpha = 0.6,
        position = ggplot2::position_dodge(width = 0.9),
        linewidth = 0.1
      ) +
      ggbeeswarm::geom_quasirandom(
        dodge.width = 0.9,
        color = "black",
        alpha = 0.5,
        size = 0.01
      ) +
      ggplot2::labs(
        title = paste("Cluster", cluster_name),
        x = "Sample",
        y = if (log_scale) "Log10(Expression)" else "Expression"
      ) +
      ggplot2::theme_minimal(base_size = 6) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(size = 12, face = "bold"),
        legend.position = "bottom"
      )

    if (log_scale) {
      p <- p + ggplot2::scale_y_log10()
    }

    return(p)
  })

  # Assign names to plot list
  names(plotlist) <- names(split_data)

  # Return plot list if not saving
  if (!plot) {
    message("Returning plot list (plot = FALSE).")
    return(plotlist)
  }

  # Ensure output directory exists
  if (!dir.exists("plots")) {
    message("Creating output directory 'plots/'...")
    dir.create("plots")
  }

  # Construct filename
  file_name <- file.path(
    "plots",
    paste0(file_prefix, if (log_scale) "_log10" else "_linear", ".pdf")
  )

  # Save to rasterized PDF
  message(sprintf("Saving plots to rasterized PDF: %s", file_name))
  rasterpdf::raster_pdf(file = file_name, width = 8, height = 6, units = "in", res = 300)

  for (cluster_name in names(plotlist)) {
    message(sprintf("  → Printing plot for cluster '%s'...", cluster_name))
    print(plotlist[[cluster_name]])
  }

  dev.off()
  message("Finished saving all plots.")
  message(sprintf("Saved %s-scale expression plots to: %s", if (log_scale) "log10" else "linear", file_name))

  return(invisible(file_name))
}
