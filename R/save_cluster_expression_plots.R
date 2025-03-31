#' Save Cluster Expression Plots to Rasterized PDF or Return Plot List
#'
#' This function creates one expression plot per cluster from a list of data frames
#' and either saves the plots to a rasterized multi-page PDF or returns them as a list
#' of ggplot objects. Each plot shows expression values across samples, using
#' boxplots with quasirandom jitter. The function supports log10 transformation of the
#' y-axis and optional multithreading for faster plot generation.
#'
#' @param split_data A named list of data frames, one per cluster. Each data frame
#'   must contain the columns: \code{sample}, \code{expression}, and \code{assay}.
#' @param log_scale Logical. If \code{TRUE}, a log10-transformed y-axis is used.
#'   Default: \code{FALSE} (linear y-axis).
#' @param file_prefix Character. Prefix for the output PDF filename. The plot will be saved as
#'   \code{plots/<file_prefix>_linear.pdf} or \code{plots/<file_prefix>_log10.pdf} depending on the y-axis scale.
#' @param plot Logical. If \code{TRUE} (default), saves the plots to a rasterized PDF.
#'   If \code{FALSE}, returns the list of ggplot objects without saving.
#' @param threads Integer. Number of threads to use for generating plots. If greater than 1,
#'   enables parallel processing via \code{future.apply}. Default: 1 (no parallelism).
#' #' @param y_max Optional numeric. Maximum y-axis limit for linear plots to trim outliers.
#'   Only used when \code{log_scale = FALSE}. Default: \code{NULL} (no limit).
#'
#' @return A named list of ggplot objects (one per cluster). If \code{plot = TRUE},
#'   the list is returned invisibly after saving to PDF.
#'
#' @details
#' The function uses either \code{pbapply::pblapply()} or \code{future.apply::future_lapply()}
#' depending on the value of \code{threads}. Plot rendering to PDF is done with
#' \code{ragg::raster_pdf()} for high-performance rasterized output.
#'
#' If the output directory \code{plots/} does not exist, it is created automatically.
#'
#' @examples
#' # Generate and return plots without saving
#' plots <- save_cluster_expression_plots(split_data, plot = FALSE)
#'
#' # Save rasterized plots and capture the list invisibly
#' plots <- save_cluster_expression_plots(split_data, plot = TRUE, threads = 4)
#'
#' @export
save_cluster_expression_plots <- function(split_data,
                                          log_scale = FALSE,
                                          file_prefix = "cluster_expression",
                                          save_as_pdf = TRUE,
                                          threads = 1,
                                          y_max = NULL) {
  # Required packages
  required_pkgs <- c("ggplot2", "ggbeeswarm", if (save_as_pdf) "rasterpdf")

  invisible(lapply(required_pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("The '%s' package is required but not installed.", pkg), call. = FALSE)
    }
  }))

  # Check for empty input
  if (length(split_data) == 0 || all(lengths(split_data) == 0)) {
    warning("No data provided: 'split_data' is empty. Returning NULL.")
    return(NULL)
  }

  # Select plot generator with progress and threading support
  cluster_plotter <- NULL

  if (threads > 1) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("You must install the 'future.apply' package to use multithreading (threads > 1).")
    }
    message(sprintf("Using future.apply with %d threads...", threads))
    future::plan(future::multisession, workers = threads)
    on.exit(future::plan(future::sequential), add = TRUE)  # Reset plan after function exits
    cluster_plotter <- future.apply::future_lapply
  } else if (requireNamespace("pbapply", quietly = TRUE)) {
    message("Using pbapply with progress bar...")
    cluster_plotter <- pbapply::pblapply
  } else {
    message("Using base lapply (no progress or parallelism)...")
    cluster_plotter <- lapply
  }

  # Static plot components
  boxplot_layer <- ggplot2::geom_boxplot(
    ggplot2::aes(fill = assay),
    outlier.shape = NA,
    alpha = 0.6,
    position = ggplot2::position_dodge(width = 0.9),
    linewidth = 0.1
  )

  jitter_layer <- ggbeeswarm::geom_quasirandom(
    dodge.width = 0.9,
    color = "black",
    alpha = 0.5,
    size = 0.01
  )

  theme_base <- ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.position = "bottom"
    )

  # Generate plots
  message("Generating plots for each cluster...")

  plotlist <- cluster_plotter(names(split_data), function(cluster_name) {
    cluster_df <- split_data[[cluster_name]]

    p <- ggplot2::ggplot(cluster_df, ggplot2::aes(x = sample, y = expression)) +
      boxplot_layer +
      jitter_layer +
      ggplot2::labs(
        title = paste("Cluster", cluster_name),
        x = "Sample",
        y = if (log_scale) "Log10(Expression)" else "Expression"
      ) +
      theme_base

    if (log_scale) {
      p <- p + ggplot2::scale_y_log10()
    } else if (!is.null(y_max)) {
      p <- p + ggplot2::coord_cartesian(ylim = c(0, y_max))
    }

    return(p)
  })

  names(plotlist) <- names(split_data)

  # If not saving, return plotlist directly
  if (!save_as_pdf) {
    return(plotlist)
  }

  # Ensure output directory exists
  if (!dir.exists("plots")) {
    if (requireNamespace("fs", quietly = TRUE)) {
      fs::dir_create("plots")
    } else {
      dir.create("plots", recursive = TRUE)
    }
    message("Created 'plots/' directory.")
  }

  # Define output file
  file_name <- file.path(
    "plots",
    paste0(file_prefix, if (log_scale) "_log10" else "_linear", ".pdf")
  )

  message(sprintf("Saving plots to rasterized PDF: %s", file_name))

  rasterpdf::raster_pdf(file = file_name, width = 8, height = 6, units = "in", res = 300)

    # pdf(file = file_name, width = 8, height = 6)

  for (cluster_name in names(plotlist)) {
    message(sprintf("  → Printing plot for cluster '%s'...", cluster_name))
    print(plotlist[[cluster_name]])
  }

  rasterpdf::dev.off()

  message("✅ Finished saving all plots.")
  message(sprintf("✅ Saved %s-scale expression plots to: %s", if (log_scale) "log10" else "linear", file_name))

  # Still return plotlist invisibly
  return(invisible(plotlist))
}
