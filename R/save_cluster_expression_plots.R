save_cluster_expression_plots <- function(split_data,
                                          log_scale = FALSE,
                                          file_prefix = "cluster_expression",
                                          save_as_pdf = TRUE,
                                          threads = 1,
                                          point = 0.1,
                                          y_max = NULL,
                                          plot_dir = "plots",
                                          combine_clusters = FALSE,
                                          fill_by_clusters = FALSE) {
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

  # Merge all data if needed
  if (combine_clusters || fill_by_clusters) {
    long_df <- dplyr::bind_rows(split_data, .id = "cluster")

    p <- ggplot2::ggplot(long_df, ggplot2::aes(x = sample, y = expression)) +
      ggbeeswarm::geom_quasirandom(
        ggplot2::aes_string(fill = if (fill_by_clusters) "cluster" else "assay"),
        dodge.width = 0.9,
        color = "black",
        alpha = 0.5,
        size = 0.01
      ) +
      ggplot2::geom_boxplot(
        ggplot2::aes_string(fill = if (fill_by_clusters) "cluster" else "assay"),
        outlier.shape = NA,
        alpha = 0.6,
        position = ggplot2::position_dodge(width = 0.9),
        linewidth = 0.1
      ) +
      ggplot2::labs(
        title = "Combined Cluster Expression",
        x = "Sample",
        y = if (log_scale) "Log10(Expression)" else "Expression"
      ) +
      ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(size = 12, face = "bold"),
        legend.position = "bottom"
      )

    if (log_scale) {
      p <- p + ggplot2::scale_y_log10()
    } else if (!is.null(y_max)) {
      p <- p + ggplot2::coord_cartesian(ylim = c(0, y_max))
    }

    # Save single plot
    if (save_as_pdf) {
      if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
      out_file <- file.path(
        plot_dir,
        paste0(file_prefix, "_combined", if (log_scale) "_log10" else "_linear", ".pdf")
      )
      message(sprintf("Saving combined plot to: %s", out_file))
      rasterpdf::raster_pdf(file = out_file, width = 10, height = 6, res = 300)
      print(p)
      rasterpdf::dev.off()
    }

    return(invisible(list(combined = p)))
  }

  # Per-cluster plotting mode
  cluster_plotter <- if (threads > 1 && requireNamespace("future.apply", quietly = TRUE)) {
    message(sprintf("Using future.apply with %d threads...", threads))
    future::plan(future::multisession, workers = threads)
    on.exit(future::plan(future::sequential), add = TRUE)
    future.apply::future_lapply
  } else if (requireNamespace("pbapply", quietly = TRUE)) {
    message("Using pbapply with progress bar...")
    pbapply::pblapply
  } else {
    message("Using base lapply (no progress or parallelism)...")
    lapply
  }

  message("Generating plots for each cluster...")

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
        size = point,
      ) +
      ggplot2::labs(
        title = paste("Cluster", cluster_name),
        x = "Sample",
        y = if (log_scale) "Log10(Expression)" else "Expression"
      ) +
      ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(size = 12, face = "bold"),
        legend.position = "bottom"
      )

    if (log_scale) {
      p <- p + ggplot2::scale_y_log10()
    } else if (!is.null(y_max)) {
      p <- p + ggplot2::coord_cartesian(ylim = c(0, y_max))
    }

    return(p)
  })

  names(plotlist) <- names(split_data)

  if (!save_as_pdf) return(plotlist)

  # Save multiple plots to PDF
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  out_file <- file.path(
    plot_dir,
    paste0(file_prefix, if (log_scale) "_log10" else "_linear", ".pdf")
  )

width_factor <- 1
if (fill_by_clusters) {
  n_fill <- length(unique(long_df$cluster))
  width_factor <- n_fill / 3
}
rasterpdf::raster_pdf(file = out_file, width = 10 * width_factor, height = 6, res = 300)

  for (cluster_name in names(plotlist)) {
    message(sprintf("  → Printing plot for cluster '%s'...", cluster_name))
    print(plotlist[[cluster_name]])
  }
  rasterpdf::dev.off()

  message("✅ Finished saving all plots.")
  return(invisible(plotlist))
}
