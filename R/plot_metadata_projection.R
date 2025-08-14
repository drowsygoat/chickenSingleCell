#' Plot Seurat Dimensional Reductions Colored by Metadata
#'
#' Plots dimensional reductions (e.g., UMAP, t-SNE) from a Seurat object,
#' coloring cells by a selected metadata column. Supports saving to PDF.
#'
#' @param seurat_obj A Seurat object.
#' @param meta_col Metadata column to color by.
#' @param reduction One or more reductions (e.g., "umap", "tsne"). If NULL, use all.
#' @param pt.size Point size.
#' @param alpha Point transparency.
#' @param verbose Logical; print progress messages.
#' @param save_as_pdf Logical; whether to save the plot(s) to PDF.
#' @param pdf_dir Directory to save PDF (default = "plots").
#'
#' @return A ggplot object (single reduction) or named list of ggplots (multiple reductions).
#' 
#' 
# this one uses GGPLOT!!!!
plot_metadata_projection <- function(seurat_obj,
                                     meta_col,
                                     reduction = NULL,
                                     pt.size = 1,
                                     alpha = 0.9,
                                     verbose = TRUE,
                                     save_as_pdf = TRUE,
                                     pdf_dir = "plots") {
  # Check metadata column
  if (!meta_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste0("❌ Metadata column '", meta_col, "' not found in Seurat object."))
  }

  # Use all reductions if not specified
  available_reductions <- Seurat::Reductions(seurat_obj)
  if (is.null(reduction)) {
    reduction <- available_reductions
    if (verbose) message("📌 No reduction specified. Using all: ", paste(reduction, collapse = ", "))
  } else {
    missing <- setdiff(reduction, available_reductions)
    if (length(missing) > 0) {
      stop("❌ These reductions are missing: ", paste(missing, collapse = ", "))
    }
  }

  # Create output directory if saving
  if (save_as_pdf && !dir.exists(pdf_dir)) {
    dir.create(pdf_dir, recursive = TRUE)
    if (verbose) message("📂 Created PDF output directory: ", pdf_dir)
  }

  # Collect plots
  plot_list <- list()

  for (red in reduction) {
    if (verbose) message("🔍 Plotting '", meta_col, "' on reduction: '", red, "'")

    plot_data <- Seurat::FetchData(seurat_obj, vars = c(meta_col, paste0(red, "_1"), paste0(red, "_2")))
    colnames(plot_data) <- c("meta", "x", "y")
    is_continuous <- is.numeric(plot_data$meta)

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(ggplot2::aes(color = meta), size = pt.size, alpha = alpha) +
      ggplot2::labs(
        title = paste0(red, " - ", meta_col),
        x = paste0(toupper(red), "_1"),
        y = paste0(toupper(red), "_2"),
        color = meta_col
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(aspect.ratio = 1)

    p <- if (is_continuous) {
      p + ggplot2::scale_color_viridis_c()
    } else {
      p + ggplot2::scale_color_manual(
        values = Seurat::DiscretePalette(length(unique(plot_data$meta)), "polychrome")
      )
    }

    plot_list[[red]] <- p
  }

  # Save to PDF
  if (save_as_pdf) {
    pdf_file <- file.path(pdf_dir, paste0("projection_", meta_col, ".pdf"))
    n_per_page <- 4
    if (verbose) message("💾 Saving plot(s) to PDF: ", pdf_file)

    grDevices::pdf(pdf_file, width = 10, height = 10)
    for (i in seq(1, length(plot_list), by = n_per_page)) {
      patch <- patchwork::wrap_plots(plot_list[i:min(i + n_per_page - 1, length(plot_list))])
      print(patch)
    }
    grDevices::dev.off()
  }

  # Return either a single plot or a named list
  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  } else {
    return(plot_list)
  }
}
