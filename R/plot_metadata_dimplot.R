#' Plot Seurat Reductions Colored by Metadata (DimPlot or FeaturePlot)
#'
#' Uses Seurat::DimPlot (for categorical metadata) or FeaturePlot (for continuous)
#' to visualize projections, optionally saving plots to PDF with 1:1 aspect ratio.
#'
#' @param seurat_obj A Seurat object.
#' @param meta_col Metadata column to color by.
#' @param reduction Character vector of reductions to use (e.g., "umap", "tsne").
#' @param pt.size Point size for cells.
#' @param verbose Logical; print progress messages.
#' @param save_as_pdf Logical; whether to save plots to PDF.
#' @param pdf_dir Directory where PDFs should be saved.
#'
#' @return A ggplot object or list of ggplot objects.
plot_metadata_dimplot <- function(seurat_obj,
                                   meta_col,
                                   reduction = NULL,
                                   pt.size = 1,
                                   verbose = TRUE,
                                   save_as_pdf = TRUE,
                                   pdf_dir = "plots") {
  if (!meta_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste0("❌ Metadata column '", meta_col, "' not found in Seurat object."))
  }

  if (is.null(reduction)) {
    reduction <- Seurat::Reductions(seurat_obj)
    if (verbose) message("📌 Using all reductions: ", paste(reduction, collapse = ", "))
  }

  missing <- setdiff(reduction, Seurat::Reductions(seurat_obj))
  if (length(missing) > 0) stop("❌ Reduction(s) not found: ", paste(missing, collapse = ", "))

  is_continuous <- is.numeric(seurat_obj@meta.data[[meta_col]])
  plot_list <- list()

  if (save_as_pdf && !dir.exists(pdf_dir)) {
    dir.create(pdf_dir, recursive = TRUE)
    if (verbose) message("📂 Created output directory: ", pdf_dir)
  }

  for (red in reduction) {
    if (verbose) message("🔍 Plotting '", meta_col, "' using ", if (is_continuous) "FeaturePlot" else "DimPlot", " on reduction: '", red, "'")

    if (is_continuous) {
      p <- Seurat::FeaturePlot(
        seurat_obj,
        features = meta_col,
        reduction = red,
        pt.size = pt.size
      )
    } else {
      Idents(seurat_obj) <- meta_col
      p <- Seurat::DimPlot(
        seurat_obj,
        reduction = red,
        group.by = meta_col,
        pt.size = pt.size
      )
    }

    # Add 1:1 aspect ratio and title
    p <- p +
      ggplot2::labs(title = paste0(red, " - ", meta_col)) +
      ggplot2::theme(aspect.ratio = 1)

    plot_list[[red]] <- p
  }

  if (save_as_pdf && length(plot_list) > 0) {
    pdf_file <- file.path(pdf_dir, paste0("DimPlot_", meta_col, ".pdf"))
    n_per_page <- 4

    grDevices::pdf(pdf_file, width = 10, height = 10)
    for (i in seq(1, length(plot_list), by = n_per_page)) {
      patch <- patchwork::wrap_plots(plot_list[i:min(i + n_per_page - 1, length(plot_list))])
      print(patch)
    }
    grDevices::dev.off()
    message("💾 Saved PDF to: ", pdf_file)
  }

  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  } else {
    return(plot_list)
  }
}
