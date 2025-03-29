#' Generate Violin Plots for Numeric Metadata Columns Specific to an Assay
#'
#' For all numeric metadata columns ending in "_<assay>", return violin plots
#' grouped by cluster or save as a cumulative PDF.
#'
#' @param seurat_obj A Seurat object with cluster info.
#' @param assay Assay name used for suffix matching (e.g. "SCT").
#' @param cluster_col Column name for grouping (e.g. "clusters_sct").
#' @param sample_name Sample name (used for plot titles and filenames).
#' @param output_dir Directory to save plots (used if return_plots = FALSE).
#' @param return_plots If TRUE (default), return the plots; otherwise save as a PDF and return NULL.
#'
#' @return Named list of ggplot objects if return_plots = TRUE, otherwise NULL.
plot_vln_metadata_by_assay <- function(seurat_obj,
                                       assay = "SCT",
                                       cluster_col = "clusters_sct",
                                       sample_name = "sample",
                                       output_dir = "plots",
                                       return_plots = TRUE) {
  require(Seurat)
  require(ggplot2)
  require(dplyr)

  suffix_pattern <- paste0("_", assay, "$")
  meta_matches <- grep(suffix_pattern, colnames(seurat_obj@meta.data), value = TRUE)

  if (length(meta_matches) == 0) {
    warning("No metadata columns matched assay suffix pattern: ", suffix_pattern)
    return(if (return_plots) list() else invisible(NULL))
  }

  # Filter to numeric columns only
  numeric_meta <- meta_matches[
    sapply(meta_matches, function(col) is.numeric(seurat_obj@meta.data[[col]]))
  ]

  if (length(numeric_meta) == 0) {
    warning("No numeric metadata columns matched for assay: ", assay)
    return(if (return_plots) list() else invisible(NULL))
  }

  message("🎻 Preparing violin plots for numeric metadata: ", paste(numeric_meta, collapse = ", "))
  plot_list <- list()

  for (var in numeric_meta) {
    base <- sub(suffix_pattern, "", var)
    p <- VlnPlot(seurat_obj, features = var, group.by = cluster_col, pt.size = 0.1) +
      ggtitle(paste("Violin Plot:", base, "-", sample_name))
    plot_list[[paste0("vln_", base)]] <- p
  }

  if (return_plots) {
    return(plot_list)
  } else {
    sample_dir <- file.path(output_dir, sample_name)
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    vln_pdf_path <- file.path(sample_dir, paste0("Cumulative_ViolinPlots_", sample_name, "_", assay, ".pdf"))

    pdf(vln_pdf_path, width = 11, height = 8.5)
    for (i in seq(1, length(plot_list), by = 4)) {
      print(cowplot::plot_grid(plotlist = plot_list[i:min(i + 3, length(plot_list))], ncol = 2))
    }
    dev.off()
    message("✅ Saved cumulative violin plot PDF to: ", vln_pdf_path)
    return(invisible(NULL))
  }
}
