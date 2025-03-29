#' Plot Top Marker Genes Across Clusters
#'
#' Identifies top 6 marker genes (based on avg_log2FC) across clusters
#' and either returns FeaturePlot objects or saves a 2x3 panel PDF.
#'
#' @param seurat_obj Processed Seurat object with clusters and UMAP.
#' @param cluster_col Metadata column with cluster labels (e.g., "clusters_sct").
#' @param use.assay Assay to use for marker detection (e.g. "SCT").
#' @param output_dir Directory to save the plot (if return_plots = FALSE).
#' @param sample_name Sample name for filename.
#' @param umap_name Reduction name for UMAP (default: "umap_sct").
#' @param return_plots If TRUE, returns the FeaturePlot objects instead of saving (default: FALSE).
#'
#' @return List of ggplot objects if return_plots = TRUE, otherwise NULL (saves PDF).
plot_top_marker_genes <- function(seurat_obj,
                                  cluster_col,
                                  use.assay = "SCT",
                                  output_dir = "plots",
                                  sample_name = "sample",
                                  umap_name = "umap_sct",
                                  return_plots = FALSE) {
  require(Seurat)
  require(cowplot)
  require(viridis)
  require(dplyr)

  message("🎯 Finding top marker genes...")
  Idents(seurat_obj) <- seurat_obj[[cluster_col, drop = TRUE]]

  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, assay = use.assay)
  top6_genes <- markers %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = avg_log2FC) %>%
    ungroup() %>%
    top_n(n = 6, wt = avg_log2FC) %>%
    arrange(desc(avg_log2FC)) %>%
    pull(gene) %>%
    unique()

  message("🖼 Plotting top 6 marker genes: ", paste(top6_genes, collapse = ", "))

  feature_plots <- lapply(top6_genes, function(gene) {
    FeaturePlot(seurat_obj, features = gene, reduction = umap_name) +
      ggtitle(paste("UMAP:", gene, "-", sample_name)) +
      scale_color_viridis_c()
  })

  if (return_plots) {
    return(feature_plots)
  } else {
    sample_dir <- file.path(output_dir, sample_name)
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    marker_pdf_path <- file.path(sample_dir, paste0("Top6_MarkerGenes_", sample_name, ".pdf"))

    pdf(marker_pdf_path, width = 11, height = 8.5)
    print(plot_grid(plotlist = feature_plots, ncol = 3))
    dev.off()

    message("✅ Saved marker gene panel to: ", marker_pdf_path)
    return(invisible(NULL))
  }
}
