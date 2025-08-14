#' Cluster and Explore an Integrated Seurat Object
#'
#' Performs dimensionality reduction, clustering, and visualizations on an integrated Seurat object.
#' Generates cumulative cluster plots as well as individual per-sample highlight plots.
#'
#' @param seurat_obj Integrated Seurat object.
#' @param dims PCA dimensions to use (default: 1:30).
#' @param resolution Clustering resolution (default: 1).
#' @param use.assay Assay to use (default: "integrated").
#' @param tsne_perplexity Perplexity for t-SNE (default: 30).
#' @param output_dir Output directory for plots (default: "cluster_and_explore_integrated_plots").
#' @param save_metadata_csv Logical; whether to save metadata (default: FALSE).
#'
#' @return Updated Seurat object with dimensional reductions and clusters.
#' @export
cluster_and_explore_integrated <- function(seurat_obj,
                                           dims = 1:30,
                                           resolution = 1,
                                           use.assay = "integrated",
                                           tsne_perplexity = 30,
                                           output_dir = "cluster_and_explore_integrated_plots",
                                           save_metadata_csv = FALSE) {
  require(Seurat)
  require(ggplot2)
  require(patchwork)

  DefaultAssay(seurat_obj) <- use.assay
  sample_name <- seurat_obj@project.name %||% "merged"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  message("📦 Processing integrated object: ", sample_name)

  # Dimensionality reduction
  seurat_obj <- RunPCA(seurat_obj, npcs = max(dims), verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = dims, reduction.name = "umap", reduction.key = "UMAP_")
  seurat_obj <- RunTSNE(seurat_obj, dims = dims, reduction.name = "tsne", reduction.key = "TSNE_")

  # Clustering
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  seurat_obj$clusters <- seurat_obj$seurat_clusters

  # Save metadata
  if (save_metadata_csv) {
    metadata_file <- file.path(output_dir, paste0("metadata_", sample_name, ".csv"))
    write.csv(seurat_obj@meta.data, metadata_file)
    message("📄 Metadata saved to: ", metadata_file)
  }

  # --- Cumulative plots ---
  p_umap_all <- DimPlot(seurat_obj, reduction = "umap", group.by = "clusters", label = TRUE) +
    ggtitle("UMAP: All Clusters")

  p_tsne_all <- DimPlot(seurat_obj, reduction = "tsne", group.by = "clusters", label = TRUE) +
    ggtitle("t-SNE: All Clusters")

  ggsave(file.path(output_dir, "UMAP_All_Clusters.pdf"), p_umap_all, width = 7, height = 5)
  
  ggsave(file.path(output_dir, "tSNE_All_Clusters.pdf"), p_tsne_all, width = 7, height = 5)

  # --- Per-sample highlight plots (one page) ---
  sample_col <- if ("sample" %in% colnames(seurat_obj@meta.data)) "sample" else "orig.ident"

  samples <- unique(seurat_obj[[sample_col]][, 1])


  umap_list <- lapply(samples, function(s) {
    DimPlot(subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj[[sample_col]][,1] == s]),

    reduction = "umap", group.by = "clusters") +
    ggtitle(paste("UMAP:", s))
  })

  tsne_list <- lapply(samples, function(s) {
    DimPlot(subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj[[sample_col]][,1] == s]),
          reduction = "tsne", group.by = "clusters") +
    ggtitle(paste("t-SNE:", s))
  })


  ggsave(file.path(output_dir, "UMAP_Per_Sample.pdf"), wrap_plots(umap_list), width = 12, height = 2.5 * ceiling(length(umap_list)/2))

  ggsave(file.path(output_dir, "tSNE_Per_Sample.pdf"), wrap_plots(tsne_list), width = 12, height = 2.5 * ceiling(length(tsne_list)/2))

  message("✅ Dimensionality reduction and clustering complete.")
  return(seurat_obj)
}


  # tsne_list <- lapply(samples, function(s) {
  #   DimPlot(seurat_obj, reduction = "tsne", group.by = "clusters",
  #           cells.highlight = colnames(seurat_obj)[seurat_obj[[sample_col]][,1] == s]) +
  #     ggtitle(paste("t-SNE:", s))


  # umap_list <- lapply(samples, function(s) {
  #   DimPlot(seurat_obj, reduction = "umap", group.by = "clusters",
  #           cells.highlight = colnames(seurat_obj)[seurat_obj[[sample_col]][,1] == s]) +
  #     ggtitle(paste("UMAP:", s))
  # })