#' Perform Clustering and Dimensionality Reduction with Optional Subsetting
#'
#' This function runs PCA, UMAP, t-SNE, and clustering on a Seurat object.
#' It optionally downsamples the object to a fraction of cells for faster runtime.
#' Individual plots are saved as PDFs and RDS files in a sample-specific subfolder.
#'
#' @param seurat_obj A Seurat object with RNA or SCT assay.
#' @param dims Numeric vector of dimensions (PCs) to use (default: 1:30).
#' @param resolution Clustering resolution (default: 0.5).
#' @param use.assay Name of assay to set as default (default: "SCT").
#' @param tsne_perplexity Perplexity for t-SNE (default: 30).
#' @param output_dir Directory to save plots (default: "plots").
#' @param sample_name Optional name for output; defaults to project name.
#' @param subset_fraction Fraction of cells to randomly retain (0–1); default = 1 (all cells).
#' @param subset_seed Random seed for reproducibility (default: 42).
#'
#' @return A Seurat object with dimensionality reductions and clustering metadata.
#'
#' @export
cluster_and_explore <- function(seurat_obj,
                                dims = 1:30,
                                resolution = 0.5,
                                use.assay = "SCT",
                                tsne_perplexity = 30,
                                output_dir = "plots",
                                sample_name = NULL,
                                subset_fraction = 1.0,
                                subset_seed = 42) {
  DefaultAssay(seurat_obj) <- use.assay

  if (is.null(sample_name)) {
    sample_name <- seurat_obj@project.name
  }

  # Subset randomly if desired
  if (subset_fraction < 1) {
    set.seed(subset_seed)
    total_cells <- ncol(seurat_obj)
    n_subset <- floor(total_cells * subset_fraction)
    selected_cells <- sample(colnames(seurat_obj), n_subset)
    seurat_obj <- subset(seurat_obj, cells = selected_cells)
    message("Subset to ", n_subset, " cells (", round(subset_fraction * 100, 1), "%)")
  }

  # Output folder per sample
  sample_dir <- file.path(output_dir, sample_name)
  dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)

  message("Processing sample: ", sample_name)

  # Dimensionality reduction
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = dims)
  seurat_obj <- RunTSNE(seurat_obj, dims = dims, perplexity = tsne_perplexity)
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)

  # Plot objects
  p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
    ggtitle(paste0("UMAP: Clusters - ", sample_name))
  p2 <- DimPlot(seurat_obj, reduction = "tsne", group.by = "seurat_clusters", label = TRUE) +
    ggtitle(paste0("t-SNE: Clusters - ", sample_name))
  # p3 <- DimPlot(seurat_obj, reduction = "umap", group.by = "doublet_status") +
  #   ggtitle(paste0("UMAP: Doublets - ", sample_name))
  p4 <- DimPlot(seurat_obj, reduction = "umap", group.by = "high_mito") +
    ggtitle(paste0("UMAP: High Mito - ", sample_name))
  p5 <- FeaturePlot(seurat_obj, features = "percent.mt", reduction = "umap", label = FALSE) +
    ggtitle(paste0("UMAP: % Mito Content - ", sample_name)) +
    scale_color_viridis_c()

  plot_list <- list(p1 = p1, p2 = p2, p4 = p4, p5 = p5)

  # Save individual plots
  for (plot_name in names(plot_list)) {
    pdf_path <- file.path(sample_dir, paste0(plot_name, "_", sample_name, ".pdf"))
    rds_path <- file.path(sample_dir, paste0(plot_name, "_", sample_name, ".rds"))
    ggsave(filename = pdf_path, plot = plot_list[[plot_name]], width = 7, height = 5)
    saveRDS(plot_list[[plot_name]], file = rds_path)
    message("Saved ", plot_name, " to: ", pdf_path)
  }

  return(seurat_obj)
}
