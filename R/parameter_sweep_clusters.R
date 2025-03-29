#' Run clustering parameter sweep and visualize with clustree
#'
#' @param seurat_obj Seurat object with PCA run (or other reduction).
#' @param dims_list List of PC ranges (e.g. list(1:10, 1:20))
#' @param resolutions Vector of resolutions for FindClusters
#' @param k_list Vector of k.param values (for FindNeighbors)
#' @param reduction Dimensionality reduction to use (default: "pca")
#' @param clustree_prefix Prefix used in metadata columns for clustering (default: "sweep_")
#' @param clustree_plot_file Optional path to save clustree plot (PDF or PNG)
#'
#' @return List of Seurat objects (one per tested setting)
#'
#' @export
parameter_sweep_clusters <- function(seurat_obj,
                                     dims_list = list(1:10, 1:30, 1:50),
                                     resolutions = c(0.4, 1.2, 2),
                                     k_list = c(20, 40, 60),
                                     reduction = "pca",
                                     clustree_prefix = "sweep_",
                                     clustree_plot_file = NULL) {

  if (!requireNamespace("clustree", quietly = TRUE)) {
    stop("Please install the 'clustree' package.")
  }

  results <- list()
  original_idents <- Idents(seurat_obj)

  for (dims in dims_list) {
    dims_str <- paste0("dims", min(dims), "-", max(dims))
    for (k in k_list) {
      for (res in resolutions) {
        id <- paste(dims_str, paste0("k", k), paste0("res", res), sep = "_")

        message("Running: ", id)

        obj <- seurat_obj
        if (!reduction %in% names(obj@reductions)) {
          obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)
        }

        obj <- FindNeighbors(obj,
                             reduction = reduction,
                             dims = dims,
                             k.param = k,
                             verbose = FALSE)

        obj <- FindClusters(obj,
                            resolution = res,
                            algorithm = 1,
                            verbose = FALSE)

        cluster_col <- paste0(clustree_prefix, id)
        obj[[cluster_col]] <- obj$seurat_clusters

        Idents(obj) <- original_idents

        results[[id]] <- obj
      }
    }
  }

  # Take one Seurat object and add all clusterings as metadata
  seurat_combined <- seurat_obj
  for (id in names(results)) {
    colname <- paste0(clustree_prefix, id)
    seurat_combined[[colname]] <- results[[id]][[colname]]
  }

  # Plot clustree
  # clustree_plot <- clustree::clustree(seurat_combined, prefix = clustree_prefix)
  # print(clustree_plot)

  # Optionally save
  if (!is.null(clustree_plot_file)) {
    message("Saving clustree plot to: ", clustree_plot_file)
    ggsave(clustree_plot_file, plot = clustree_plot, width = 10, height = 8)
  }

  return(invisible(results))
}
