#' Calculate Aggregated (Pseudobulk) Expression per Cluster per Sample
#'
#' Computes pseudobulk (summed) expression for each cluster within each sample.
#' Handles either a single Seurat object or a list of Seurat objects.
#' Uses `AggregateExpression()` and supports RNA, SCTransform, or other assays.
#'
#' @param obj_list A Seurat object or a list of Seurat objects.
#' @param assay The assay to use (e.g., `"RNA"` or `"SCT"`). If NULL, uses the object's active assay.
#' @param slot The expression slot to use. Default is `"data"`.
#' @param cluster_col Name of the metadata column containing cluster identities. Default is `"seurat_clusters"`.
#' @param normalization.method Method for normalization. Options:
#'   \describe{
#'     \item{"LogNormalize"}{(default) Divide feature counts by total counts per cell, multiply by `scale.factor`, then `log1p`}
#'     \item{"CLR"}{Centered log-ratio transformation}
#'     \item{"RC"}{Relative counts; divide by total counts per cell, multiply by `scale.factor`, no log}
#'     \item{"none" or NULL}{No normalization (recommended for SCT residuals)}
#'   }
#' @param scale.factor Scaling factor used in normalization. Default is `10000`.
#' @param ... Additional arguments passed to `AggregateExpression()`.
#'
#' @return A named list of matrices. Each matrix contains aggregated gene expression (genes x clusters) for a sample.
#' @export
#'
#' @examples
#' # Default log-normalization (RNA)
#' result <- AverageExpressionPerClusterPerSample(seurat_obj)
#'
#' # Pseudobulk counts without normalization (e.g., for SCT)
#' result <- AverageExpressionPerClusterPerSample(seurat_obj, assay = "SCT", normalization.method = "none")
#'
#' # Counts per million
#' result <- AverageExpressionPerClusterPerSample(seurat_obj, normalization.method = "RC", scale.factor = 1e6)
AverageExpressionPerClusterPerSample <- function(obj_list,
                                                 assay = NULL,
                                                 slot = "data",
                                                 cluster_col = "seurat_clusters",
                                                 normalization.method = c("LogNormalize", "CLR", "RC", "none"),
                                                 scale.factor = 10000,
                                                 ...) {
                                                  
  normalization.method <- match.arg(normalization.method)

  if (inherits(obj_list, "Seurat")) {
    obj_list <- list(obj_list)
  }

  result_list <- unlist(
    lapply(seq_along(obj_list), function(i) {
      obj <- obj_list[[i]]

      if (is.null(obj)) {
        message(sprintf("Object %d is NULL. Skipping.", i))
        return(setNames(list(NULL), paste0("obj", i)))
      }

      this_assay <- if (is.null(assay)) DefaultAssay(obj) else assay
      message(sprintf("Object %d: using assay = '%s'", i, this_assay))

      sample_ids <- infer_sample_id(obj)
      obj[["sample_id"]] <- sample_ids
      
      if (length(unique(sample_ids)) > 1) {
        message(sprintf("Object %d has multiple samples. Will compute per-sample aggregates.", i))
      }

      Idents(obj) <- cluster_col
      clusters <- levels(obj)
      message(sprintf("Object %d has %d clusters.", i, length(clusters)))

      # Determine normalization
      norm_method <- if (tolower(normalization.method) %in% c("none", "null")) {
        message(sprintf("Skipping normalization for object %d (normalization.method = 'none')", i))
        NULL
      } else {
        normalization.method
      }

      sample_names <- unique(obj$sample_id)

      lapply(sample_names, function(sample) {
        sample_obj <- subset(obj, subset = sample_id == sample)
        if (length(unique(sample_obj[[cluster_col]])) < 1) return(NULL)

             if (!this_assay %in% Assays(sample_obj)) {
        stop(sprintf("Assay '%s' not found in object %d", this_assay, i))
      }

      if (!cluster_col %in% colnames(sample_obj@meta.data)) {
        stop(sprintf("Cluster column '%s' not found in object %d", cluster_col, i))
      }

      if (all(is.na(sample_obj[[cluster_col]]))) {
        warning(sprintf("All cluster values are NA in object %d (sample %s)", i, sample))
        return(NULL)
      }

        agg_exp <- AggregateExpression(
          sample_obj,
          assays = this_assay,
          slot = slot,
          group.by = cluster_col,
          normalization.method = norm_method,
          scale.factor = scale.factor,
          return.seurat = FALSE,
          verbose = FALSE
        )[[this_assay]]

        name <- paste0(sample, "_obj", i, "_", this_assay)
        setNames(list(agg_exp), name)
      })
    }),
    recursive = FALSE
  )

  return(flatten(result_list))
  
}