#' Integrate Multiple Seurat Objects with Mitochondrial Regression
#'
#' A universal function to preprocess and integrate multiple Seurat objects,
#' optionally regressing out mitochondrial content instead of filtering it out.
#' Supports flexible quality control, normalization, and integration via CCA or RPCA.
#'
#' @param seurat_list A list of \code{Seurat} objects to integrate.
#' @param min_features Minimum number of detected features per cell (default = 200).
#' @param mito_threshold Threshold for annotating high mitochondrial cells (default = 10).
#' @param assay_type Assay name to use for normalization and regression (default = "RNA").
#' @param regress_mito Logical; whether to regress out mitochondrial percentage (default = TRUE).
#' @param regress_vars Character vector of metadata columns to regress out during scaling.
#' @param nfeatures Number of variable features to select (default = 3000).
#' @param dims Number of PCA dimensions to use for integration (default = 1:30).
#' @param anchor_reduction Reduction method for finding anchors ("cca" or "rpca").
#' @param return_unintegrated Logical; if TRUE, returns merged but unintegrated object (default = FALSE).
#' @param diet Logical; whether to slim Seurat objects using DietSeurat (default = TRUE).
#' @param diet_args A named list of arguments passed to \code{DietSeurat} if \code{diet = TRUE}.
#' @param verbose Logical; whether to print progress messages (default = TRUE).
#'
#' @return An integrated \code{Seurat} object with PCA run on the integrated assay,
#'         or a merged unintegrated object if \code{return_unintegrated = TRUE}.
#'
#' @importFrom Seurat RenameCells NormalizeData FindVariableFeatures ScaleData
#' @importFrom Seurat FindIntegrationAnchors IntegrateData RunPCA DietSeurat
#' @export
IntegrateSeuratObjects <- function(
  seurat_list,
  min_features = 200,
  mito_threshold = 10,
  assay_type = "RNA",
  regress_mito = TRUE,
  regress_vars = c("nFeature_RNA"),
  nfeatures = 2000,
  dims = 1:30,
  anchor_reduction = "cca",
  return_unintegrated = FALSE,
  diet = TRUE,
  diet_args = list(assays = "RNA", counts = TRUE, data = FALSE, scale.data = FALSE),
  verbose = TRUE
) {
  require(Seurat)
  require(dplyr)

  if (length(seurat_list) < 2) {
    stop("Integration requires at least two Seurat objects.")
  }

  # Use pbapply::pblapply if available, otherwise fall back to lapply
  apply_fun <- if (requireNamespace("pbapply", quietly = TRUE)) {
    function(X, FUN) pbapply::pblapply(X, FUN)
  } else {
    function(X, FUN) lapply(X, FUN)
  }

  # Preprocess each object
  seurat_list <- apply_fun(seq_along(seurat_list), function(i) {

    obj <- seurat_list[[i]]

    DefaultAssay(obj) <- assay_type

    # Optionally slim object
    if (diet) {
      obj <- do.call(DietSeurat, c(list(obj), diet_args))
    }


    obj <- label_high_mito(obj, threshold = mito_threshold, assay.type = assay_type)

    percent_mt_col <- paste0("percent.mt_", assay_type)
    if (regress_mito && !(percent_mt_col %in% regress_vars)) {
      regress_vars <- unique(c(regress_vars, percent_mt_col))
    }

    obj <- subset(obj, subset = nFeature_RNA >= min_features)

    obj <- NormalizeData(obj, assay = assay_type)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures, assay = assay_type)
    obj <- ScaleData(obj, vars.to.regress = regress_vars, verbose = FALSE, assay = assay_type)

    return(obj)
  })

  # MEMORY - CAREFUL!!!!
  if (return_unintegrated) {
    if (verbose) message("Returning merged object without integration.")
    merged <- merge(x = seurat_list[[1]], y = seurat_list[-1])
    return(merged)
  }

  if (verbose) message("Finding integration anchors...")
  anchors <- FindIntegrationAnchors(
    object.list = seurat_list,
    normalization.method = "LogNormalize",
    anchor.features = nfeatures,
    reduction = anchor_reduction,
    dims = dims
  )

  if (verbose) message("Integrating data...")
  integrated <- IntegrateData(anchorset = anchors, dims = dims)

  DefaultAssay(integrated) <- "integrated"
  integrated <- ScaleData(integrated, verbose = FALSE)
  integrated <- RunPCA(integrated, verbose = FALSE)

  if (verbose) message("Integration complete.")
  return(integrated)
}


    # # Infer sample ID
    # sample_ids <- infer_sample_id(seurat_obj)
    # unique_samples <- unique(sample_ids)
    # if (length(unique_samples) != 1) {
    #   stop("Multiple samples detected: ", paste(unique_samples, collapse = ", "))
    # }

    # sample_name <- unique_samples