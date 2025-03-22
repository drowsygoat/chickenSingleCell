#' Normalize RNA Assay Using SCTransform
#'
#' This function applies Seurat's \code{SCTransform} normalization method
#' to the RNA assay of a Seurat object. It performs variance stabilization
#' and normalization based on regularized negative binomial regression.
#'
#' This step should be applied after quality control and doublet removal,
#' and before dimensionality reduction or integration.
#'
#' @param seurat_obj A \code{Seurat} object containing an RNA assay.
#'
#' @return A \code{Seurat} object with the normalized data stored in the "SCT" assay.
#'
#' @examples
#' seurat_obj <- normalize_with_sct(seurat_obj)
#'
#' @export
normalize_with_sct <- function(seurat_obj) {
  message("Running SCTransform on RNA assay...")
  seurat_obj <- SCTransform(seurat_obj, assay = "RNA", verbose = TRUE)
  return(seurat_obj)
}