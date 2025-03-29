#' Normalize RNA Assay Using SCTransform with Technical Regressors
#'
#' This function applies Seurat's \code{SCTransform} normalization method
#' to the RNA assay of a Seurat object. It performs variance stabilization
#' using regularized negative binomial regression and optionally regresses out
#' technical covariates such as mitochondrial content and batch.
#'
#' Regressed variables should already be present in \code{seurat_obj@meta.data},
#' e.g. \code{percent.mt}, \code{batch}, or cell cycle scores.
#'
#' @param seurat_obj A \code{Seurat} object with raw counts in the "RNA" assay.
#' @param vars.to.regress Character vector of metadata columns to regress out (default: c("percent.mt"))
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#'
#' @return A \code{Seurat} object with SCTransform-normalized data stored in the "SCT" assay.
#'
#' @examples
#' seurat_obj <- normalize_with_sct(seurat_obj, vars.to.regress = c("percent.mt", "batch"))
#'
#' @export
normalize_with_sct <- function(seurat_obj,
                               vars.to.regress = c("percent.mt"),
                               verbose = TRUE) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object.")
  }

  if (!"RNA" %in% Assays(seurat_obj)) {
    stop("RNA assay not found in Seurat object.")
  }

  missing_vars <- setdiff(vars.to.regress, colnames(seurat_obj@meta.data))
  if (length(missing_vars) > 0) {
    stop("The following vars.to.regress are missing from metadata: ",
         paste(missing_vars, collapse = ", "))
  }

  message("Running SCTransform on RNA assay with regression of: ",
          paste(vars.to.regress, collapse = ", "))

  seurat_obj <- SCTransform(
    seurat_obj,
    assay = "RNA",
    vars.to.regress = vars.to.regress,
    verbose = verbose
  )

  return(seurat_obj)
}
