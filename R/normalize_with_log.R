#' Normalize RNA Assay Using LogNormalization
#'
#' This function applies Seurat's \code{NormalizeData} method with LogNormalization
#' to the RNA assay of a Seurat object. It optionally regresses out technical
#' covariates (e.g. mitochondrial content, batch) using \code{ScaleData}.
#'
#' Regressed variables should already be present in \code{seurat_obj@meta.data},
#' e.g. \code{percent.mt}, \code{batch}, or cell cycle scores.
#'
#' @param seurat_obj A \code{Seurat} object with raw counts in the "RNA" assay.
#' @param assay Character, name of the assay to normalize (default: "RNA")
#' @param vars.to.regress Character vector of metadata columns to regress out
#'        (default: c("percent.mt")); set to NULL to skip regression
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#'
#' @return A \code{Seurat} object with normalized and optionally scaled data,
#'         or \code{NULL} if input is malformed.
#'
#' @examples
#' seurat_obj <- normalize_with_log(seurat_obj, vars.to.regress = c("percent.mt", "batch"))
#'
#' @export
normalize_with_log <- function(seurat_obj,
                               assay = "RNA",
                               vars.to.regress = c("percent.mt"),
                               verbose = TRUE) {
  if (!inherits(seurat_obj, "Seurat")) {
    warning("Input must be a Seurat object. Returning NULL.")
    return(NULL)
  }

  if (!assay %in% names(seurat_obj@assays)) {
    stop(paste("Assay", assay, "not found in Seurat object."))
  }

  assay_obj <- seurat_obj[[assay]]
  counts_data <- tryCatch({
    if (inherits(assay_obj, "Assay5")) {
      assay_obj@layers[["counts"]]
    } else {
      assay_obj@counts
    }
  }, error = function(e) NULL)

  if (!inherits(counts_data, "Matrix") && !is.matrix(counts_data)) {
    warning("Assay '", assay, "' does not contain a proper 'counts' matrix. Skipping normalization.")
    return(NULL)
  }

  DefaultAssay(seurat_obj) <- assay

  if (!is.null(vars.to.regress)) {
    missing_vars <- setdiff(vars.to.regress, colnames(seurat_obj@meta.data))
    if (length(missing_vars) > 0) {
      stop("The following vars.to.regress are missing from metadata: ",
           paste(missing_vars, collapse = ", "))
    }
  }

  if (verbose) message("Running LogNormalization on assay: ", assay)

  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize",
                              layer = if (inherits(assay_obj, "Assay5")) "counts" else NULL,
                              verbose = verbose)

  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = verbose)

  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = vars.to.regress, verbose = verbose)

  return(seurat_obj)
}
