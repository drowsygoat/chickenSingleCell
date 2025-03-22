#' Detect Doublets in a Seurat Object Using scDblFinder
#'
#' This function detects doublets in a Seurat object using the \code{scDblFinder} package.
#' It converts the Seurat object to a \code{SingleCellExperiment}, performs doublet detection,
#' and annotates the Seurat metadata with doublet scores, classifications, and a simplified status.
#'
#' Doublets are not removed, but are labeled in the metadata column \code{doublet_status}
#' as either "Singlet" or "Doublet".
#'
#' @param seurat_obj A \code{Seurat} object, typically after QC filtering.
#' @param assay.type Character string specifying the assay to use (default: "RNA").
#'
#' @return A \code{Seurat} object with additional metadata columns:
#' \itemize{
#'   \item \code{doublet_score} – numeric score from scDblFinder.
#'   \item \code{predicted_doublet} – classification ("doublet" or "singlet").
#'   \item \code{doublet_status} – simplified label ("Doublet" or "Singlet").
#' }
#'
#' @examples
#' seurat_obj <- detect_doublets(seurat_obj)
#'
#' @importFrom scDblFinder scDblFinder
#' @importFrom SeuratObject as.SingleCellExperiment
#' @export
detect_doublets <- function(seurat_obj, assay.type = "RNA") {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object.")
  }

  if (!assay.type %in% Assays(seurat_obj)) {
    stop("Assay '", assay.type, "' not found in Seurat object.")
  }

  message("Running scDblFinder on: ", seurat_obj@project.name)

  # Convert Seurat object to SingleCellExperiment
  sce <- as.SingleCellExperiment(seurat_obj, assay = assay.type)

  # Run scDblFinder
  sce <- scDblFinder(sce)

  # Transfer results back to Seurat metadata
  meta <- seurat_obj@meta.data
  meta$doublet_score <- sce$scDblFinder.score
  meta$predicted_doublet <- sce$scDblFinder.class
  meta$doublet_status <- ifelse(meta$predicted_doublet == "doublet", "Doublet", "Singlet")

  seurat_obj@meta.data <- meta

  return(seurat_obj)
}
