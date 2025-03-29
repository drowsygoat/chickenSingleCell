#' Run scDblFinder on a Seurat object and add assay-specific metadata
#'
#' This function runs doublet detection on a specified assay (e.g., "RNA" or "SCT")
#' using the \code{scDblFinder} package and stores the results in the Seurat object's metadata
#' under assay-specific column names.
#'
#' @param seurat_obj A \code{Seurat} object.
#' @param assay.type A character string specifying which assay to use (e.g., "RNA", "SCT").
#'
#' @return A \code{Seurat} object with added metadata:
#' \itemize{
#'   \item \code{doublet_score_<assay>}
#'   \item \code{predicted_doublet_<assay>}
#'   \item \code{doublet_status_<assay>}
#' }
#'
#' @importFrom scDblFinder scDblFinder
#' @importFrom SeuratObject as.SingleCellExperiment
#' @export
find_doublets_scDblFinder <- function(seurat_obj, assay.type = "RNA") {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object.")
  }

  if (!assay.type %in% names(seurat_obj@assays)) {
    stop("Assay '", assay.type, "' not found in Seurat object.")
  }

  message("Running scDblFinder on assay: ", assay.type)

  # Convert to SingleCellExperiment using specified assay
  sce <- as.SingleCellExperiment(seurat_obj, assay = assay.type)

  # Run scDblFinder
  sce <- scDblFinder(sce)

  # Extract metadata from SCE
  barcodes <- colnames(sce)
  doublet_score <- as.vector(colData(sce)$scDblFinder.score)
  predicted_doublet <- as.vector(colData(sce)$scDblFinder.class)
  doublet_status <- ifelse(predicted_doublet == "doublet", "Doublet", "Singlet")

  # Create metadata data.frame with appropriate column names
  meta <- data.frame(row.names = barcodes)
  meta[[paste0("doublet_score_", assay.type)]] <- doublet_score
  meta[[paste0("predicted_doublet_", assay.type)]] <- predicted_doublet
  meta[[paste0("doublet_status_", assay.type)]] <- doublet_status

  # Add to Seurat metadata
  seurat_obj <- AddMetaData(seurat_obj, metadata = meta)

  return(seurat_obj)
}
