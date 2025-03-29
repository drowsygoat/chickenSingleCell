#' Label High Mitochondrial Content Cells in a Seurat Object
#'
#' This function labels cells with high mitochondrial gene expression in a specified assay
#' of a Seurat object. It calculates the percentage of mitochondrial genes if not already present,
#' and adds a metadata column flagging cells with mitochondrial content above the specified threshold.
#'
#' @param seurat_obj A \code{Seurat} object.
#' @param threshold A numeric value (default = 10). Cells with mitochondrial percentage above this threshold will be labeled as "High".
#' @param assay.type A character string specifying which assay to use (e.g., "RNA", "SCT").
#'
#' @return A \code{Seurat} object with an additional metadata column:
#' \itemize{
#'   \item \code{percent.mt_<assay>} if not already present
#'   \item \code{high_mito_<assay>} with values "High" or "Normal"
#' }
#'
#' @importFrom Seurat PercentageFeatureSet AddMetaData Cells
#' @export
label_high_mito <- function(seurat_obj, threshold = 10, assay.type = "RNA") {
  
  # Ensure input is a Seurat object
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object.")
  }
  
  # Create assay-specific percent.mt metadata if it doesn't exist
  percent_mt_col <- paste0("percent.mt_", assay.type)
  if (!percent_mt_col %in% colnames(seurat_obj@meta.data)) {
    seurat_obj[[percent_mt_col]] <- PercentageFeatureSet(
      seurat_obj, 
      pattern = "^ND1|^ND3|^ND4|^ND5|^ND6|^ATP6$|^ATP8$|^CYTB|^COII|^COX3", 
      assay = assay.type
    )
  }

  # Extract cell names (barcodes)
  barcodes <- Cells(seurat_obj)
  
  # Determine high mito status based on threshold
  high_mito_status <- ifelse(seurat_obj[[percent_mt_col]][, 1] > threshold, "High", "Normal")
  
  # Create metadata dataframe with assay-specific column name
  meta <- data.frame(row.names = barcodes)
  meta[[paste0("high_mito_", assay.type)]] <- high_mito_status
  
  # Add to Seurat object
  seurat_obj <- AddMetaData(seurat_obj, metadata = meta)
  
  return(seurat_obj)
}
