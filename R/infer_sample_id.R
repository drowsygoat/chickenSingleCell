#' Infer Sample ID from a Seurat Object
#'
#' Attempts to infer sample identity for each cell in a Seurat object.
#' First checks metadata columns (`orig.ident`, `sample.*`, `Sample.*`, `SAMPLE.*`).
#' If not found, attempts to parse the sample name from cell name prefixes using `_` or `#` as separators.
#'
#' @param seurat_obj A Seurat object.
#'
#' @return A character vector of sample IDs, one per cell.
#' @export
infer_sample_id <- function(seurat_obj) {
  meta <- seurat_obj@meta.data
  sample_col <- grep("^orig\\.ident$|^sample.*$|^Sample.*$|^SAMPLE.*$", colnames(meta), value = TRUE)
  
  if (length(sample_col) == 1) {
    return(meta[[sample_col]])
  } else {
    # Try to extract sample name from cell names
    cell_names <- colnames(seurat_obj)
    if (all(grepl("_", cell_names))) {
      return(sapply(strsplit(cell_names, "_"), `[`, 1))
    } else if (all(grepl("#", cell_names))) {
      return(sapply(strsplit(cell_names, "#"), `[`, 1))
    } else {
      stop("Sample ID could not be inferred from metadata or cell names.")
    }
  }
}
