#' Randomly subset a fraction of cells from a Seurat object
#'
#' @param seurat_obj A Seurat object
#' @param fraction Numeric (0 < fraction ≤ 1): fraction of cells to keep
#' @param seed Random seed for reproducibility
#'
#' @return A subsetted Seurat object (or original if fraction = 1 or <100 cells), or NULL if input is not a Seurat object
#' @export
subset_seurat <- function(seurat_obj, fraction = 1.0, seed = 42) {
  # Check input class
  if (!inherits(seurat_obj, "Seurat")) {
    warning("❌ Input is not a Seurat object. Returning NULL.")
    return(NULL)
  }

  # Validate other inputs
  stopifnot(
    is.numeric(fraction), fraction > 0, fraction <= 1,
    is.numeric(seed)
  )

  total_cells <- length(Cells(seurat_obj))

  if (total_cells < 100) {
    message("⚠️ Fewer than 100 cells (", total_cells, "). Returning full object.")
    return(seurat_obj)
  }

  if (fraction < 1) {
    set.seed(seed)
    num_cells <- floor(fraction * total_cells)
    selected_cells <- sample(Cells(seurat_obj), num_cells)
    seurat_obj <- subset(seurat_obj, cells = selected_cells)
    message("✅ Subset to ", num_cells, " cells (", round(fraction * 100, 1), "%)")
  } else {
    message("Fraction = 1: returning full Seurat object (", total_cells, " cells)")
  }

  return(seurat_obj)
}
