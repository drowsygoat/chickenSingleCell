# Function to randomly subset a fraction of cells from a Seurat object
subset_random_cells <- function(seurat_obj, fraction = 0.2, seed = 123) {
  # Validate inputs
  stopifnot(is(seurat_obj, "Seurat"), 
            is.numeric(fraction), fraction > 0, fraction <= 1, 
            is.numeric(seed))

  # Set seed for reproducibility
  set.seed(seed)

  # Get total number of cells
  total_cells <- length(Cells(seurat_obj))

  # Compute number of cells to sample
  num_cells <- round(fraction * total_cells)

  # Randomly sample cell names
  selected_cells <- sample(Cells(seurat_obj), num_cells)

  # Subset Seurat object
  seurat_subset <- subset(seurat_obj, cells = selected_cells)

  return(seurat_subset)
}


# Example usage
# seurat_subset <- subset_random_cells(seurat_obj, fraction = 0.3)
