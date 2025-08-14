prepare_intersection_column_seurat <- function(seurat_obj, columns, n = NULL, intersection_colname = "intersection") {
  stopifnot(inherits(seurat_obj, "Seurat"))
  stopifnot(length(columns) == 2)
  
  md <- seurat_obj@meta.data
  
  if (!all(columns %in% colnames(md))) {
    stop("Missing metadata columns: ", paste(setdiff(columns, colnames(md)), collapse = ", "))
  }

  md[[intersection_colname]] <- paste(md[[columns[1]]], md[[columns[2]]], sep = "_")

  if (!is.null(n)) {
    counts <- table(md[[intersection_colname]])
    valid_groups <- names(counts[counts >= n])
    keep_cells <- rownames(md)[md[[intersection_colname]] %in% valid_groups]
    seurat_obj <- subset(seurat_obj, cells = keep_cells)
  }

  seurat_obj@meta.data[[intersection_colname]] <- md[rownames(seurat_obj@meta.data), intersection_colname]
  return(seurat_obj)
}