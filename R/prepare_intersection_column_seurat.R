prepare_intersection_column_seurat <- function(
  seurat_obj,
  columns,
  n = NULL,
  intersection_colname = "intersection"
) {
  stopifnot(inherits(seurat_obj, "Seurat"))
  stopifnot(length(columns) == 2)

  md <- seurat_obj@meta.data

  if (!all(columns %in% colnames(md))) {
    stop("Missing metadata columns: ",
         paste(setdiff(columns, colnames(md)), collapse = ", "))
  }

  # Build intersection (handles factors, NAs become "NA")
  intersection <- interaction(md[[columns[1]]], md[[columns[2]]],
                              sep = "_", drop = TRUE)
  md[[intersection_colname]] <- as.character(intersection)

  # Attach intersection
  seurat_obj <- SeuratObject::AddMetaData(
    seurat_obj,
    metadata = setNames(list(md[[intersection_colname]]), intersection_colname)
  )

  n_before <- ncol(seurat_obj)

  # If no filtering requested
  if (is.null(n)) {
    n_groups <- length(unique(md[[intersection_colname]]))
    message(sprintf("Total cells: %d | groups: %d", n_before, n_groups))
    return(seurat_obj)
  }

  # Count and filter by threshold
  counts <- table(md[[intersection_colname]])
  valid_groups <- names(counts[counts >= n])

  if (length(valid_groups) == 0L) {
    message(
      sprintf("No '%s' groups have at least n = %d cells. Returning NULL.",
              intersection_colname, n)
    )
    return(NULL)
  }

  keep_cells <- rownames(md)[md[[intersection_colname]] %in% valid_groups]

  if (length(keep_cells) == 0L) {
    message("After applying the threshold, keep_cells is empty. Returning NULL.")
    return(NULL)
  }

  # Subset and reattach the intersection
  seurat_obj <- subset(seurat_obj, cells = keep_cells)
  seurat_obj[[intersection_colname]] <- md[Cells(seurat_obj), intersection_colname, drop = TRUE]

  n_after <- ncol(seurat_obj)
  n_groups <- length(unique(seurat_obj[[intersection_colname]][,1]))

  message(sprintf(
    "Cells before: %d | after: %d | groups preserved: %d",
    n_before, n_after, n_groups
  ))

  return(seurat_obj)
}
