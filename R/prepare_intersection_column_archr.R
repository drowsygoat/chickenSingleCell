prepare_intersection_column_archr <- function(archr_proj, columns, n = NULL, intersection_colname = "intersection") {
  stopifnot(inherits(archr_proj, "ArchRProject"))
  stopifnot(length(columns) == 2)
  
  md <- as.data.frame(archr_proj@cellColData)

  if (!all(columns %in% colnames(md))) {
    stop("Missing metadata columns: ", paste(setdiff(columns, colnames(md)), collapse = ", "))
  }

  md[[intersection_colname]] <- paste(md[[columns[1]]], md[[columns[2]]], sep = "_")

  if (!is.null(n)) {
    counts <- table(md[[intersection_colname]])
    valid_groups <- names(counts[counts >= n])
    keep_cells <- rownames(md)[md[[intersection_colname]] %in% valid_groups]
    archr_proj <- subsetCells(archr_proj, cellNames = keep_cells)
  }

  archr_proj@cellColData[[intersection_colname]] <- md[rownames(archr_proj@cellColData), intersection_colname]
  return(archr_proj)
}