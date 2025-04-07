#' Merge ArchR Metadata into a Seurat Object
#'
#' Merges one or more metadata columns from an ArchR-exported dataframe (e.g., cellColData)
#' into a Seurat object. Matches cells by sample and barcode, handling different ID separators ("_" vs "#").
#'
#' @param seurat_obj A Seurat object.
#' @param archr_df A data.frame with ArchR cell names as rownames and metadata columns.
#' @param mode "intersect" keeps only matched cells, "add" preserves all Seurat cells.
#' @param absent_value Value assigned to unmatched cells in "add" mode.
#' @param archr_value_cols Optional. Character vector of column names in `archr_df` to merge. If NULL, merges all columns.
#' @param dry_run If TRUE, only prints merge stats and doesn't modify the Seurat object.
#'
#' @return Modified Seurat object (or NULL if no matches found / dry_run = TRUE)
add_ArchR_meta_to_Seurat <- function(seurat_obj, archr_df,
                                     mode = c("intersect", "add"),
                                     absent_value = "absent",
                                     archr_value_cols = NULL,
                                     dry_run = FALSE) {
  mode <- match.arg(mode)

  # Parse Seurat cells
  seurat_df <- data.frame(cell = Seurat::Cells(seurat_obj), stringsAsFactors = FALSE)
  seurat_df$sample  <- sub("_.*", "", seurat_df$cell)
  seurat_df$barcode <- sub(".*?_", "", seurat_df$cell)

  # Parse ArchR rownames
  archr_ids <- rownames(archr_df)
  archr_df$sample  <- sub("#.*", "", archr_ids)
  archr_df$barcode <- sub(".*#", "", archr_ids)

  # Select metadata columns to merge
  if (is.null(archr_value_cols)) {
    value_cols <- setdiff(colnames(archr_df), c("sample", "barcode"))
  } else {
    stopifnot(all(archr_value_cols %in% colnames(archr_df)))
    value_cols <- archr_value_cols
  }

  # Filter to common samples
  common_samples <- intersect(unique(seurat_df$sample), unique(archr_df$sample))
  if (length(common_samples) == 0) {
    warning("❗ No overlapping samples between Seurat and ArchR. Returning NULL.")
    return(NULL)
  }
  seurat_df <- seurat_df[seurat_df$sample %in% common_samples, ]
  archr_df  <- archr_df[archr_df$sample %in% common_samples, ]

  # Merge metadata
  merged <- merge(seurat_df, archr_df[, c("sample", "barcode", value_cols)],
                  by = c("sample", "barcode"), all.x = (mode == "add"))
  if (nrow(merged) == 0) {
    warning("❗ No matching cells found between Seurat and ArchR. Returning NULL.")
    return(NULL)
  }

  merged$cell <- paste0(merged$sample, "_", merged$barcode)
  rownames(merged) <- merged$cell

  # Merge stats
  stats <- aggregate(!is.na(merged[[value_cols[1]]]), by = list(Sample = merged$sample), FUN = function(m) {
    total <- length(m)
    matched <- sum(m)
    sprintf("%d/%d matched (%.1f%%)", matched, total, 100 * matched / total)
  })
  names(stats)[2] <- "Matched"

  if (dry_run) {
    message("✅ Dry run - merge summary per sample:")
    print(stats, row.names = FALSE)
    return(invisible(NULL))
  }

  message("✅ Merge summary per sample:")
  print(stats, row.names = FALSE)

  # Final metadata: all selected value_cols
  new_meta <- merged[, c("cell", value_cols), drop = FALSE]
  rownames(new_meta) <- new_meta$cell
  new_meta$cell <- NULL
  new_meta[is.na(new_meta)] <- absent_value

  # Subset Seurat if intersect mode
  if (mode == "intersect") {
    seurat_obj <- subset(seurat_obj, cells = rownames(new_meta))
  }

  # Add new metadata
  seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = new_meta)

  return(seurat_obj)
}
