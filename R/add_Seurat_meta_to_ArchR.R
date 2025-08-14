add_Seurat_meta_to_ArchR <- function(archr_proj, seurat_input,
                                     mode = c("intersect", "add"),
                                     absent_value = "absent",
                                     seurat_value_cols = NULL,
                                     dry_run = FALSE) {
  stopifnot(inherits(archr_proj, "ArchRProject"))

  mode <- match.arg(mode)
  message("Mode: ", mode, ifelse(mode == "add", " (adding Seurat metadata to ArchR)", " (subsetting ArchR to matched Seurat cells)"))

  # ── Step 1: Extract Seurat metadata ─────────────────────────────
  if (inherits(seurat_input, "Seurat")) {
    seurat_list <- list(seurat_input)
  } else if (is.list(seurat_input)) {
    # Filter out NULLs
    seurat_list <- Filter(Negate(is.null), seurat_input)
    if (!all(sapply(seurat_list, inherits, "Seurat"))) {
      stop("All non-NULL elements in the list must be Seurat objects.")
    }
  } else {
    stop("`seurat_input` must be a Seurat object or a list of Seurat objects.")
  }

  if (length(seurat_list) == 0) {
    stop("No valid Seurat objects found in input list.")
  }

  seurat_meta_list <- lapply(seurat_list, function(seurat_obj) {
    df <- data.frame(cell = Seurat::Cells(seurat_obj), seurat_obj@meta.data, stringsAsFactors = FALSE)
    df$sample  <- sub("_.*", "", df$cell)
    df$barcode <- sub(".*?_", "", df$cell)
    df
  })

  seurat_df <- do.call(rbind, seurat_meta_list)
  rownames(seurat_df) <- seurat_df$cell

  # Select metadata columns from Seurat
  if (is.null(seurat_value_cols)) {
    value_cols <- setdiff(colnames(seurat_df), c("cell", "sample", "barcode", "orig.ident"))
  } else {
    stopifnot(all(seurat_value_cols %in% colnames(seurat_df)))
    value_cols <- seurat_value_cols
  }

  # ── Step 2: Prepare ArchR metadata ──────────────────────────────
  archr_df <- as.data.frame(archr_proj@cellColData)
  archr_df$sample  <- sub("#.*", "", rownames(archr_df))
  archr_df$barcode <- sub(".*#", "", rownames(archr_df))

  # ── Step 3: Filter to shared samples ────────────────────────────
  common_samples <- intersect(unique(seurat_df$sample), unique(archr_df$sample))
  if (length(common_samples) == 0) {
    warning("❗ No overlapping samples between Seurat and ArchR. Returning NULL.")
    return(NULL)
  }

  seurat_df <- seurat_df[seurat_df$sample %in% common_samples, ]
  archr_df  <- archr_df[archr_df$sample %in% common_samples, ]

  # ── Step 4: Merge ───────────────────────────────────────────────
  merged <- merge(archr_df, seurat_df[, c("sample", "barcode", value_cols)],
                  by = c("sample", "barcode"), all.x = (mode == "add"))

  if (nrow(merged) == 0) {
    warning("❗ No matching cells found between ArchR and Seurat. Returning NULL.")
    return(NULL)
  }

  merged$cell <- paste0(merged$sample, "#", merged$barcode)
  rownames(merged) <- merged$cell

  # ── Step 5: Merge stats ─────────────────────────────────────────
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

  # ── Step 6: Add metadata to ArchRProject ────────────────────────
  new_meta <- merged[, value_cols, drop = FALSE]
  new_meta[is.na(new_meta)] <- absent_value

  if (mode == "intersect") {
    archr_proj <- subsetCells(archr_proj, cellNames = rownames(new_meta))
  }

  for (col in value_cols) {
    archr_proj@cellColData[[col]] <- NA
    matching_cells <- intersect(rownames(archr_proj@cellColData), rownames(new_meta))
    archr_proj@cellColData[matching_cells, col] <- new_meta[matching_cells, col]
  }

  return(archr_proj)
}
