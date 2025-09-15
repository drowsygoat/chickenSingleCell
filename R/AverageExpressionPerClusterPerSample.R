AverageExpressionPerClusterPerSample <- function(obj_list,
                                                 assay = NULL,
                                                 layer = "data",
                                                 cluster_col = "seurat_clusters",
                                                 scale.factor = 10000,
                                                 return.se = FALSE,
                                                 se.save.path = NULL) {

  if (inherits(obj_list, "Seurat")) {
    obj_list <- list(obj_list)
  }

  result_list <- unlist(
    lapply(seq_along(obj_list), function(i) {
      obj <- obj_list[[i]]

      if (is.null(obj)) {
        message(sprintf("Object %d is NULL. Skipping.", i))
        return(setNames(list(NULL), paste0("obj", i)))
      }

      this_assay <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay
      message(sprintf("Object %d: using assay = '%s'", i, this_assay))

      sample_ids <- infer_sample_id(obj)
      obj[["sample_id"]] <- sample_ids

      if (length(unique(sample_ids)) > 1) {
        message(sprintf("Object %d has multiple samples. Will compute per-sample aggregates.", i))
      }

      Seurat::Idents(obj) <- cluster_col
      clusters <- levels(obj)
      message(sprintf("Object %d has %d clusters.", i, length(clusters)))

      sample_names <- unique(obj$sample_id)

      lapply(sample_names, function(sample) {
        sample_obj <- subset(obj, subset = sample_id == sample)

        # ✅ Safety check: skip if too few cells
        if (ncol(sample_obj) < 2) {
          message(sprintf("⚠️ Skipping sample '%s' in object %d: too few cells (%d)", sample, i, ncol(sample_obj)))
          return(setNames(list(NULL), paste0(sample, "_obj", i, "_", this_assay)))
        }

        # ✅ Drop empty cluster levels before pseudobulk
        if (!cluster_col %in% colnames(sample_obj@meta.data)) {
          stop(sprintf("Cluster column '%s' not found in object %d", cluster_col, i))
        }

        # Set idents from the column (works for numeric/character/factor)
        cl_vals <- sample_obj[[cluster_col]][, 1, drop = TRUE]
        Seurat::Idents(sample_obj) <- cl_vals

        # Drop empty levels only if it's a factor
        if (is.factor(cl_vals)) cl_vals <- droplevels(cl_vals)

        # Count unique, non-NA clusters
        n_nonempty_clusters <- length(unique(cl_vals[!is.na(cl_vals)]))

        if (n_nonempty_clusters == 0L) {
        message(sprintf("⚠️ Skipping sample '%s' in object %d: no non-empty clusters.", sample, i))
        return(setNames(list(NULL), paste0(sample, "_obj", i, "_", this_assay)))
        }

        if (!this_assay %in% names(sample_obj@assays)) {
          stop(sprintf("Assay '%s' not found in object %d", this_assay, i))
        }

        if (all(is.na(sample_obj[[cluster_col]]))) {
          warning(sprintf("All cluster values are NA in object %d (sample %s)", i, sample))
          return(NULL)
        }

        message(sprintf("Processing sample: %s (object %d)", sample, i))
        message(sprintf("  Cells: %d", ncol(sample_obj)))
        message(sprintf("  Assay: %s, Cluster col: %s", this_assay, cluster_col))

        message("  Cluster table:")
        print(table(sample_obj[[cluster_col]]))

        message("  Sample ID table:")
        print(table(sample_obj$sample_id))

        name <- paste0(sample, "_obj", i, "_", this_assay)

        result <- tryCatch({
          agg_raw <- Seurat::PseudobulkExpression(
            sample_obj,
            assays = this_assay,
            layer = layer,
            group.by = cluster_col,
            scale.factor = scale.factor,
            return.seurat = FALSE,
            verbose = FALSE
          )

          message("  ✅ PseudobulkExpression succeeded")
          agg_exp <- agg_raw[[this_assay]]

          # 🔒 Guard: empty result (no clusters produced columns)
          if (is.null(agg_exp) || ncol(agg_exp) == 0) {
            message(sprintf("  ⚠️ No clusters produced columns for sample '%s' in object %d. Skipping.", sample, i))
            return(setNames(list(NULL), name))
          }

          setNames(list(agg_exp), name)
        }, error = function(e) {
          message("  ❌ PseudobulkExpression failed:")
          message(e$message)
          message("  🔍 Inspecting structure of sample_obj and metadata:")
          print(str(sample_obj@meta.data))
          setNames(list(NULL), name)
        })

        return(result)
      })
    }),
    recursive = FALSE
  )

  # If you use purrr elsewhere, keep this; otherwise base R is fine.
  result_list <- purrr::list_flatten(result_list)

  if (!return.se) {
    return(result_list)
  }

  # Convert matrices to SummarizedExperiment objects
  se_list <- lapply(names(result_list), function(name) {
    mat <- result_list[[name]]
    if (is.null(mat) || ncol(mat) == 0) return(NULL)  # 🔒 guard against 0-column matrices

    # ✅ Fallback for missing gene names
    if (is.null(rownames(mat)) || any(rownames(mat) == "")) {
      warning(sprintf("⚠️ Matrix '%s' has missing gene names. Assigning fallback names.", name))
      rownames(mat) <- paste0("gene_", seq_len(nrow(mat)))
    }

    parts <- strsplit(name, "_")[[1]]
    sample <- parts[1]
    assay_name <- parts[3]
    cluster_ids <- colnames(mat)

    new_colnames <- paste(sample, cluster_ids, assay_name, sep = "_")
    colnames(mat) <- new_colnames

    col_data <- data.frame(
      cluster = cluster_ids,
      sample_id = sample,
      assay = assay_name,
      row.names = new_colnames,
      stringsAsFactors = FALSE
    )

    row_data <- data.frame(
      gene = rownames(mat),
      row.names = rownames(mat),
      stringsAsFactors = FALSE
    )

    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = mat),
      colData = col_data,
      rowData = row_data
    )
  })

  se_list <- Filter(Negate(is.null), se_list)

  if (length(se_list) == 0) {
    warning("No valid SummarizedExperiments to merge. Returning NULL.")
    return(NULL)
  }

  all_genes <- unique(unlist(lapply(se_list, SummarizedExperiment::rownames)))

  pad_se <- function(se, all_genes) {
    mat <- SummarizedExperiment::assay(se)

    if (!inherits(mat, "dgCMatrix")) {
      mat <- Matrix::as(mat, "dgCMatrix")
    }

    missing_genes <- setdiff(all_genes, rownames(mat))

    if (length(missing_genes) > 0) {
      zero_mat <- Matrix::Matrix(
        0, nrow = length(missing_genes), ncol = ncol(mat),
        dimnames = list(missing_genes, colnames(mat)), sparse = TRUE
      )
      mat <- rbind(mat, zero_mat)
    }

    mat <- mat[all_genes, , drop = FALSE]

    col_data <- as.data.frame(SummarizedExperiment::colData(se))
    row_data <- data.frame(gene = all_genes, row.names = all_genes, stringsAsFactors = FALSE)

    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = mat),
      colData = col_data,
      rowData = row_data
    )
  }

  se_list_padded <- lapply(se_list, pad_se, all_genes = all_genes)

  merged_se <- tryCatch({
    do.call(SummarizedExperiment::cbind, se_list_padded)
  }, error = function(e) {
    stop("❌ Failed to merge SummarizedExperiment objects: ", e$message)
  })

  if (!is.null(se.save.path)) {
    dir.create(dirname(se.save.path), showWarnings = FALSE, recursive = TRUE)
    saveRDS(merged_se, se.save.path)
    message(sprintf("📁 Saved SummarizedExperiment to: %s", se.save.path))
  }

  if (return.se && !is.null(merged_se)) {
    covariate_dir_mean <- file.path(dirname(se.save.path), "cluster_summary_covariates_mean")
    covariate_dir_median <- file.path(dirname(se.save.path), "cluster_summary_covariates_median")

    dir.create(covariate_dir_mean, showWarnings = FALSE, recursive = TRUE)
    dir.create(covariate_dir_median, showWarnings = FALSE, recursive = TRUE)

    mat <- SummarizedExperiment::assay(merged_se)
    col_data <- as.data.frame(SummarizedExperiment::colData(merged_se))

    col_data$cluster <- as.character(col_data$cluster)
    col_data$sample_id <- as.character(col_data$sample_id)

    all_clusters <- sort(unique(col_data$cluster))
    all_samples <- sort(unique(col_data$sample_id))

    for (clust in all_clusters) {
      cols_in_cluster <- which(col_data$cluster == clust)
      sub_mat <- mat[, cols_in_cluster, drop = FALSE]
      sub_coldata <- col_data[cols_in_cluster, , drop = FALSE]
      sample_groups <- split(seq_along(cols_in_cluster), sub_coldata$sample_id)

      median_vals <- sapply(all_samples, function(sid) {
        if (!sid %in% names(sample_groups)) return(NA_real_)
        vals <- as.numeric(sub_mat[, sample_groups[[sid]], drop = FALSE])
        vals <- vals[vals != 0]
        if (length(vals) == 0) return(NA_real_)
        stats::median(vals, na.rm = TRUE)
      })

      mean_vals <- sapply(all_samples, function(sid) {
        if (!sid %in% names(sample_groups)) return(NA_real_)
        vals <- as.numeric(sub_mat[, sample_groups[[sid]], drop = FALSE])
        vals <- vals[vals != 0]
        if (length(vals) == 0) return(NA_real_)
        mean(vals, na.rm = TRUE)
      })

      median_df <- as.data.frame(t(median_vals))
      rownames(median_df) <- paste0("cluster_", clust)
      median_file <- file.path(covariate_dir_median, paste0("cluster_", clust, "_median.txt"))
      median_df_out <- cbind(V1 = rownames(median_df), median_df)
      utils::write.table(median_df_out, file = median_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

      mean_df <- as.data.frame(t(mean_vals))
      rownames(mean_df) <- paste0("cluster_", clust)
      mean_file <- file.path(covariate_dir_mean, paste0("cluster_", clust, "_mean.txt"))
      mean_df_out <- cbind(V1 = rownames(mean_df), mean_df)
      utils::write.table(mean_df_out, file = mean_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

      message(sprintf("✅ Saved cluster summary: %s (mean + median)", clust))
    }
  }

  return(merged_se)
}
