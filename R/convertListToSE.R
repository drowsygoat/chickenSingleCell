library(SummarizedExperiment)
library(Matrix)
library(tibble)
library(dplyr)

convertListToSE <- function(expr_list, assay_name = "counts") {
  stopifnot(length(expr_list) > 0)

  standardized_list <- list()

  for (sample_id in names(expr_list)) {
    mat <- expr_list[[sample_id]]

    if (inherits(mat, "dgCMatrix")) {
      mat <- as.matrix(mat)
    } else if (is.data.frame(mat)) {
      if (!is.null(rownames(mat))) {
        mat <- as.matrix(mat)
      } else {
        genes <- mat[[1]]
        mat <- as.matrix(mat[, -1])
        rownames(mat) <- genes
      }
    } else if (tibble::is_tibble(mat)) {
      genes <- mat[[1]]
      mat <- as.matrix(mat[, -1])
      rownames(mat) <- genes
    } else {
      warning(paste("Skipping unsupported matrix format in list element:", sample_id))
      next
    }

    if (anyDuplicated(rownames(mat))) {
      stop(paste("Duplicated feature names in sample:", sample_id))
    }

    colnames(mat) <- paste0(sample_id, "_", colnames(mat))
    standardized_list[[sample_id]] <- mat
  }

  if (length(standardized_list) == 0) {
    stop("No valid matrices found in the input list.")
  }

  # Auto-intersect feature names
  feature_sets <- lapply(standardized_list, rownames)
  common_features <- Reduce(intersect, feature_sets)

  if (length(common_features) == 0) {
    stop("No common features (genes) found across samples.")
  }

  standardized_list <- lapply(standardized_list, function(mat) mat[common_features, , drop = FALSE])
  combined_mat <- do.call(cbind, standardized_list)

  cluster_ids <- colnames(combined_mat)

  # Extract SampleID (e.g., ID10) and Cluster (e.g., C1) while ignoring parts like "obj" and "SCT"
  col_parts <- strsplit(cluster_ids, "_")
  sample_ids <- sapply(col_parts, function(x) x[1])
  cluster_labels <- sapply(col_parts, function(x) grep("^C[0-9]+$", x, value = TRUE)[1])

  coldata <- tibble(
    Sample = cluster_ids,
    SampleID = sample_ids,
    Cluster = cluster_labels
  )

  se <- SummarizedExperiment(
    assays = setNames(list(Matrix::Matrix(combined_mat, sparse = TRUE)), assay_name),
    colData = coldata,
    rowData = DataFrame(Feature = common_features)
  )
  return(se)
}