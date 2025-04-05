IntegrateSeuratObjects <- function(
      seurat_list,
      mito_threshold = 10,
      min_features = 200,
      assay_type = "RNA",
      regress_mito = TRUE,
      regress_vars = c("nFeature_RNA"),
      nFeatures = 2000,
      dims = 1:30,
      integration_method = "cca",  # One of: cca, rpca, harmony, mnn, scvi
      return_unintegrated = FALSE,
      diet = TRUE,
      diet_args = list(assays = "RNA", counts = TRUE, data = FALSE, scale.data = FALSE),
      k.weight = 100,
      log_warnings = TRUE,
      verbose = TRUE
    ) {
  require(Seurat)
  require(dplyr)

  if (length(seurat_list) < 2) {
    stop("Integration requires at least two Seurat objects.")
  }

  apply_fun <- if (requireNamespace("pbapply", quietly = TRUE)) {
    function(X, FUN) pbapply::pblapply(X, FUN)
  } else {
    function(X, FUN) lapply(X, FUN)
  }

  seurat_list <- Filter(function(x) inherits(x, "Seurat"), seurat_list)

  seurat_list <- apply_fun(seq_along(seurat_list), function(i) {
    obj <- seurat_list[[i]]

    if (ncol(obj) < 200) return(NULL)

    DefaultAssay(obj) <- assay_type

    if (diet) {
      obj <- do.call(DietSeurat, c(list(obj), diet_args))
    }

    obj <- label_high_mito(obj, threshold = mito_threshold, assay.type = assay_type)

    percent_mt_col <- paste0("percent.mt_", assay_type)

    if (regress_mito && !(percent_mt_col %in% regress_vars)) {
      regress_vars <- unique(c(regress_vars, percent_mt_col))
    }

    obj <- subset(obj, subset = nFeature_RNA >= min_features)

    obj <- NormalizeData(obj, assay = assay_type)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nFeatures)
    obj <- ScaleData(obj, vars.to.regress = regress_vars, verbose = FALSE)

    return(obj)
  })

  seurat_list <- Filter(Negate(is.null), seurat_list)
  if (length(seurat_list) < 2) stop("At least two valid Seurat objects are required.")

  if (return_unintegrated) return(do.call(merge, seurat_list))

  if (verbose) message("Merging and running PCA...")
  merged <- do.call(merge, seurat_list)
  cat("abc")
  merged <- RunPCA(merged, npcs = max(dims), verbose = FALSE)

  # Choose integration method
  if (verbose) message("Running integration using method: ", integration_method)
  method_map <- list(
    cca = CCAIntegration,
    rpca = RPCAIntegration,
    harmony = HarmonyIntegration,
    mnn = SeuratWrappers::FastMNNIntegration,
    scvi = SeuratWrappers::scVIIntegration
  )
  if (!integration_method %in% names(method_map)) {
    stop("Unsupported integration method: ", integration_method)
  }

  integration_args <- list(
    object = merged,
    method = method_map[[integration_method]],
    verbose = verbose
  )

  if (integration_method %in% c("cca", "rpca", "harmony")) {
    integration_args$orig.reduction <- "pca"
    integration_args$new.reduction <- paste0("integrated.", integration_method)
  } else if (integration_method == "mnn") {
    integration_args$new.reduction <- "integrated.mnn"
  } else if (integration_method == "scvi") {
    integration_args$new.reduction <- "integrated.scvi"
  }

  obj <- do.call(IntegrateLayers, integration_args)

  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, reduction = integration_args$new.reduction, verbose = FALSE)

  if (log_warnings) {
    w <- warnings()
    if (length(w) > 0) {
      message("\u26A0\uFE0F Integration warnings:")
      print(w)
    }
  }

  if (verbose) message("Integration complete.")
  
  return(obj)
}
