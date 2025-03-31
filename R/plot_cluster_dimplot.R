plot_cluster_dimplot <- function(seurat_obj,
                                 group_by,
                                 highlight = NULL,
                                 reduction = NULL,
                                 label = TRUE,
                                 pt.size = 1,
                                 pdf_dir = "plots",
                                 save_as_pdf = TRUE) {

  make_error_plot <- function(msg, sample_name = "unknown") {
    p <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::geom_text(
        ggplot2::aes(x = 0.5, y = 0.5, label = msg),
        size = 6, hjust = 0.5, vjust = 0.5
      ) +
      ggplot2::ggtitle(paste0("Error - ", sample_name))
    return(p)
  }

  tryCatch({
    # Check Seurat object
    if (!inherits(seurat_obj, "Seurat")) {
      stop("Input is not a Seurat object.")
    }

    # Check metadata column
    if (!group_by %in% colnames(seurat_obj@meta.data)) {
      stop(paste0("Metadata column '", group_by, "' not found in Seurat object."))
    }

    # Infer sample ID
    sample_ids <- infer_sample_id(seurat_obj)
    unique_samples <- unique(sample_ids)
    if (length(unique_samples) != 1) {
      stop("Multiple samples detected: ", paste(unique_samples, collapse = ", "))
    }
    sample_name <- unique_samples

    # Determine reductions
    all_reductions <- Seurat::Reductions(seurat_obj)
    if (is.null(reduction)) {
      reduction <- all_reductions
      message("📌 No reduction specified. Using all: ", paste(reduction, collapse = ", "))
    } else {
      missing <- setdiff(reduction, all_reductions)
      if (length(missing) > 0) {
        stop("Reduction(s) not found: ", paste(missing, collapse = ", "))
      }
    }

    # Create temp metadata column
    is_continuous <- is.numeric(seurat_obj@meta.data[[group_by]])
    temp_col <- paste0("highlight_", group_by)
    seurat_tmp <- seurat_obj

    if (is_continuous || is.null(highlight)) {
      seurat_tmp[[temp_col]] <- seurat_tmp@meta.data[[group_by]]
    } else {
      seurat_tmp[[temp_col]] <- ifelse(seurat_tmp@meta.data[[group_by]] %in% highlight,
                                       seurat_tmp@meta.data[[group_by]], "Other")
      seurat_tmp[[temp_col]] <- factor(seurat_tmp[[temp_col]],
                                       levels = c(sort(unique(highlight)), "Other"))
    }

    # Generate plots
    plot_list <- list()
    for (red in reduction) {
      p <- if (is_continuous) {
        Seurat::FeaturePlot(
          seurat_tmp,
          features = temp_col,
          reduction = red,
          pt.size = pt.size
        )
      } else {
        Seurat::DimPlot(
          seurat_tmp,
          group.by = temp_col,
          reduction = red,
          label = label,
          pt.size = pt.size
        )
      }

      p <- p +
        ggplot2::labs(
          title = paste0(red, " - ", group_by, " (", sample_name, ")"),
          color = group_by
        ) +
        ggplot2::theme(aspect.ratio = 1)

      plot_list[[red]] <- p
    }

    # Save to PDF
    if (save_as_pdf && length(plot_list) > 0) {
      if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
      tag <- if (is.null(highlight)) "all" else paste(highlight, collapse = "_")
      filename <- file.path(pdf_dir, paste0("DimPlot_", sample_name, "_", group_by, "_", tag, ".pdf"))

      grDevices::pdf(filename, width = 20, height = 15)
      for (i in seq(1, length(plot_list), by = 4)) {
        patch <- patchwork::wrap_plots(plot_list[i:min(i + 3, length(plot_list))])
        print(patch)
      }
      grDevices::dev.off()
      message("💾 Saved PDF to: ", filename)
    }

    # Return plots
    if (length(plot_list) == 1) {
      return(invisible(plot_list[[1]]))
    } else {
      combined <- patchwork::wrap_plots(plot_list)
      return(invisible(combined))
    }

  }, error = function(e) {
    msg <- paste("❌", conditionMessage(e))
    sample_name <- tryCatch({
      sample_ids <- infer_sample_id(seurat_obj)
      unique(sample_ids)
    }, error = function(e2) {
      "unknown"
    })

    err_plot <- make_error_plot(msg, sample_name = sample_name)

    if (save_as_pdf) {
      if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
      filename <- file.path(pdf_dir, paste0("DimPlot_", sample_name, "_", group_by, "_ERROR.pdf"))
      grDevices::pdf(filename, width = 8, height = 8)
      print(err_plot)
      grDevices::dev.off()
      message("💾 Saved error PDF to: ", filename)
    }

    return(invisible(err_plot))
  })
}
