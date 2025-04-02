# 

NOW WE WORK IN SOMETHING  ESLSE


#' Plot Integrated Seurat Dimensionality Reductions with Features or Metadata
#'
#' Generates `FeaturePlot` or `DimPlot` visualizations for dimensionality reduction results
#' (e.g., UMAP, t-SNE) from a Seurat object, grouped by metadata or gene expression.
#' Supports subsetting, highlighting, automatic feature selection from DE results, and PDF export.
#'
#' @param seurat_obj A Seurat object.
#' @param group_by Character. Name of the metadata column to group/color cells by. If `NULL` and `feature` is also `NULL`, defaults to `Idents(seurat_obj)`.
#' @param highlight Character vector of values in `group_by` to highlight. All other values will be grouped under "Other".
#' @param feature Character. Gene name (or vector of genes) to plot expression with `FeaturePlot`. If `NULL`, uses `group_by` metadata instead.
#' @param auto_feature_from_de Logical. If `TRUE` and `feature = NULL`, attempts to auto-select a top DE gene (based on `de_source`).
#' @param cluster_id Character or numeric. Cluster ID to filter DE genes by when selecting feature automatically.
#' @param n_genes Integer. Number of top DE genes to plot when `auto_feature_from_de = TRUE` and `cluster_id` is provided.
#' @param de_source Character. Where to extract DE results from, e.g., "misc$de_results".
#' @param reduction Character vector. Which reductions to plot (e.g., "umap", "tsne"). If `NULL`, all reductions are used.
#' @param split_by Character. Metadata column to facet/split plots by.
#' @param pt.size Numeric. Point size for cells.
#' @param subset_expr Expression. A logical expression (using metadata columns) to subset the object before plotting, e.g., `percent.mt < 10`.
#' @param pdf_dir Character. Directory where PDF will be saved if `save_as_pdf = TRUE`.
#' @param save_as_pdf Logical. Whether to save the plots to a PDF.
#' @param color_gradient Character vector of length 2. Colors for gradient if plotting a numeric metadata column.
#' @param show_legend Logical. Whether to display the legend in the output plots.
#' @param verbose Logical. If `TRUE`, prints informative messages during execution.
#'
#' @return A `ggplot` object (if 1 plot), or a `patchwork` object combining multiple plots. Also optionally saves to PDF.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_integrated_dimplot(seurat_obj)
#' plot_integrated_dimplot(seurat_obj, group_by = "cell_type")
#' plot_integrated_dimplot(seurat_obj, feature = "Sox2")
#' plot_integrated_dimplot(seurat_obj, feature = "Sox2", split_by = "sample")
#' plot_integrated_dimplot(seurat_obj, group_by = "cell_type", highlight = c("Neuron", "Astrocyte"))
#' plot_integrated_dimplot(seurat_obj, group_by = "cell_type", subset_expr = percent.mt < 10)
#' plot_integrated_dimplot(seurat_obj, auto_feature_from_de = TRUE, cluster_id = "1", n_genes = 3)
#' }
plot_integrated_dimplot2 <- function(seurat_obj,
                                    group_by = NULL,
                                    highlight = NULL,
                                    feature = NULL,
                                    auto_feature_from_de = FALSE,
                                    cluster_id = NULL,
                                    n_genes = 1,
                                    de_source = "misc$de_results",
                                    reduction = NULL,
                                    split_by = NULL,
                                    pt.size = 1,
                                    subset_expr = NULL,
                                    pdf_dir = "plot_integrated_dimplot2",
                                    save_as_pdf = TRUE,
                                    color_gradient = c("lightgrey", "blue"),
                                    show_legend = TRUE,
                                    verbose = TRUE) {

  make_error_plot <- function(msg) {
    ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::geom_text(
        ggplot2::aes(x = 0.5, y = 0.5, label = msg),
        size = 6, hjust = 0.5, vjust = 0.5
      ) +
      ggplot2::ggtitle("Error")
  }

  if (!auto_feature_from_de && !is.null(cluster_id)) {
    warning("⚠️ 'cluster_id' is ignored because auto_feature_from_de = FALSE.")
  }

  tryCatch({
    if (!inherits(seurat_obj, "Seurat")) stop("Input is not a Seurat object.")
    if (verbose) message("✅ Seurat object check passed")

    if (!is.null(subset_expr)) {
      if (verbose) message("🔍 Subsetting Seurat object using expression: ", deparse(substitute(subset_expr)))
      tryCatch({
        seurat_obj <- subset(seurat_obj, subset = subset_expr)
      }, error = function(e) {
        stop("❌ Failed to subset Seurat object using `subset_expr`: ", conditionMessage(e))
      })
    }

    if (is.null(feature) && is.null(group_by)) {
      group_by <- "seurat_idents"
      seurat_obj[[group_by]] <- Idents(seurat_obj)
      if (verbose) message("📌 Defaulting to Idents(seurat_obj) as group_by.")
    }

    if (is.null(feature) && auto_feature_from_de) {
      if (verbose) message("🔍 Trying to auto-select feature(s) from DE results...")
      de_results <- tryCatch(eval(parse(text = paste0("seurat_obj@", de_source))), error = function(e) NULL)
      top_genes <- NULL

      if (!is.null(de_results) && "gene" %in% colnames(de_results)) {
        if (!is.null(cluster_id) && "cluster" %in% colnames(de_results)) {
          top_genes <- unique(de_results[de_results$cluster == cluster_id, "gene"])[1:n_genes]
          if (verbose) message("🔎 Using top ", n_genes, " DE gene(s) from cluster ", cluster_id, ": ", paste(top_genes, collapse = ", "))
        } else {
          top_genes <- head(unique(de_results$gene), n_genes)
          if (verbose) message("✨ Using top ", n_genes, " DE gene(s) from DE table: ", paste(top_genes, collapse = ", "))
        }
      }

      if (is.null(top_genes) && length(VariableFeatures(seurat_obj)) > 0) {
        top_genes <- head(VariableFeatures(seurat_obj), n_genes)
        if (verbose) message("⚠️ DE results not found. Using top variable gene(s): ", paste(top_genes, collapse = ", "))
      }

      if (!is.null(top_genes) && all(top_genes %in% rownames(seurat_obj))) {
        feature <- top_genes
      } else {
        warning("❌ Could not auto-select gene(s). You may need to set `feature` manually.")
      }
    }

    all_reductions <- Seurat::Reductions(seurat_obj)
    if (is.null(reduction)) {
      reduction <- all_reductions
      if (verbose) message("📌 No reduction specified. Using: ", paste(reduction, collapse = ", "))
    } else {
      missing <- setdiff(reduction, all_reductions)
      if (length(missing) > 0) stop("Reduction(s) not found: ", paste(missing, collapse = ", "))
      if (verbose) message("✅ Using specified reductions: ", paste(reduction, collapse = ", "))
    }

    plots <- list()
    split_levels <- if (!is.null(split_by)) unique(seurat_obj[[split_by]][, 1]) else NA
    if (verbose && !is.null(split_by)) message("📌 Split by: ", split_by, " with levels: ", paste(split_levels, collapse = ", "))

    for (red in reduction) {
      if (verbose) message("🔄 Processing reduction: ", red)

      genes_to_plot <- if (!is.null(feature)) as.character(feature) else NULL

      if (!is.null(genes_to_plot)) {
        for (plot_feature in genes_to_plot) {
          if (!(plot_feature %in% rownames(seurat_obj))) next

          if (is.na(split_levels[1])) {
            p <- Seurat::FeaturePlot(
              seurat_obj,
              features = plot_feature,
              reduction = red,
              pt.size = pt.size
            ) +
              ggplot2::labs(title = paste0(red, " - ", plot_feature), color = "Expression") +
              ggplot2::theme(aspect.ratio = 1)
            if (!show_legend) p <- p + ggplot2::theme(legend.position = "none")
            plots[[paste0(red, "_", plot_feature)]] <- p
          } else {
            for (split_val in split_levels) {
              subset_cells <- colnames(seurat_obj)[seurat_obj[[split_by]][, 1] == split_val]
              subset_obj <- subset(seurat_obj, cells = subset_cells)
              p <- Seurat::FeaturePlot(
                subset_obj,
                features = plot_feature,
                reduction = red,
                pt.size = pt.size
              ) +
                ggplot2::labs(title = paste0(red, " - ", plot_feature, " (", split_by, ": ", split_val, ")"), color = "Expression") +
                ggplot2::theme(aspect.ratio = 1)
              if (!show_legend) p <- p + ggplot2::theme(legend.position = "none")
              plots[[paste0(red, "_", plot_feature, "_", split_val)]] <- p
            }
          }
        }
      } else {
        if (!group_by %in% colnames(seurat_obj@meta.data)) {
          stop("Column '", group_by, "' not found in metadata.")
        }

        is_numeric <- is.numeric(seurat_obj@meta.data[[group_by]])
        temp_col <- paste0("highlight_", group_by)
        meta_col <- if (is.data.frame(seurat_obj[[group_by]])) seurat_obj[[group_by]][, 1] else seurat_obj[[group_by]]
        seurat_tmp <- seurat_obj

        if (!is.null(highlight) && !is_numeric) {
          valid_highlights <- highlight[highlight %in% unique(meta_col)]
          if (verbose) message("✨ Highlighting: ", paste(valid_highlights, collapse = ", "))
          seurat_tmp[[temp_col]] <- ifelse(meta_col %in% valid_highlights, meta_col, "Other")
          seurat_tmp[[temp_col]] <- factor(seurat_tmp[[temp_col]], levels = c(sort(valid_highlights), "Other"))
        } else {
          seurat_tmp[[temp_col]] <- meta_col
        }

        if (is.na(split_levels[1])) {
          if (is_numeric) {
            p <- Seurat::FeaturePlot(
              seurat_tmp,
              features = group_by,
              reduction = red,
              pt.size = pt.size
            ) +
              ggplot2::labs(title = paste0(red, " - ", group_by), color = group_by) +
              ggplot2::theme(aspect.ratio = 1) +
              ggplot2::scale_color_gradient(low = color_gradient[1], high = color_gradient[2])
          } else {
            p <- Seurat::DimPlot(
              seurat_tmp,
              group.by = temp_col,
              reduction = red,
              pt.size = pt.size,
              label = TRUE
            ) +
              ggplot2::labs(title = paste0(red, " - ", group_by), color = group_by) +
              ggplot2::theme(aspect.ratio = 1)
          }
          if (!show_legend) p <- p + ggplot2::theme(legend.position = "none")
          plots[[red]] <- p
        } else {
          for (split_val in split_levels) {
            subset_cells <- colnames(seurat_tmp)[seurat_tmp[[split_by]][, 1] == split_val]
            subset_obj <- subset(seurat_tmp, cells = subset_cells)

            if (is_numeric) {
              p <- Seurat::FeaturePlot(
                subset_obj,
                features = group_by,
                reduction = red,
                pt.size = pt.size
              ) +
                ggplot2::labs(title = paste0(red, " - ", group_by, " (", split_by, ": ", split_val, ")"), color = group_by) +
                ggplot2::theme(aspect.ratio = 1) +
                ggplot2::scale_color_gradient(low = color_gradient[1], high = color_gradient[2])
            } else {
              p <- Seurat::DimPlot(
                subset_obj,
                group.by = temp_col,
                reduction = red,
                pt.size = pt.size,
                label = TRUE
              ) +
                ggplot2::labs(title = paste0(red, " - ", group_by, " (", split_by, ": ", split_val, ")"), color = group_by) +
                ggplot2::theme(aspect.ratio = 1)
            }
            if (!show_legend) p <- p + ggplot2::theme(legend.position = "none")
            plots[[paste0(red, "_", split_val)]] <- p
          }
        }
      }
    }

    if (save_as_pdf && length(plots) > 0) {
      if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
      tag <- if (!is.null(feature)) paste0(feature, collapse = "_") else if (!is.null(highlight)) paste(highlight, collapse = "_") else group_by
      suffix <- paste0("_", tag, if (!is.null(split_by)) paste0("_splitby_", split_by) else "")
      filename <- file.path(pdf_dir, paste0("IntegratedPlot", suffix, ".pdf"))

      plot_count <- length(plots)
      if (plot_count <= 3) {
        ncol <- plot_count
        plot_width <- 10
        plot_height <- 10
      } else {
        ncol <- 3
        plot_width <- 7
        plot_height <- 7
      }
      nrow <- ceiling(plot_count / ncol)

      patch <- patchwork::wrap_plots(plots, ncol = ncol)
      ggplot2::ggsave(filename, plot = patch, width = plot_width * ncol, height = plot_height * nrow, limitsize = FALSE)
      if (verbose) message("💾 Saved PDF to: ", filename)
    }

    if (length(plots) == 1) {
      return(invisible(plots[[1]]))
    } else {
      return(invisible(patchwork::wrap_plots(plots, ncol = 3)))
    }
  }, error = function(e) {
    msg <- paste("Error:", conditionMessage(e))
    message(msg)

    err_plot <- make_error_plot(msg)
    if (save_as_pdf) {
      if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
      filename <- file.path(pdf_dir, "IntegratedPlot_ERROR.pdf")
      ggplot2::ggsave(filename, err_plot, width = 8, height = 8)
      if (verbose) message("💾 Saved error PDF to: ", filename)
    }
    return(invisible(err_plot))
  })
}
