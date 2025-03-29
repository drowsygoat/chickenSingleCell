#' Cluster and Explore a Single Seurat Object
#'
#' Performs PCA, UMAP, t-SNE, clustering, and metadata visualization for a single sample.
#'
#' @param seurat_obj Seurat object.
#' @param dims PCA dimensions to use (default: 1:30).
#' @param resolution Clustering resolution (default: 0.5).
#' @param use.assay Assay to use (e.g. "SCT", "RNA") (default: "SCT").
#' @param tsne_perplexity Perplexity for t-SNE (default: 30).
#' @param output_dir Output directory (default: "plots").
#' @param sample_name Optional sample name; defaults to project name.
#' @param subset_fraction Fraction of cells to retain (0–1; default: 1).
#' @param subset_seed Random seed (default: 42).
#' @param save_metadata_csv If TRUE, saves metadata as CSV (default: TRUE).
#' @param save_heatmap If TRUE, saves DimHeatmap of top PCs (default: TRUE).
#'
#' @return Processed Seurat object.
#' @export
cluster_and_explore <- function(seurat_obj,
                                dims = 1:30,
                                resolution = 1,
                                use.assay = "SCT",
                                tsne_perplexity = 30,
                                output_dir = "plots",
                                sample_name = NULL,
                                subset_fraction = 1.0,
                                subset_seed = 42,
                                save_metadata_csv = TRUE,
                                save_heatmap = TRUE) {

  require(Seurat)
  require(ggplot2)
  require(cowplot)
  require(viridis)
  require(grid)

  DefaultAssay(seurat_obj) <- use.assay
  
  if (is.null(sample_name)) sample_name <- seurat_obj@project.name

  seurat_obj <- subset_seurat(seurat_obj, subset_fraction, subset_seed)

  sample_dir <- file.path(output_dir, sample_name)
  individual_dir <- file.path(sample_dir, "plots_individual")
  dir.create(individual_dir, showWarnings = FALSE, recursive = TRUE)

  message("📦 Processing sample: ", sample_name)

  # ---- PCA ----
  message("🧠 Running PCA...")
  seurat_obj <- RunPCA(seurat_obj, npcs = max(dims), verbose = FALSE)
  if (save_heatmap) {
    pdf(file = "dummy.pdf")
    p_heatmap <- DimHeatmap(seurat_obj, dims = dims[1:9], cells = 500, balanced = TRUE)
    dev.off()
    ggsave(file.path(individual_dir, paste0("heatmap_top_pcs_", sample_name, ".pdf")), p_heatmap, width = 7, height = 6)
  }

  umap_name <- paste0("umap_", tolower(use.assay))
  tsne_name <- paste0("tsne_", tolower(use.assay))

  # ---- UMAP and tSNE ----
  message("🔀 Running UMAP...")
  seurat_obj <- RunUMAP(seurat_obj, dims = dims, reduction.name = umap_name, reduction.key = "UMAP_")

  message("🔀 Running t-SNE...")
  seurat_obj <- RunTSNE(seurat_obj, dims = dims, perplexity = tsne_perplexity,
                        reduction.name = tsne_name, reduction.key = "TSNE_")

  # ---- Clustering ----
  message("🔗 Finding neighbors...")
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims, verbose = FALSE)

  message("🧩 Finding clusters...")
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution, verbose = FALSE)

  cluster_label <- paste0("clusters_", tolower(use.assay))
  seurat_obj[[cluster_label]] <- seurat_obj$seurat_clusters

  # ---- Metadata export ----
  if (save_metadata_csv) {
    metadata_path <- file.path(sample_dir, paste0("metadata_", sample_name, ".csv"))
    write.csv(seurat_obj@meta.data, file = metadata_path)
    message("📄 Saved metadata: ", metadata_path)
  }

  # ---- Plotting Setup ----
  plot_list <- list()
  cluster_plots <- list()

  # Cluster plots
  for (red in c(umap_name, tsne_name)) {
    short <- gsub(".*_", "", red)
    dim_short <- gsub("_.*", "", red)
    key <- paste0("clusters_", dim_short, "_", short)

    p <- DimPlot(seurat_obj, reduction = red, group.by = cluster_label, label = TRUE) +
      ggtitle(paste(toupper(dim_short), ": Clusters -", sample_name, "(", use.assay, ")"))

    cluster_plots[[key]] <- p
    plot_list[[key]] <- p
  }

  # Metadata plots
  suffix_pattern <- paste0("_", use.assay, "$")
  meta_matches <- grep(suffix_pattern, colnames(seurat_obj@meta.data), value = TRUE)

  for (var_full in meta_matches) {
    var_base <- sub(suffix_pattern, "", var_full)
    is_cat <- is.factor(seurat_obj[[var_full]][[1]]) || is.character(seurat_obj[[var_full]][[1]])

    for (red in c(umap_name, tsne_name)) {
      dim_short <- gsub("_.*", "", red)
      red_short <- gsub(".*_", "", red)

      plot_key <- paste0("meta_", var_base, "_", dim_short, "_", red_short)
      title <- paste(toupper(dim_short), ":", var_base, "(", use.assay, ") -", sample_name)

      p <- if (is_cat) {
        DimPlot(seurat_obj, reduction = red, group.by = var_full, label = TRUE) + ggtitle(title)
      } else {
        FeaturePlot(seurat_obj, features = var_full, reduction = red) +
          ggtitle(title) +
          scale_color_viridis_c()
      }

      plot_list[[plot_key]] <- p
    }
  }

  # Save all individual plots
  for (plot_name in names(plot_list)) {
    base_filename <- paste0(plot_name, "_", sample_name)
    pdf_path <- file.path(individual_dir, paste0(base_filename, ".pdf"))
    rds_path <- file.path(individual_dir, paste0(base_filename, ".rds"))
    ggsave(pdf_path, plot_list[[plot_name]], width = 7, height = 5)
    saveRDS(plot_list[[plot_name]], file = rds_path)
    message("📸 Saved: ", pdf_path)
  }

  # ---- Violin + Marker Plots ----
  vln_plots <- plot_vln_metadata_by_assay(seurat_obj,
                                          assay = use.assay,
                                          cluster_col = cluster_label,
                                          sample_name = sample_name,
                                          return_plots = TRUE)

  # 🐛 Debug: Print cluster sizes
  cluster_sizes <- table(seurat_obj[[cluster_label]][, 1])
  message("🔎 Cluster sizes (", cluster_label, "):")
  print(cluster_sizes)

  # Safety check
  if (length(unique(Idents(seurat_obj))) <= 1) {
    warning("⚠️ Only one cluster found — skipping marker gene detection.")
    marker_plots <- list()
  } else {
    # Run marker detection with error capture
    marker_plots <- tryCatch({
      plot_top_marker_genes(seurat_obj = seurat_obj,
                            cluster_col = cluster_label,
                            use.assay = use.assay,
                            output_dir = output_dir,
                            sample_name = sample_name,
                            umap_name = umap_name,
                            return_plots = TRUE)
    }, error = function(e) {
      message("❌ ERROR in plot_top_marker_genes: ", conditionMessage(e))
      return(list())
    })
  }

  marker_plots <- plot_top_marker_genes(seurat_obj = seurat_obj,
                                        cluster_col = cluster_label,
                                        use.assay = use.assay,
                                        output_dir = output_dir,
                                        sample_name = sample_name,
                                        umap_name = umap_name,
                                        return_plots = TRUE)

  # ---- Helper: Save dynamically sized plot section ----
  save_plot_section <- function(plots, file_prefix, sample_name, assay, ncol = 2, plot_width = 5, plot_height = 5) {
    n <- length(plots)
    if (n == 0) return(NULL)
    nrow <- ceiling(n / ncol)
    width <- plot_width * ncol
    height <- plot_height * nrow

    pdf_path <- file.path(sample_dir, paste0(file_prefix, "_", sample_name, "_", assay, ".pdf"))
    pdf(pdf_path, width = width, height = height)
    print(plot_grid(plotlist = plots, ncol = ncol))
    dev.off()
    message("✅ Saved: ", pdf_path)
  }


  # ---- Assign metadata plots ----
  meta_keys  <- names(plot_list)[grepl("^meta_", names(plot_list))]
  meta_plots <- plot_list[meta_keys]

  cat(meta_keys)

  # ---- Save all sections as separate, dynamically sized PDFs ----
  save_plot_section(cluster_plots, "Cluster_Plots", sample_name, use.assay, ncol = 2)
  save_plot_section(meta_plots, "Metadata_Plots", sample_name, use.assay, ncol = 2)
  save_plot_section(vln_plots, "Violin_Plots", sample_name, use.assay, ncol = 2)
  save_plot_section(marker_plots, "Marker_Plots", sample_name, use.assay, ncol = 3)

  message("✅ All plots saved for: ", sample_name)

  return(seurat_obj)
}
