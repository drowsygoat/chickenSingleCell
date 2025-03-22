

#' Create a Seurat Object from 10X Multiome Data with Quality Filtering
#'
#' @param sample_dir Path to the 10X sample directory.
#' @param gtf_file Path to the GTF file for gene annotations.
#' @param dry_run Logical, if TRUE, checks for required files without loading data.
#' @param min_counts_atac Minimum ATAC counts per cell (default: 3000).
#' @param min_nFeature_RNA Minimum unique genes per cell (default: 1500).
#' @param remove_failed_samples Logical, if TRUE (default), samples with no passing cells are removed.
#' @param include_RNA Logical, if TRUE (default), includes RNA data in the Seurat object.
#' @param include_ATAC Logical, if TRUE (default), includes ATAC data in the Seurat object.
#' @param include_metadata Logical, if TRUE (default), adds metadata to the Seurat object.
#' @return A Seurat object with RNA and/or ATAC assays, filtered by quality metrics.
#' @export
import_to_seurat <- function(sample_dir, gtf_file, dry_run = FALSE,
                             min_counts_atac = 2000,
                             min_nFeature_RNA = 1000,
                             remove_failed_samples = TRUE,
                             include_RNA = TRUE,
                             include_ATAC = TRUE,
                             include_metadata = FALSE) {

  if (!dir.exists(sample_dir)) {
    stop("The specified sample directory does not exist: ", sample_dir)
  }

  message("Processing sample: ", sample_dir)

  # Define paths
  h5_path <- file.path(sample_dir, "outs", "filtered_feature_bc_matrix.h5")
  atac_fragments <- file.path(sample_dir, "outs", "atac_fragments.tsv.gz")
  metadata_path <- file.path(sample_dir, "outs", "per_barcode_metrics.csv")

  # Define sample name
  sample_name <- basename(sample_dir)

  # Dry run: Check file existence
  if (dry_run) {
    message("Dry run enabled. Checking for required files:")
    message("- RNA counts file: ", ifelse(file.exists(h5_path), "FOUND", "NOT FOUND"))
    message("- ATAC fragments file: ", ifelse(file.exists(atac_fragments), "FOUND", "NOT FOUND"))
    message("- Metadata file: ", ifelse(file.exists(metadata_path), "FOUND", "NOT FOUND"))
    message("- GTF file: ", ifelse(file.exists(gtf_file), "FOUND", "NOT FOUND"))
    return(invisible(NULL))
  }

  # Load RNA and ATAC counts
  message("Loading RNA and ATAC counts...")
  counts <- Read10X_h5(h5_path)

  if (!"Gene Expression" %in% names(counts) || !"Peaks" %in% names(counts)) {
    stop("Expected 'Gene Expression' and 'Peaks' data in the H5 file.")
  }

  rna_counts <- counts$`Gene Expression`
  atac_counts <- counts$Peaks

  # Create Seurat object
  seurat_obj <- NULL

  if (include_RNA) {
    message("Creating Seurat object for RNA assay...")
    seurat_obj <- CreateSeuratObject(counts = rna_counts, assay = "RNA", project = sample_name)
  }

  # Import GTF as GRanges
  if (!file.exists(gtf_file)) {
    stop("GTF file not found at: ", gtf_file)
  }

  message("Importing GTF file for gene annotations...")
  gtf <- rtracklayer::import(gtf_file, format = "gtf")
  gene.coords <- gtf[gtf$type == "gene"]
  mcols(gene.coords) <- mcols(gene.coords)[, colSums(!is.na(mcols(gene.coords))) > 0]

  # Add ATAC data
  if (include_ATAC && file.exists(atac_fragments)) {
    message("Adding ATAC assay...")
    atac_assay <- CreateChromatinAssay(
      counts = atac_counts,
      sep = c(":", "-"),
      fragments = atac_fragments,
      annotation = gene.coords
    )

    if (is.null(seurat_obj)) {
      seurat_obj <- CreateSeuratObject(assays = list(ATAC = atac_assay))
    } else {
      seurat_obj[["ATAC"]] <- atac_assay
    }
  } else if (include_ATAC) {
    warning("ATAC fragments file not found: ", atac_fragments)
  }

  # Add metadata if available
  if (include_metadata && file.exists(metadata_path)) {
    message("Adding metadata...")
    metadata <- read.csv(metadata_path, row.names = 1)
    seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)
  } else if (include_metadata) {
    warning("Metadata file not found: ", metadata_path)
  }

  # Apply Quality Control Filtering
  if (!is.null(seurat_obj)) {
    message("Applying quality control filtering...")

    if ("nFeature_RNA" %in% colnames(seurat_obj@meta.data) &&
        "nCount_ATAC" %in% colnames(seurat_obj@meta.data)) {

      cells_to_keep <- which(
        seurat_obj$nFeature_RNA >= min_nFeature_RNA &
        seurat_obj$nCount_ATAC >= min_counts_atac
      )

      if (length(cells_to_keep) == 0) {
        if (remove_failed_samples) {
          message("No cells passed filtering. Sample will be removed.")
          return(NULL)
        } else {
          warning("No cells passed filtering. Returning unfiltered Seurat object.")
        }
      } else {
        seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
        message("Seurat object created successfully with ", ncol(seurat_obj), " high-quality cells.")
      }

    } else {
      warning("Quality metrics not found in Seurat object metadata. Skipping filtering.")
    }
  }

  # Rename cells to include sample name
  seurat_obj <- RenameCells(seurat_obj, add.cell.id = sample_name)

  return(seurat_obj)
}

# ADD SCT transforrm and save
# return list of surat object and corresponding ggplots