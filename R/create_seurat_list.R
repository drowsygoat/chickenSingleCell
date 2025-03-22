#' Create Seurat Objects from 10X Multiome Data
#'
#' This script processes a directory of 10X Genomics (Cell Ranger ARC) sample folders,
#' creates Seurat objects for each sample, and stores them in a list. Each Seurat object
#' includes both RNA and ATAC assays along with per-barcode metadata if available.
#'
#' The function expects the standard 10X output directory structure, including the filtered H5 matrix,
#' ATAC fragments file, and per-barcode metrics. It also requires a GTF file for gene annotations.
#'
#' @param parent_dir A character string specifying the path to the parent directory containing sample folders.
#' @param gtf_file A character string specifying the path to the GTF file for gene annotations.
#'
#' @return A list of Seurat objects, one for each sample directory.
#'
#' @examples
#' seurat_list <- create_seurat_list(
#'   parent_dir = "/path/to/parent_dir",
#'   gtf_file = "/path/to/annotations.gtf"
#' )
#'
#' @export
#' 

create_seurat_list <- function(parent_dir, gtf_file, dry_run = FALSE) {
  if (!dir.exists(parent_dir)) {
    stop("The specified parent directory does not exist: ", parent_dir)
  }

  message("Processing samples in directory: ", parent_dir)

  # List all sample directories
  sample_dirs <- list.dirs(parent_dir, recursive = FALSE)

  # Initialize list to store Seurat objects
  seurat_list <- list()

  # Process each sample directory
  for (sample in sample_dirs) {
    seurat_tmp <- import_to_seurat(sample, gtf_file, dry_run = dry_run)
    if (is.null(seurat_tmp)) {
      message("Skipping ", sample, " due to no cells after subsetting")
      next
    }
    seurat_list[[basename(sample)]] <- seurat_tmp
  }

  message("All samples processed successfully.")
  return(seurat_list)
}

# Example usage
# gtf_file <- "/cfs/klemming/projects/snic/sllstore2017078/lech/sarek/refs/112/gtf/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf.gz"
# parent_dir <- "/cfs/klemming/projects/snic/sllstore2017078/kaczma-workingdir/RR/scAnalysis/single_cell_gal7b/count_arc"
# seurat_list <- create_seurat_list(parent_dir, gtf_file)
