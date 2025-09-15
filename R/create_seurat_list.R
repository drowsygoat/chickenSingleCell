#' Create Seurat Objects from 10x Multiome (Cell Ranger ARC) Outputs
#'
#' Iterates over sample subdirectories in a parent folder containing standard
#' 10x Genomics Cell Ranger ARC outputs and returns a named list of Seurat
#' objects (one per sample). Each object includes an RNA assay and, optionally,
#' an ATAC assay plus per-barcode metadata when available.
#'
#' @details
#' **Expected directory structure per sample (typical Cell Ranger ARC output):**
#' - Filtered 10x HDF5 matrix: `outs/filtered_feature_bc_matrix.h5`
#' - ATAC fragments: `outs/atac_fragments.tsv.gz` (and index `*.tbi`) *(required if `include_ATAC = TRUE`)*
#' - Per-barcode metrics (if available): e.g., `outs/per_barcode_metrics.csv` or similar
#'
#' For each sample directory found in `parent_dir`, this function calls
#' `import_to_seurat()` to build a multi-assay Seurat object using the provided GTF
#' for gene annotations. Samples yielding zero cells after internal subsetting are
#' skipped with an informative message.
#'
#' @param parent_dir Character path to the parent directory containing one or more
#'   **sample subdirectories** with Cell Ranger ARC outputs.
#' @param gtf_file Character path to a GTF file with gene annotations used to map
#'   features for the RNA assay (and any gene-dependent steps).
#' @param dry_run Logical; if `TRUE`, perform discovery/validation and log what would
#'   be done without creating Seurat objects. Defaults to `FALSE`.
#' @param include_ATAC Logical; if `TRUE`, add an ATAC assay using the fragments file
#'   when present. Defaults to `TRUE`.
#' @param include_metadata Logical; if `TRUE`, attempt to join per-barcode metrics
#'   (e.g., from `per_barcode_metrics.csv`) into object metadata. Defaults to `FALSE`.
#'
#' @return A **named list** of Seurat objects, one per sample directory. List names
#'   are the basenames of the sample directories. Samples with no cells after
#'   subsetting are omitted.
#'
#' @section Notes:
#' - This function assumes Cell Ranger ARC-style layout and filenames; minor
#'   variations may be supported by `import_to_seurat()`.
#' - If `include_ATAC = TRUE` but fragments are missing, the function will proceed
#'   without an ATAC assay for that sample (behavior delegated to `import_to_seurat()`).
#' - Large datasets can be memory-intensive; consider running with `dry_run = TRUE`
#'   first to validate paths and inputs.
#'
#' @seealso
#' \code{\link{import_to_seurat}} for per-sample import logic; packages \pkg{Seurat}
#' and \pkg{Signac} for downstream analysis of RNA/ATAC assays.
#'
#' @examples
#' \dontrun{
#' # Validate inputs without heavy loading
#' create_seurat_list(
#'   parent_dir = "/path/to/parent_dir",
#'   gtf_file    = "/path/to/annotations.gtf",
#'   dry_run = TRUE
#' )
#'
#' # Create RNA+ATAC objects and attach per-barcode metadata
#' seurat_list <- create_seurat_list(
#'   parent_dir = "/path/to/parent_dir",
#'   gtf_file    = "/path/to/annotations.gtf",
#'   dry_run = FALSE,
#'   include_ATAC = TRUE,
#'   include_metadata = TRUE
#' )
#' }
#'
#' @export

create_seurat_list <- function(parent_dir, gtf_file, dry_run = FALSE, include_ATAC = TRUE, include_metadata = FALSE) {

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
    seurat_tmp <- import_to_seurat(sample, gtf_file, dry_run = dry_run, include_ATAC = include_ATAC, include_metadata = include_metadata)
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
