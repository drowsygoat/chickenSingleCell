% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_ArchR_meta_to_Seurat.R
\name{add_ArchR_meta_to_Seurat}
\alias{add_ArchR_meta_to_Seurat}
\title{Merge ArchR Metadata into a Seurat Object}
\usage{
add_ArchR_meta_to_Seurat(
  seurat_obj,
  archr_df,
  mode = c("intersect", "add"),
  absent_value = "absent",
  archr_value_cols = NULL,
  dry_run = FALSE
)
}
\arguments{
\item{seurat_obj}{A Seurat object.}

\item{archr_df}{A data.frame with ArchR cell names as rownames and metadata columns.}

\item{mode}{"intersect" keeps only matched cells, "add" preserves all Seurat cells.}

\item{absent_value}{Value assigned to unmatched cells in "add" mode.}

\item{archr_value_cols}{Optional. Character vector of column names in \code{archr_df} to merge. If NULL, merges all columns.}

\item{dry_run}{If TRUE, only prints merge stats and doesn't modify the Seurat object.}
}
\value{
Modified Seurat object (or NULL if no matches found / dry_run = TRUE)
}
\description{
Merges one or more metadata columns from an ArchR-exported dataframe (e.g., cellColData)
into a Seurat object. Matches cells by sample and barcode, handling different ID separators ("_" vs "#").
}
