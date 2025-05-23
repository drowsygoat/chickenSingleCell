% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infer_sample_id.R
\name{infer_sample_id}
\alias{infer_sample_id}
\title{Infer Sample ID from a Seurat Object}
\usage{
infer_sample_id(seurat_obj)
}
\arguments{
\item{seurat_obj}{A Seurat object.}
}
\value{
A character vector of sample IDs, one per cell.
}
\description{
Attempts to infer sample identity for each cell in a Seurat object.
First checks metadata columns (\code{orig.ident}, \verb{sample.*}, \verb{Sample.*}, \verb{SAMPLE.*}).
If not found, attempts to parse the sample name from cell name prefixes using \verb{_} or \verb{#} as separators.
}
