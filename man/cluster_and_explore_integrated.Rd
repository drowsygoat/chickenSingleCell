% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_and_explore_integrated.R
\name{cluster_and_explore_integrated}
\alias{cluster_and_explore_integrated}
\title{Cluster and Explore an Integrated Seurat Object}
\usage{
cluster_and_explore_integrated(
  seurat_obj,
  dims = 1:30,
  resolution = 1,
  use.assay = "integrated",
  tsne_perplexity = 30,
  output_dir = "cluster_and_explore_integrated_plots",
  save_metadata_csv = FALSE
)
}
\arguments{
\item{seurat_obj}{Integrated Seurat object.}

\item{dims}{PCA dimensions to use (default: 1:30).}

\item{resolution}{Clustering resolution (default: 1).}

\item{use.assay}{Assay to use (default: "integrated").}

\item{tsne_perplexity}{Perplexity for t-SNE (default: 30).}

\item{output_dir}{Output directory for plots (default: "cluster_and_explore_integrated_plots").}

\item{save_metadata_csv}{Logical; whether to save metadata (default: FALSE).}
}
\value{
Updated Seurat object with dimensional reductions and clusters.
}
\description{
Performs dimensionality reduction, clustering, and visualizations on an integrated Seurat object.
Generates cumulative cluster plots as well as individual per-sample highlight plots on one page.
}
