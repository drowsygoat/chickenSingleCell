#!/usr/bin/env Rscript

require(ArchR)
require(stringr)
require(BSgenome.Gg.Ensembl.GRCg7b)
require(Seurat)
require(scran)
require(pheatmap)
require(scran)
require(magick)
require(hexbin)
require(Cairo)

addArchRGenome("hg19")
addArchRDebugging(debug = T)
addArchRLocking(locking = F)
addArchRThreads(threads = 16)

num_cores <- detectCores()
cat("Available CPU cores:", num_cores, "\n")

geneAnnotation_gal7 <- readRDS("/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/archr_gal7/geneAnnotation_gal7.rds")

ArrowFiles <- list.files(path = ".", pattern = "arrow$", recursive = TRUE, full.names = TRUE)

samples <- str_extract(ArrowFiles, "ID.{2}")

samples <- gsub("\\/","", samples)

add_leading_zero <- function(ids) {
  gsub("ID(\\d)$", "ID0\\1", ids)
}

samples <- sapply(samples, add_leading_zero, simplify = T, USE.NAMES = F)

file.exists(ArrowFiles)

geneAnnotation_gal7_no_M <- geneAnnotation_gal7
geneAnnotation_gal7_no_M[[3]][[2]] <- geneAnnotation_gal7_no_M[[3]][[2]][1:41]

atac_with_scores_ArchProj <- ArchRProject(
    ArrowFiles = ArrowFiles,
    showLogo = FALSE,
    genomeAnnotation = geneAnnotation_gal7_no_M[[3]],
    geneAnnotation = geneAnnotation_gal7[[1]],
    outputDirectory = "atac_with_scores_ArchProj",
    copyArrows = TRUE
)

atac_with_scores_ArchProj <- loadArchRProject("atac_with_scores_ArchProj")

save.image(file = "atac_with_scores_ArchProj.RData")

atac_with_scores_ArchProj <- subsetCells(
  ArchRProj = atac_with_scores_ArchProj,
  cellNames = atac_with_scores_ArchProj$cellNames[atac_with_scores_ArchProj$TSSEnrichment >= 7 & atac_with_scores_ArchProj$ReadsInTSS >= 1500]
)

atac_with_scores_ArchProj <- addIterativeLSI(
    ArchRProj = atac_with_scores_ArchProj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI_TM",
    iterations = 3,  # More iterations to refine cell identities
    clusterParams = list(
      resolution = c(0.5),  # Better balance for clustering (https://github.com/GreenleafLab/ArchR/discussions/946)
      sampleCells = NULL,  # More cells sampled per iteration
      maxClusters = NULL,  # Let ArchR determine optimal clusters
      n.start = 10
    ),
    firstSelection = "top",
    depthCol = "nFrags",
    varFeatures = 40000,  # Increased to capture diverse cell types
    dimsToUse = 1:40,  # More dimensions for better separation
    LSIMethod = 2,
    scaleDims = TRUE,
    corCutOff = 0.8,  # Higher cutoff to remove noise
    binarize = TRUE,
    outlierQuantiles = c(0.01, 0.99),  # Capture more outliers
    filterBias = TRUE,
    sampleCellsPre = 50000,  # Larger subsampling for feature selection
    projectCellsPre = FALSE,
    sampleCellsFinal = 10000,  # Use more cells for final projection
    selectionMethod = "var",
    scaleTo = 10000,
    totalFeatures = 750000,  # Increased for brain dataset complexity
    filterQuantile = 0.995,
    excludeChr = c("chrM"),  # or Exclude mitochondrial reads
    keep0lsi = FALSE,
    saveIterations = TRUE,
    UMAPParams = list(
      n_neighbors = 50,  # Higher for large dataset
      min_dist = 0.3,  # Helps separate rare cell populations
      metric = "cosine",
      verbose = TRUE,
      fast_sgd = TRUE
    ),
    nPlot = 10000,  # More cells plotted for visualization
    outDir = getOutputDirectory(atac_with_scores_ArchProj),
    threads = getArchRThreads(),  # Use all available threads
    seed = 1,
    verbose = TRUE,
    force = TRUE,
    logFile = createLogFile("addIterativeLSI_TM")
)

saveArchRProject(atac_with_scores_ArchProj)
save.image(file = "atac_with_scores_ArchProj.RData")

  atac_with_scores_ArchProj <- addHarmony(
  ArchRProj = atac_with_scores_ArchProj,
  reducedDims = "IterativeLSI_TM",
  dimsToUse = NULL,
  scaleDims = NULL,
  corCutOff = 0.75,
  name = "Harmony_TM",
  groupBy = "Sample",
  verbose = TRUE,
  force = TRUE)

save.image(file = "atac_with_scores_ArchProj.RData")
quit(save = "no", status = 0)   

atac_with_scores_ArchProj <- addClusters(
  input = atac_with_scores_ArchProj,
  reducedDims = "IterativeLSI_TM",  # Use Harmony if batch-corrected, or "IterativeLSI"
  name = "SeuratClusters_TM",
  sampleCells = NULL,       # Larger sample size for better cluster resolution
  method = "Seurat",
  dimsToUse = NULL,         # More dimensions for brain cell diversity
  scaleDims = TRUE,         # Standardize dimensions
  corCutOff = 0.8,          # Higher cutoff to reduce noise
  knnAssign = 30,           # Higher KNN for large datasets
  nOutlier = 10,            # Allow more outliers in complex tissue
  maxClusters = 100,        # Increase to capture diverse brain cell types
  testBias = TRUE,          # Ensure no bias by fragment count
  filterBias = TRUE,        # Remove cells with high sequencing bias
  biasClusters = 0.01,
  biasCol = "nFrags",
  biasQuantiles = c(0.05, 0.95),
  biasEnrich = 10,
  nPerm = 500,
    force = TRUE,
  verbose = TRUE,
  # Seurat args
  initial.membership = NULL,
  node.sizes = NULL,
  resolution = 1.2,
  algorithm = 1,
  n.start = 20,
  n.iter = 20,
  random.seed = 0,
  group.singletons = TRUE,
  temp.file.location = NULL,
  edge.file.name = NULL,
  #######
  logFile = createLogFile("addClusters")
)

atac_with_scores_ArchProj <- addClusters(
  input = atac_with_scores_ArchProj,
  reducedDims = "Harmony_TM",  # Use Harmony if batch-corrected, or "IterativeLSI"
  name = "SeuratClustersHarmony_TM",
  sampleCells = NULL,     # Larger sample size for better cluster resolution
  method = "Seurat",
  dimsToUse = NULL,         # More dimensions for brain cell diversity
  scaleDims = TRUE,         # Standardize dimensions
  corCutOff = 0.8,          # Higher cutoff to reduce noise
  knnAssign = 30,           # Higher KNN for large datasets
  nOutlier = 10,            # Allow more outliers in complex tissue
  maxClusters = 100,        # Increase to capture diverse brain cell types
  testBias = TRUE,          # Ensure no bias by fragment count
  filterBias = TRUE,        # Remove cells with high sequencing bias
  biasClusters = 0.01,
  biasCol = "nFrags",
  biasQuantiles = c(0.05, 0.95),
  biasEnrich = 10,
  nPerm = 500,
  verbose = TRUE,
  # Seurat args
  initial.membership = NULL,
  node.sizes = NULL,
  resolution = 1.2,
  algorithm = 1,
  n.start = 20,
  n.iter = 20,
  random.seed = 0,
  group.singletons = TRUE,
  temp.file.location = NULL,
  edge.file.name = NULL,
    force = TRUE,
  #######
  logFile = createLogFile("addClusters")
)

saveArchRProject(atac_with_scores_ArchProj)
save.image(file = "atac_with_scores_ArchProj.RData")

atac_with_scores_ArchProj <- addUMAP(
    ArchRProj = atac_with_scores_ArchProj, 
    reducedDims = "IterativeLSI_TM", 
    name = "UMAP_LSI_TM", 
    nNeighbors = 40,
    minDist = 0.4,
    metric = "cosine",
    dimsToUse = NULL,
    scaleDims = NULL,
    corCutOff = 0.75,
    sampleCells = NULL,
    outlierQuantile = 0.9,
    saveModel = TRUE,
    verbose = TRUE,
    seed = 1,
  force = TRUE,
    threads = 1
)

atac_with_scores_ArchProj <- addUMAP(
    ArchRProj = atac_with_scores_ArchProj, 
    reducedDims = "Harmony_TM", 
    name = "UMAP_Harmony_TM", 
    nNeighbors = 40,
    minDist = 0.4,
    metric = "cosine",
    dimsToUse = NULL,
    scaleDims = NULL,
    corCutOff = 0.75,
    sampleCells = NULL,
    outlierQuantile = 0.9,
    saveModel = TRUE,
    verbose = TRUE,
    seed = 1,
  force = TRUE,
    threads = 1
)

saveArchRProject(atac_with_scores_ArchProj)
save.image(file = "atac_with_scores_ArchProj.RData")

atac_with_scores_ArchProj <- addTSNE(
  ArchRProj = atac_with_scores_ArchProj,
  reducedDims = "IterativeLSI_TM",
  method = "RTSNE",
  name = "TSNE_LSI_TM",
  perplexity = 50,
  maxIterations = 1000,
  learningRate = 200,
  dimsToUse = NULL,
  scaleDims = NULL,
  corCutOff = 0.75,
  saveModel = FALSE,
  verbose = TRUE,
  seed = 1,
  force = TRUE,
  threads = max(floor(getArchRThreads()/2), 1)
)

atac_with_scores_ArchProj <- addTSNE(
  ArchRProj = atac_with_scores_ArchProj,
  reducedDims = "Harmony_TM",
  method = "RTSNE",
  name = "TSNE_Harmony_TM",
  perplexity = 50,
  maxIterations = 1000,
  learningRate = 200,
  dimsToUse = NULL,
  scaleDims = NULL,
  corCutOff = 0.75,
  saveModel = FALSE,
  verbose = TRUE,
  seed = 1,
  force = TRUE,
  threads = max(floor(getArchRThreads()/2), 1)
)

saveArchRProject(atac_with_scores_ArchProj)
save.image(file = "pilot_scores.RData")

load(file = "pilot_scores.RData")

p1 <- plotEmbedding(ArchRProj = atac_with_scores_ArchProj, colorBy = "cellColData", name = "Sample", embedding = "UMAP_LSI_TM")

p2 <- plotEmbedding(ArchRProj = atac_with_scores_ArchProj, colorBy = "cellColData", name = "SeuratClusters_TM", embedding = "UMAP_LSI_TM")

p3 <- plotEmbedding(ArchRProj = atac_with_scores_ArchProj, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Harmony_TM")

p4 <- plotEmbedding(ArchRProj = atac_with_scores_ArchProj, colorBy = "cellColData", name = "SeuratClustersHarmony_TM", embedding = "UMAP_Harmony_TM")

p5 <- plotEmbedding(ArchRProj = atac_with_scores_ArchProj, colorBy = "cellColData", name = "Sample", embedding = "TSNE_LSI_TM")

p6 <- plotEmbedding(ArchRProj = atac_with_scores_ArchProj, colorBy = "cellColData", name = "SeuratClusters_TM", embedding = "TSNE_LSI_TM")

p7 <- plotEmbedding(ArchRProj = atac_with_scores_ArchProj, colorBy = "cellColData", name = "Sample", embedding = "TSNE_Harmony_TM")

p8 <- plotEmbedding(ArchRProj = atac_with_scores_ArchProj, colorBy = "cellColData", name = "SeuratClustersHarmony_TM", embedding = "TSNE_Harmony_TM")

saveArchRProject(atac_with_scores_ArchProj)
save.image(file = "atac_with_scores_ArchProj.RData")

plotPDF(p1,p2,p3,p4,p5,p6,p7,p8, name = "clusters.pdf", ArchRProj = atac_with_scores_ArchProj, addDOC = FALSE, width = 10, height = 10)

library(cowplot)

# Define the embeddings to iterate through
embeddings <- c("UMAP_LSI_TM", "UMAP_Harmony_TM", "TSNE_LSI_TM", "TSNE_Harmony_TM")

# Safe plotting function with error handling
safe_plotEmbedding <- function(proj, sample, embedding) {
  cluster <- if (grepl("LSI", embedding)) {
    "SeuratClusters_TM"
  } else {
    "SeuratClustersHarmony_TM"
  }

  tryCatch({
    plotEmbedding(
      ArchRProj = proj[which(getCellColData(proj, "Sample", drop = TRUE) == sample)],
      embedding = embedding,
      colorBy = "cellColData",
      name = cluster,
      size = 1.5,
      sampleCells = NULL,
      highlightCells = NULL,
      baseSize = 8,
      labelSize = 2,
      threads = 10,
      plotAs = "points"
    ) + guides(color = FALSE, fill = FALSE)
  }, error = function(e) {
    message(paste("Error for sample:", sample, "embedding:", embedding, "cluster:", cluster))
    ggplot() + theme_void() + ggtitle(paste("Error:", sample))
  })
}

# Loop through selected embeddings
for (embedding in embeddings) {
  plotlist <- list()

  # Loop through each sample
  for (sample in unique(atac_with_scores_ArchProj@cellColData$Sample)) {
    plotlist[[sample]] <- safe_plotEmbedding(
      proj = atac_with_scores_ArchProj,
      sample = sample,
      embedding = embedding
    )
  }

  # Customize plots
  plotlist2 <- lapply(names(plotlist), function(sample) {
    plotlist[[sample]] +
      guides(color = FALSE, fill = FALSE) +
      ggtitle(paste(sample, "-", embedding)) +
      theme_ArchR(baseSize = 15) +
      theme(
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  })

  # Create grid plot
  cowp <- cowplot::plot_grid(
    plotlist = plotlist2,
    ncol = 5,
    align = 'hv',
    hjust = -1,
    vjust = -1
  )

  # Save the plot
  ggsave(
    file = paste0("clusters_per_sample_", embedding, "_rush.pdf"),
    plot = cowp,
    width = 30,
    height = 30,
    dpi = 800
  )
}

saveArchRProject(atac_with_scores_ArchProj)
save.image(file = "atac_with_scores_ArchProj.RData")

atac_with_scores_ArchProj <- addGroupCoverages(
  ArchRProj = atac_with_scores_ArchProj,
  groupBy = "SeuratClustersHarmony_TM",
  useLabels = TRUE,
  minCells = 1,
  maxCells = 2000,
  maxFragments = 5000 * 10^6, # * 200
  minReplicates = 2,
  maxReplicates = 34,
  sampleRatio = 0.8,
  kmerLength = 6,
  threads = getArchRThreads(),fff
  returnGroups = FALSE,
  parallelParam = NULL,
  force = T,
  verbose = TRUE,
  logFile = createLogFile("addGroupCoverages")
)

saveArchRProject(atac_with_scores_ArchProj, load = FALSE)
save.image(file = "atac_with_scores_ArchProj.RData")

seMarkerFeatures <- getMarkerFeatures(
    ArchRProj = atac_with_scores_ArchProj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "SeuratClustersHarmony_TM",
    maxCells = 500 * 10, #!!!
    scaleTo = 10^4,
    k = 100,
    threads = 2,
    bufferRatio = 0.8,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

saveArchRProject(atac_with_scores_ArchProj, load = FALSE)
save.image(file = "atac_with_scores_ArchProj.RData")

 atac_with_scores_ArchProj <- addReproduciblePeakSet(
  ArchRProj = atac_with_scores_ArchProj,
  groupBy = "SeuratClustersHarmony_TM",
  peakMethod = "Macs2",
  reproducibility = "1", # changed !!!!!!!!
  peaksPerCell = 500 * 100, # x100
  maxPeaks = 150000 * 100, # x100
  minCells = 25, # kept, but maybe lower??
  excludeChr = c("chrM"),
  pathToMacs2 = "macs3",
  genomeSize = 1065000000, # this is from UCSC
  shift = -75,
  extsize = 150,
  method = "q",
  cutOff = 0.1, # change?
  additionalParams = "--nomodel --nolambda",
  extendSummits = 250,
  promoterRegion = c(2000, 100),
  genomeAnnotation = geneAnnotation_gal7_no_M[[3]],
  geneAnnotation = geneAnnotation_gal7[[1]],
  plot = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE,
  verbose = TRUE,
  logFile = createLogFile("addReproduciblePeakSet")
)

saveArchRProject(atac_with_scores_ArchProj, load = FALSE)
save.image(file = "atac_with_scores_ArchProj.RData")

atac_with_scores_ArchProj <- addPeakMatrix(atac_with_scores_ArchProj)
saveArchRProject(atac_with_scores_ArchProj, load = FALSE)


# Define cluster column and extract cluster IDs
# cluster_col <- "SeuratClustersHarmony_TM"
# cluster_assignments <- getCellColData(atac_with_scores_ArchProj)[[cluster_col]]
# clusters <- unique(cluster_assignments)

atac_with_scores_ArchProj <- addCellColData(ArchRProj = atac_with_scores_ArchProj, data = paste0(atac_with_scores_ArchProj@cellColData$SeuratClustersHarmony_TM, "_x_", atac_with_scores_ArchProj@cellColData$Sample), name = "SeuratClustersHarmony_TM_x_Sample", cells = getCellNames(atac_with_scores_ArchProj), force = TRUE)

atac_with_scores_ArchProj_scale_10k_full <- getGroupSE(ArchRProj = atac_with_scores_ArchProj, useMatrix = "PeakMatrix", scaleTo = 10000, groupBy = "SeuratClustersHarmony_TM_x_Sample")

# Save the SummarizedExperiment object with clusters
saveRDS(atac_with_scores_ArchProj_scale_10k_full, "atac_with_scores_ArchProj_scale_10k_full.rds")