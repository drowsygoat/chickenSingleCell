gtf_file <- "/cfs/klemming/projects/snic/sllstore2017078/lech/sarek/refs/112/gtf/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf.gz"
parent_dir <- "/cfs/klemming/projects/snic/sllstore2017078/kaczma-workingdir/RR/scAnalysis/single_cell_gal7b/count_arc"

library(Seurat)
library(tidyverse)
library(Signac)
library(GenomicRanges)
library(harmony)  # Optional batch correction
library(GenomeInfoDb)  # For chromosome annotation
library(BSgenome.Gg.Ensembl.GRCg7b)  # Chicken genome

devtools::load_all()

### --------------------

# seurat_list <- create_seurat_list(parent_dir, gtf_file, dry_run = F)
# save.image(file = "workshop_trimmed.RData")
# quit(save = "no", status = 0)

### --------------------
load("workshop_trimmed.RData")

# for (i in seq_along(seurat_list)) {
#   obj <- seurat_list[[i]]
#   cat("\nChecking object", i, "\n")
  
#   # Show some example cell names and ident names
#   cat("Head of colnames:\n")
#   print(head(colnames(obj)))
  
#   cat("Head of names(Idents):\n")
#   print(head(names(Idents(obj))))
  
#   # Check if they are identical
#   match <- identical(colnames(obj), names(Idents(obj)))
#   cat("Do colnames and Idents names match?:", match, "\n")
  
#   # Check for missing values
#   if (any(is.na(names(Idents(obj))))) {
#     cat("Warning: NA values found in Idents names!\n")
#   }
# }

import_test <- import_to_seurat(sample_dir = "/cfs/klemming/projects/snic/sllstore2017078/kaczma-workingdir/RR/scAnalysis/single_cell_gal7b/count_arc/ID5/", gtf_file = gtf_file)

import_test <- normalize_with_sct(import_test)

label_high_mito <- function(import_test, threshold = 10) {
  if (!"percent.mt" %in% colnames(import_test@meta.data)) {
    import_test[["percent.mt"]] <- PercentageFeatureSet(import_test, pattern = "^ND1|^ND3|^ND4|^ND5|^ND6|^ATP6$|^ATP8$|^CYTB|^COII|^COX3", assay = "RNA")
  }
  import_test$high_mito <- ifelse(import_test$percent.mt > threshold, "High", "Normal")
  return(import_test)
}

import_test <- label_high_mito(import_test)

import_test <- cluster_and_explore(import_test)






import_test_2[["percent.mt"]] <- PercentageFeatureSet(import_test_2, pattern = "^MT-|^ND1|^ND3|^ND4|^ND5|^ND6|^ATP6$|^ATP8$|^CYTB|^COII|^COX3")

import_test_2 <- NucleosomeSignal(import_test_2, assay = "ATAC")
import_test_2 <- TSSEnrichment(import_test_2, assay = "ATAC")


seurat_merged[["percent.mt"]] <- PercentageFeatureSet(seurat_merged, pattern = "^MT-|^ND1|^ND3|^ND4|^ND5|^ND6|^ATP6$|^ATP8$|^CYTB|^COII|^COX3")

seurat_merged <- seurat_list[[4]]  # Start with the first object


import_test@

for (i in 2:length(seurat_list[4:5])) {
  seurat_merged <- merge(seurat_merged, seurat_list[[i]])
}

pdf()
VlnPlot(seurat_merged, features = c("percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()

  # Load RNA and ATAC counts
  message("Loading RNA and ATAC counts...")
  counts <- Read10X_h5("/cfs/klemming/projects/snic/sllstore2017078/kaczma-workingdir/RR/scAnalysis/single_cell_gal7b/count_arc/ID4/outs/filtered_feature_bc_matrix.h5")

  if (!"Gene Expression" %in% names(counts) || !"Peaks" %in% names(counts)) {
    stop("Expected 'Gene Expression' and 'Peaks' data in the H5 file.")
  }

  rna_counts <- counts$`Gene Expression`
  atac_counts <- counts$Peaks

  # Create Seurat object
  seurat_obj <- NULL

  # Include RNA data
  if (include_RNA) {
    message("Creating Seurat object for RNA assay...")
    seurat_obj <- CreateSeuratObject(counts = rna_counts, assay = "RNA")
  }

  # Import GTF as GRanges
  if (!file.exists(gtf_file)) {
    stop("GTF file not found at: ", gtf_file)
  }



# seurat_merged <- Reduce(function(x, y) merge(x, y), seurat_list)

seurat_merged <- seurat_list[[1]]  # Start with the first object

for (i in 2:length(seurat_list)) {
  seurat_merged <- merge(seurat_merged, seurat_list[[i]])
}

save.image(file = "workshop_trimmed.RData")
quit(save = "no", status = 0)

seurat_subset <- subset_random_cells(seurat_merged, fraction = 0.01)

Harmony()

# ✅ Step 1: Normalize RNA with SCTransform (Better for Brain Data)
seurat_merged <- SCTransform(seurat_merged, assay = "RNA", verbose = FALSE)
seurat_merged <- FindVariableFeatures(seurat_merged, selection.method = "vst", nfeatures = 3000)
seurat_merged <- RunPCA(seurat_merged, npcs = 50)

# ✅ (Optional) Batch Correction
seurat_merged <- RunHarmony(seurat_merged, group.by.vars = "batch", reduction = "pca")

# ✅ Step 2: Process ATAC Data (TF-IDF & Dimensionality Reduction)
seurat_merged <- RunTFIDF(seurat_merged, assay = "ATAC", method = 2)
seurat_merged <- FindTopFeatures(seurat_merged, assay = "ATAC", min.cutoff = 20)
seurat_merged <- RunSVD(seurat_merged, reduction.name = "lsi", n = 50)

# ✅ Step 3: Ensure Proper Chicken Genome Format
seqlevelsStyle(seurat_merged) <- "UCSC"  # Ensures ATAC peaks match genome

# ✅ Step 4: Find Multiome Anchors for RNA-ATAC Integration
transfer_anchors <- FindTransferAnchors(
  reference = seurat_merged, 
  query = seurat_merged, 
  reduction = "cca", 
  features = VariableFeatures(seurat_merged), 
  dims = 1:50
)

# ✅ Step 5: Transfer RNA Labels to ATAC Cells
seurat_merged <- MapQuery(
  anchorset = transfer_anchors,
  reference = seurat_merged,
  query = seurat_merged,
  refdata = list(celltype = "predicted.celltype"),
  reduction.model = "umap"
)

# ✅ Step 6: Multi-modal Clustering
seurat_merged <- FindMultiModalNeighbors(
  seurat_merged, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 1:50)
)

# ✅ Step 7: UMAP and Clustering
seurat_merged <- RunUMAP(seurat_merged, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
seurat_merged <- FindClusters(seurat_merged, graph.name = "wsnn", resolution = 0.8)

# ✅ Step 8: Visualization
DimPlot(seurat_merged, reduction = "wnn.umap", group.by = "seurat_clusters")
DimPlot(seurat_merged, reduction = "wnn.umap", group.by = "predicted.celltype")

# ------------------------------------------------------------
# ✅ Step 9: Peak-to-Gene Link Analysis (Chicken Genome Adjusted)
# ------------------------------------------------------------

# Convert chromosome names if needed
seqlevelsStyle(seurat_merged) <- "UCSC"

# ✅ Find Peak-to-Gene Links
seurat_merged <- LinkPeaks(
  object = seurat_merged,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  genes.use = rownames(seurat_merged),  # Use all expressed genes
  method = "pearson"
)

# ✅ Visualize a Specific Gene & Its Linked Peaks (Example: Neurogenesis marker)
gene_of_interest <- "NEUROD1"  # Chicken neurodevelopment marker
CoveragePlot(
  object = seurat_merged,
  region = gene_of_interest,
  features = gene_of_interest,
  assay = "ATAC",
  expression.assay = "RNA",
  extend.upstream = 50000,
  extend.downstream = 50000
)

# ✅ Explore Peak-Gene Links
head(seurat_merged@assays$ATAC@links)  # Show linked peaks

# ✅ Plot peak-to-gene correlations
LinkedFeatures(seurat_merged, features = gene_of_interest, assay = "ATAC")



save.image(file = "workshop.RData")



grep()

for (i in names(seurat_list)) {
  seurat_list[[i]][["percent.mt"]] <- PercentageFeatureSet(
    seurat_list[[i]], pattern = "^MT-|^ND1|^ND3|^ND4|^ND5|^ND6|^ATP6$|^ATP8$|^CYTB|^COII|^COX3"
  )
  
  # Plot QC metrics
  VlnPlot(seurat_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
}

### --------------------