gtf_file <- "/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/sarek/refs/112/gtf/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf.gz"

parent_dir <- "/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/single_cell_gal7b/count_arc"

devtools::load_all("/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/seurat_gal7/chickenSingleCell")

library(Seurat)
library(tidyverse)
library(Signac)
library(GenomicRanges)
library(harmony)  # Optional batch correction
library(GenomeInfoDb)
library(BSgenome.Gg.Ensembl.GRCg7b)  # Chicken genome
library(ggbeeswarm)

options(future.globals.maxSize = 90 * 1024^3) 

### --------------------

seurat_list <- create_seurat_list(parent_dir, gtf_file, dry_run = F,include_ATAC = FALSE, include_metadata = FALSE)

seurat_list <- lapply(seurat_list, label_high_mito)

seurat_list <- lapply(seurat_list, normalize_with_sct, vars.to.regress = c("percent.mt_RNA"))

clusters_atac_file <- "/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/archr_gal7/gal7_hypo_with_gene_scores/clusters_atac_harmony_TM.rds"

atac_clusters <- readRDS(file = clusters_atac_file)

seurat_list <- lapply(seurat_list, add_ArchR_meta_to_Seurat, archr_df = atac_clusters, atac_clusters, mode = "intersect", dry_run = FALSE)

se <- AverageExpressionPerClusterPerSample2(seurat_list, cluster_col = "SeuratClustersHarmony_TM", assay = "RNA", layer = "counts", return.se = TRUE, se.save.path = "results/pseudobulk_expression_test/log_norm_atac_clusters_rush.rds")

