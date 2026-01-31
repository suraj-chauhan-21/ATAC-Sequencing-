# ATAC-Sequencing-
ATAC-seq: “Who left the DNA doors open?” 
The setting  Imagine your genome as a huge apartment building. 
Each room = a gene  
Locked doors = gene OFF 
Open doors = gene ON or ready to be ON  
The problem: we cannot directly see which doors are open.  
So we use a special molecular spy.
The key player: Tn5 transposase 
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(tidyverse)

# 1. Data Loading and Object Creation
counts <- Read10X_h5('data/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5')

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "data/atac_v1_pbmc_10k_fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

metadata <- read.csv('data/atac_v1_pbmc_10k_singlecell.csv', header = T, row.names = 1)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  meta.data = metadata,
  assay = 'ATAC'
)

# 2. Gene Annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
Annotation(pbmc) <- annotations

# 3. Computing Quality Control Metrics
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

# 4. Filtering
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC > 3000 &
    nCount_ATAC < 30000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)

# 5. Normalization and Dimensionality Reduction
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

# 6. Clustering and UMAP Visualization
# Starting from Dim 2 to exclude sequencing depth technical variation
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, algorithm = 3)

# 7. Plotting
DimPlot(object = pbmc, label = TRUE) + NoLegend()
