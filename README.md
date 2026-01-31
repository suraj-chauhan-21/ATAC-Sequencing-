******# ATAC-Sequencing-
ATAC-seq: “Who left the DNA doors open?” 
The setting  Imagine your genome as a huge apartment building. 
Each room = a gene  
Locked doors = gene OFF 
Open doors = gene ON or ready to be ON  
The problem: we cannot directly see which doors are open.  
So we use a special molecular spy.
The key player: Tn5 transposase 

# ==============================================================================
# Single-Cell ATAC-Seq Pipeline (Signac/Seurat)
# ==============================================================================

# 1. INITIALIZATION & DATA LOADING ---------------------------------------------
# We need Signac for ATAC-specific math and EnsDb for the gene "Address Book."

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75) 
library(tidyverse)

# Load the count matrix (The summary table of open doors)
counts <- Read10X_h5('data/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5')

# Create the Chromatin Assay 
# 'fragments' is the most important file; it's the raw list of every "tag" left by Tn5.
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "data/atac_v1_pbmc_10k_fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

# Build the Seurat Object (The "Building Manager" that holds all our data)
metadata <- read.csv(file = 'data/atac_v1_pbmc_10k_singlecell.csv', header = T, row.names = 1)
pbmc <- CreateSeuratObject(counts = chrom_assay, meta.data = metadata, assay = 'ATAC')


# 2. ANNOTATION ----------------------------------------------------------------
# Why? Raw data gives coordinates (e.g. chr1:100-200). 
# Annotation tells us if that coordinate is a bedroom, a kitchen, or a hallway (a gene).

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations)) # Match UCSC naming
Annotation(pbmc) <- annotations


# 3. QUALITY CONTROL (The "Bouncers") ------------------------------------------
# We filter out low-quality cells to ensure our "map" is accurate.

# Nucleosome Signal: Checks if DNA is breaking at the right intervals (around spools).
pbmc <- NucleosomeSignal(pbmc)

# TSS Enrichment: Checks if there is a high "pile-up" of tags at gene front doors.
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# Blacklist Ratio: Ignores "trash" regions of the genome that always show up.
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

# Filtering criteria:
pbmc <- subset(x = pbmc,
                subset = nCount_ATAC > 3000 & nCount_ATAC < 30000 &
                  pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 &
                  nucleosome_signal < 4 & TSS.enrichment > 3)


# 4. NORMALIZATION & DIMENSIONAL REDUCTION ------------------------------------
# Why TF-IDF? DNA accessibility is binary (open/closed). 
# TF-IDF (Term Frequency-Inverse Document Frequency) highlights the "rare" 
# open doors that define a specific cell type.

pbmc <- RunTFIDF(pbmc)              # Normalize
pbmc <- FindTopFeatures(pbmc)       # Pick the most informative "doors"
pbmc <- RunSVD(pbmc)                # Linear reduction (LSI)

# 5. CLUSTERING & VISUALIZATION ------------------------------------------------
# We group cells with similar "open doors" into clusters.
# We start at Dim 2 because Dim 1 usually just represents total "sequencing depth."

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, algorithm = 3)

# Plot the neighborhood map
DimPlot(object = pbmc, label = TRUE) + NoLegend()
