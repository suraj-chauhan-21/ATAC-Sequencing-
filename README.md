# ATAC-Sequencing-
ATAC-seq: “Who left the DNA doors open?” 
The setting  Imagine your genome as a huge apartment building. 
Each room = a gene  
Locked doors = gene OFF 
Open doors = gene ON or ready to be ON  
The problem: we cannot directly see which doors are open.  
So we use a special molecular spy.
The key player: Tn5 transposase 

# ==============================================================================
# SINGLE-CELL ATAC-SEQ PROCESSING PIPELINE
# Analogy: "The Apartment Building" - Mapping which doors (genes) are open.
# ==============================================================================

# 1. SETUP AND LIBRARIES
# ----------------------
# Concept: Gathering our tools. We need Signac (for ATAC) and Seurat (for single-cell).
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75) # The "Address Book" for human genes
library(tidyverse)

# 2. DATA LOADING
# ---------------
# Concept: Reading the "Locksmith's Notes" (Fragment file).
# Why? The fragment file tells us exactly where the Tn5 enzyme cut the DNA.
counts <- Read10X_h5('data/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5')

# Creating a ChromatinAssay: This tells the computer that these aren't RNA counts,
# but physical regions of DNA that were "open."
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "data/atac_v1_pbmc_10k_fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

# Creating the Seurat Object: This is our "Building Manager." It holds all the info.
metadata <- read.csv(file = 'data/atac_v1_pbmc_10k_singlecell.csv', header = T, row.names = 1)
pbmc <- CreateSeuratObject(counts = chrom_assay, meta.data = metadata, assay = 'ATAC')

# 3. GENE ANNOTATION
# ------------------
# Concept: Adding the Map. 
# Why? Raw data only gives coordinates (e.g., Chromosome 1: 100-200). 
# Annotations tell us "Chromosome 1: 100-200 is the door to the T-cell Receptor gene."
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations)) # Match UCSC naming style
Annotation(pbmc) <- annotations

# 4. QUALITY CONTROL (QC)
# -----------------------
# Concept: Evicting "Bad" Cells.
# Why? We only want healthy cells. 

# A. Nucleosome Signal: DNA is wrapped around histones (spools). 
#    If the DNA is breaking randomly (not at spools), the cell is likely dying.
pbmc <- NucleosomeSignal(pbmc)

# B. TSS Enrichment: The "Transcription Start Site" is the gene's front door.
#    In a good experiment, the locksmith (Tn5) should be very busy at these doors.
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# C. Blacklist Ratio: Some parts of the genome are "sticky" and grab enzymes 
#    randomly. We want to ignore these "trash cans" of the genome.
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

# VISUALIZING QC
# Why? To see where to draw the line. We want high TSS scores and low Nucleosome signals.
VlnPlot(pbmc, features = c('TSS.enrichment', 'nucleosome_signal'), ncol = 2)

# FILTERING
pbmc <- subset(x = pbmc,
                subset = nCount_ATAC > 3000 & nCount_ATAC < 30000 &
                  pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 &
                  nucleosome_signal < 4 & TSS.enrichment > 3)

# 5. NORMALIZATION & MATH (TF-IDF)
# --------------------------------
# Concept: Highlighting the unique features.
# Why? Simple math doesn't work for ATAC because DNA is binary (the door is either 
# open or closed). TF-IDF (Term Frequency-Inverse Document Frequency) lowers the 
# volume on "boring" doors open in every cell and boosts "interesting" doors.
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc) # Dimensionality reduction (LSI)

# 6. CLUSTERING & MAPPING (UMAP)
# ------------------------------
# Concept: Grouping similar apartments into neighborhoods.
# Why? We want to see if all "T-cells" have the same doors open.
# We skip Dim 1 because it often represents "cell size/depth" rather than biology.
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, algorithm = 3)

# FINAL VISUALIZATION
# The UMAP plot: Each dot is a cell. Clusters are different cell types!
DimPlot(object = pbmc, label = TRUE) + NoLegend()
