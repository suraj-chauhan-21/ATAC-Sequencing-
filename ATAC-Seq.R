###############################################################

# 01_scATAC_signac_pipeline.R

# Single-cell ATAC-seq Analysis Pipeline (10x Genomics + Signac)

#

# Author: Suraj Chauhan

# Description:

# Loads 10x scATAC data → QC → LSI → Clustering → Save object

#

# Expected directory structure:

# project/

# ├── 01_scATAC_signac_pipeline.R

# ├── data/

# └── results/   (auto-created)

###############################################################

# -----------------------------#

# 0. Setup Paths

# -----------------------------#

project_dir <- getwd()
data_dir    <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")

dir.create(results_dir, showWarnings = FALSE)

cat("Project directory:", project_dir, "\n")

# -----------------------------#

# 1. Check Required Packages

# -----------------------------#

required_pkgs <- c(
"Signac",
"Seurat",
"EnsDb.Hsapiens.v75",
"tidyverse",
"patchwork"
)

for (pkg in required_pkgs) {
if (!requireNamespace(pkg, quietly = TRUE)) {
stop(paste("Missing package:", pkg,
"\nInstall it before running the pipeline."))
}
}

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(tidyverse)
library(patchwork)

cat("All packages loaded successfully.\n")

# -----------------------------#

# 2. Define Input Files

# -----------------------------#

frag_file   <- file.path(data_dir, "atac_v1_pbmc_10k_fragments.tsv.gz")
matrix_file <- file.path(data_dir, "atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
meta_file   <- file.path(data_dir, "atac_v1_pbmc_10k_singlecell.csv")

if (!file.exists(frag_file)) stop("Fragment file not found.")
if (!file.exists(matrix_file)) stop("Matrix file not found.")
if (!file.exists(meta_file)) stop("Metadata file not found.")

# -----------------------------#

# 3. Load 10x Peak Matrix

# -----------------------------#

cat("Reading 10x peak matrix...\n")
counts <- Read10X_h5(matrix_file)

# -----------------------------#

# 4. Create Chromatin Assay

# -----------------------------#

cat("Creating ChromatinAssay...\n")

chrom_assay <- CreateChromatinAssay(
counts = counts,
sep = c(":", "-"),
fragments = frag_file,
min.cells = 10,
min.features = 200
)

# -----------------------------#

# 5. Load Metadata and Create Object

# -----------------------------#

metadata <- read.csv(meta_file, row.names = 1)

pbmc <- CreateSeuratObject(
counts = chrom_assay,
assay = "ATAC",
meta.data = metadata
)

cat("Seurat object created.\n")

# -----------------------------#

# 6. Add Gene Annotation

# -----------------------------#

cat("Adding gene annotations...\n")

annotations <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v75)

# Convert Ensembl → UCSC naming

seqlevels(annotations) <- paste0("chr", seqlevels(annotations))

Annotation(pbmc) <- annotations

# -----------------------------#

# 7. Compute QC Metrics

# -----------------------------#

cat("Computing QC metrics...\n")

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc, fast = FALSE)

pbmc$pct_reads_in_peaks <-
pbmc$peak_region_fragments / pbmc$passed_filters * 100

pbmc$blacklist_ratio <-
pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

# -----------------------------#

# 8. QC Filtering (Editable Thresholds)

# -----------------------------#

cat("Filtering low-quality cells...\n")

min_fragments <- 3000
max_fragments <- 30000
min_TSS       <- 3
max_nuc       <- 4
max_blacklist <- 0.05
min_peak_pct  <- 15

pbmc <- subset(
pbmc,
subset =
nCount_ATAC > min_fragments &
nCount_ATAC < max_fragments &
pct_reads_in_peaks > min_peak_pct &
blacklist_ratio < max_blacklist &
nucleosome_signal < max_nuc &
TSS.enrichment > min_TSS
)

cat("Remaining cells:", ncol(pbmc), "\n")

# -----------------------------#

# 9. Dimensional Reduction (LSI)

# -----------------------------#

cat("Running TF-IDF normalization...\n")

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
pbmc <- RunSVD(pbmc)

DepthCor(pbmc)

# -----------------------------#

# 10. Clustering + UMAP

# -----------------------------#

cat("Running clustering...\n")

lsi_dims <- 2:30   # Exclude LSI_1 (depth effect)

pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = lsi_dims)
pbmc <- FindNeighbors(pbmc, reduction = "lsi", dims = lsi_dims)
pbmc <- FindClusters(pbmc, algorithm = 3)

# -----------------------------#

# 11. Save Results

# -----------------------------#

cat("Saving processed object...\n")

saveRDS(pbmc, file = file.path(results_dir, "pbmc_atac_processed.rds"))

# Save UMAP figure

png(file.path(results_dir, "UMAP_clusters.png"), width = 1800, height = 1400, res = 200)
print(DimPlot(pbmc, label = TRUE) + NoLegend())
dev.off()

# -----------------------------#

# 12. Save Session Info (Reproducibility)

# -----------------------------#

writeLines(
capture.output(sessionInfo()),
file.path(results_dir, "sessionInfo.txt")
)

cat("Pipeline completed successfully.\n")
###############################################################




# ATAC-Sequencing

ATAC-seq: “Who left the DNA doors open?” 
The setting  Imagine your genome as a huge apartment building. 
Each room = a gene  
Locked doors = gene OFF 
Open doors = gene ON or ready to be ON  
The problem: we cannot directly see which doors are open.  

![image alt](https://github.com/suraj-chauhan-21/ATAC-Sequencing-/blob/8a8a1f6c99040e46ad86cc98fd900c4e2b192d40/Tn5_Transposase_in_ATAC-seq.webp)

So we use a special molecular spy.
The key player: Tn5 transposase 
# script to process single-cell ATAC-Seq data
Vignette: https://stuartlab.org/signac/articles/pbmc_vignette
setwd("~/Desktop/demo/single_cell_ATACSeq")

# Why these packages? 
 Signac: The extension of Seurat designed specifically for chromatin data.
EnsDb.Hsapiens.v75: Provides the genomic coordinates (genes/exons) for the hg19 genome.

<img width="990" height="571" alt="image" src="https://github.com/user-attachments/assets/dfe0df32-4c24-45c2-95af-54d7bb1b8826" />


# install packages
remotes::install_github("stuart-lab/signac", ref="develop")
install.packages("Matrix", type = "source")
install.packages("irlba", type = "source")
BiocManager::install("EnsDb.Hsapiens.v75")

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(tidyverse)

#  What is a fragment file? 
 It is a tab-indexed file containing every single Tn5 integration site recorded. 
We need this because the 'peak matrix' only counts reads in specific regions; 
 the fragment file allows us to calculate QC metrics like TSS enrichment from scratch.
frag.file <- read.delim('data/atac_v1_pbmc_10k_fragments.tsv.gz', header = F, nrows = 10)
head(frag.file)


# 1. Read in data 

 H5 files contain the sparse matrix of counts (Cells x Peaks).
counts <- Read10X_h5('data/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5')
counts[1:10,1:10]

Why CreateChromatinAssay? 
Unlike RNA-seq, ATAC data requires genomic ranges. This step links the count matrix 
with the fragment file and defines which genome assembly (e.g., hg19) is being used.
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "data/atac_v1_pbmc_10k_fragments.tsv.gz",
  min.cells = 10,     # Tip: Removes rare peaks that might be noise.
  min.features = 200  # Tip: Removes empty droplets or low-quality cells early.
)

str(chrom_assay)

# Metadata contains per-cell statistics (e.g., total fragments) generated by CellRanger.
metadata <- read.csv(file = 'data/atac_v1_pbmc_10k_singlecell.csv', header = T, row.names = 1)
View(metadata)


Create a Seurat Object: 
This is the container that holds counts, metadata, and later, clusters.
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  meta.data = metadata,
  assay = 'ATAC'
)

str(pbmc)


 Adding Gene Annotation 
 Why? The ATAC matrix only tells us about coordinates (e.g., chr1:100-200). 
Adding annotations allows us to see which genes are near those open regions.

pbmc@assays$ATAC@annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

Change to UCSC style (chr1) because EnsDb uses Ensembl style (1). 
If styles don't match, you won't be able to map peaks to genes.
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))

Annotation(pbmc) <- annotations
pbmc@assays$ATAC@annotation


# 2. Computing QC 


Why NucleosomeSignal? DNA wraps around nucleosomes. 
Successful ATAC-seq should show a "ladder" pattern of fragments (mononucleosomal, dinucleosomal).
pbmc <- NucleosomeSignal(pbmc)

 Why TSSEnrichment? Transcription Start Sites (TSS) are usually very open. 
High signal at TSS vs. background is the "gold standard" for ATAC-seq quality.
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

 
Why Blacklist Ratio?
Certain regions of the genome (blacklist) produce high signal 
regardless of cell type (usually due to repetitive elements). High ratios indicate noise.
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

View(pbmc@meta.data)

# Visualizing QC 
Use these plots to define your "cutoff" lines for filtering.
colnames(pbmc@meta.data)
a1 <- DensityScatter(pbmc, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
a2 <- DensityScatter(pbmc, x = 'nucleosome_signal', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

a1 | a2

VlnPlot(object = pbmc, 
        features = c('nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'nucleosome_signal', 'blacklist_ratio', 'pct_reads_in_peaks'),
        pt.size = 0.1,
        ncol = 6)


# Filtering poor quality cells 
 Tip: These thresholds are "vignette" defaults. Always adjust based on your VlnPlots above.
pbmc <- subset(x = pbmc,
                subset = nCount_ATAC > 3000 &
                 nCount_ATAC < 30000 &
                 pct_reads_in_peaks > 15 & 
                 blacklist_ratio < 0.05 &
                 nucleosome_signal < 4 &
                 TSS.enrichment > 3)


# 3. Normalization and linear dimensional reduction 

Why RunTFIDF? ATAC data is binary (open or closed). 
TF-IDF normalizes for total library size and for how common a peak is across cells.
pbmc <- RunTFIDF(pbmc) 

Why FindTopFeatures? We only use the most variable peaks for clustering to reduce noise.
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0') 

Why RunSVD? (Latent Semantic Indexing - LSI)
#This is the "PCA equivalent" for ATAC-seq. It compresses the sparse peak data.
pbmc <- RunSVD(pbmc) 

Tip: Use DepthCor to see if the first LSI component correlates with sequencing depth. 
If it does, you should exclude LSI_1 from downstream clustering (dims = 2:30).
DepthCor(pbmc)


# 4. Non-linear dimensional reduction and Clustering 

 We use dims 2:30 because LSI component 1 often captures technical variation (depth).
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)

Algorithm 3 is SLM (Smart Local Moving), which is robust for large datasets.
pbmc <- FindClusters(object = pbmc, algorithm = 3)

Visualize the result!
DimPlot(object = pbmc, label = TRUE) + NoLegend()

<img width="749" height="207" alt="image" src="https://github.com/user-attachments/assets/78680284-d739-4914-82b7-4896f739d144" />


# script to perform differential peak accesibility analysis and visualize genomic regions using single-cell ATAC-Seq data
# Vignette: https://stuartlab.org/signac/articles/pbmc_vignette
# continued from: PART 1 link: https://youtu.be/yEKZJVjc5DY?si=cm0okOcJQMwkCvPo
# setwd("~/Desktop/demo/single_cell_ATACSeq")

# install packages
# remotes::install_github("stuart-lab/signac", ref="develop")
# install.packages("Matrix", type = "source")
# install.packages("irlba", type = "source")
# BiocManager::install("EnsDb.Hsapiens.v75")

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(tidyverse)
library(SingleR)

# Pre-processed ATAC data ----------------------------------------------------
pbmc
DimPlot(object = pbmc, label = TRUE) + NoLegend()

# Create a gene activity matrix ------------------------------------------------
gene.activities <- GeneActivity(pbmc)
gene.activities[1:10,1:10]

# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc@assays
pbmc <- NormalizeData(object = pbmc,
              assay = 'RNA',
              normalization.method = 'LogNormalize',
              scale.factor = median(pbmc$nCount_RNA))


# to interpret ATAC-Seq clusters, visualizing activity of canonical marker genes
# assuming a general correspondence between gene body/promoter accessibility and gene expression which may not always be the case

DefaultAssay(pbmc) <- 'RNA'

FeaturePlot(pbmc, features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
            pt.size = 0.1,
            max.cutoff = 'q95',
            ncol = 3)


# Integrating with scRNA-Seq data. ---------------------------------------------
# Link to pre-processed RNA-Seq data: https://signac-objects.s3.amazonaws.com/pbmc_10k_v3.rds

# Load the pre-processed scRNA-seq data for PBMCs
pbmc_rna <- readRDS('data/pbmc_10k_v3.rds')
pbmc_rna <- UpdateSeuratObject(pbmc_rna)

View(pbmc_rna@meta.data)

# plot them before integrating
p1 <- DimPlot(pbmc, reduction = 'umap') + NoLegend() + ggtitle('scATAC-Seq')
p2 <- DimPlot(pbmc_rna, reduction = 'umap', group.by = 'celltype', repel = TRUE, label = TRUE) + ggtitle('scRNA-Seq') + NoLegend()

p1 | p2

# ** Should have the prior knowledge of cell types expected in your query dataset when using ref dataset


# ....Transfer Anchors by Seurat --------------
# Identify anchors

transfer.anchors <- FindTransferAnchors(reference = pbmc_rna,
                    query = pbmc,
                    reduction = 'pcaproject') # CCA is very slow


predicted.labels <- TransferData(anchorset = transfer.anchors,
             refdata = pbmc_rna$celltype,
             weight.reduction = pbmc[['lsi']],
             dims = 2:30)
head(predicted.labels)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)
View(pbmc@meta.data)


plot1 <- DimPlot(pbmc, 
        reduction = 'umap',
        group.by = 'predicted.id',
        label = TRUE,
        repel = TRUE) + NoLegend() + ggtitle('scATAC-Seq')

plot2 <- DimPlot(pbmc_rna, 
                 reduction = 'umap',
                 group.by = 'celltype',
                 label = TRUE,
                 repel = TRUE) + NoLegend() + ggtitle('scRNA-Seq')

plot1 | plot2



# Finding differentially accessible peaks between cell types -------------------
Idents(pbmc) <- pbmc$predicted.id


# change back to working with peaks instead of gene activities
DefaultAssay(pbmc) <- 'ATAC'

da_peaks <- FindMarkers(object = pbmc,
            ident.1 = 'CD4 Naive',
            ident.2 = 'CD14+ Monocytes',
            test.use = 'LR',
            latent.vars = 'nCount_ATAC')

head(da_peaks)

da_plot1 <- VlnPlot(object = pbmc,
        features = rownames(da_peaks)[1],
        pt.size = 0.1,
        idents = c('CD4 Naive','CD14+ Monocytes'))

da_plot2 <- FeaturePlot(object = pbmc,
            features = rownames(da_peaks)[1],
            pt.size = 0.1)

da_plot1 | da_plot2


# fold change between two groups of cells
fc <- FoldChange(object = pbmc, ident.1 = 'CD4 Naive', ident.2 = 'CD14+ Monocytes')
# order by fold change
fc <- fc[order(fc$avg_log2FC, decreasing = TRUE),]
head(fc)



# plotting genomic regions -----------------------------------


# set plotting order

levels(pbmc) <- unique(pbmc$predicted.id)

CoveragePlot(object = pbmc,
             region = rownames(da_peaks)[1],
             extend.upstream = 40000,
             extend.downstream = 20000)


# create interactive version of these plots?
CoverageBrowser(pbmc, region = 'CD8A')

**
