# ATAC-Sequencing

## Overview

This repository contains a reproducible workflow for analyzing **10x Genomics single-cell ATAC-seq (scATAC-seq)** data using the **Signac + Seurat** framework in R.

The goal is to identify open chromatin regions, perform quality control, reduce dimensionality using LSI, cluster cells, and integrate chromatin accessibility with gene activity for biological interpretation.

---

## Conceptual Background

ATAC-seq measures **chromatin accessibility**.

* Each genomic region behaves like a regulatory “door”.
* **Closed chromatin** → transcriptionally inactive.
* **Open chromatin** → regulatory potential or active transcription.
* The Tn5 transposase inserts sequencing adapters preferentially into open regions, allowing us to map accessible DNA genome-wide.

This enables inference of:

* Regulatory landscapes
* Cell identity
* Transcription factor activity
* Epigenomic heterogeneity at single-cell resolution

---

## Workflow Summary

The analysis is divided into two scripts:

| Step | Script                        | Description                                                |
| ---- | ----------------------------- | ---------------------------------------------------------- |
| 1    | `01_scATAC_signac_pipeline.R` | Preprocessing, QC, LSI, clustering                         |
| 2    | `02_label_transfer_and_DA.R`  | Gene activity, RNA integration, differential accessibility |

---

## Repository Structure

```
ATAC-Sequencing/
│
├── 01_scATAC_signac_pipeline.R
├── 02_label_transfer_and_DA.R
├── data/                         # Input files go here
├── results/                      # Auto-generated outputs
└── README.md
```

---

## Input Data Requirements

Place the following 10x Genomics outputs inside the `data/` directory:

```
data/
├── atac_v1_pbmc_10k_fragments.tsv.gz
├── atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5
├── atac_v1_pbmc_10k_singlecell.csv
└── pbmc_10k_v3.rds   (scRNA reference for label transfer)
```

---

## Installation

Install required packages in R:

```r
install.packages(c("Seurat", "tidyverse", "patchwork"))
install.packages("BiocManager")
BiocManager::install("EnsDb.Hsapiens.v75")

# Signac (development version)
install.packages("remotes")
remotes::install_github("stuart-lab/signac", ref = "develop")
```

---

## Running the Pipeline

From the project directory:

### Step 1 — Preprocess scATAC Data

```bash
Rscript 01_scATAC_signac_pipeline.R
```

This performs:

* Fragment parsing
* Chromatin assay construction
* TSS enrichment & nucleosome QC
* TF-IDF normalization
* Latent Semantic Indexing (LSI)
* UMAP embedding
* Clustering

Outputs saved in:

```
results/
├── pbmc_atac_processed.rds
├── UMAP_clusters.png
└── sessionInfo.txt
```

---

### Step 2 — Biological Interpretation

```bash
Rscript 02_label_transfer_and_DA.R
```

This performs:

* Gene activity inference
* Integration with scRNA-seq reference
* Cell type annotation via label transfer
* Differential accessibility analysis
* Coverage visualization

Outputs:

```
results/
├── pbmc_atac_annotated.rds
├── UMAP_predicted_labels.png
├── DA_peaks_CD4_vs_Mono.csv
└── Coverage_top_DA_peak.png
```

---

## Key Methods Implemented

### Quality Control Metrics

* **TSS Enrichment** — signal-to-noise indicator
* **Nucleosome Signal** — fragment periodicity check
* **Blacklist Ratio** — removal of artefactual regions
* **Fragments in Peaks (%)** — library specificity

### Dimensional Reduction

TF-IDF + SVD (Latent Semantic Indexing) is used instead of PCA because ATAC-seq data is:

* Sparse
* Binary-like
* Region-based rather than gene-based

### Integration Strategy

Gene activity scores derived from accessibility are aligned to a reference scRNA-seq atlas to transfer biological labels.

---

## Customization

QC thresholds can be tuned inside Script-01:

```
min_fragments
max_fragments
min_TSS
max_nuc
min_peak_pct
```

These should be dataset-specific.

---

## Reproducibility Features

* No hard-coded paths
* Session info automatically saved
* Relative directory structure
* Modular execution
* Compatible with HPC or local systems

---

## Recommended Citation of Methods

If using this workflow, please cite:

* Signac framework
* Seurat single-cell integration methodology
* 10x Genomics scATAC technology description

---

## Notes

This repository is designed as:

* A learning resource for chromatin accessibility analysis
* A reproducible starting point for new datasets
* A modular template adaptable to multi-omic studies

---

## Future Extensions

Planned additions:

* Motif enrichment analysis
* Peak-to-gene linkage
* Pseudotime chromatin dynamics
* Multiome integration support

---

## Contact

For questions, suggestions, or collaboration inquiries, feel free to reach out.

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




**
