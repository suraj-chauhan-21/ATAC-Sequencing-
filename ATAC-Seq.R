# =============================================================================
# ATAC-Seq Analysis Pipeline — Script 01: Preprocessing, QC & Clustering
# =============================================================================
# Author  : suraj-chauhan-21
# Framework: Signac + Seurat (R)
# Dataset : 10x Genomics PBMC 10k scATAC-seq (hg19)
# Usage   : Rscript ATAC-Seq.R
#           (or source() interactively)
# =============================================================================

# ── 0. Load config & libraries ──────────────────────────────────────────────

# All tunable parameters live in config.yaml — never edit thresholds here.
library(yaml)
cfg <- yaml::read_yaml("config.yaml")

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(EnsDb.Hsapiens.v75)
  library(tidyverse)
  library(patchwork)
})

set.seed(42)                          # reproducible UMAP / clustering
dir.create("results", showWarnings = FALSE)

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  flush.console()
}

log_msg("Pipeline started.")

# ── 1. Read input data ───────────────────────────────────────────────────────
#
# WHY H5?  The filtered_peak_bc_matrix.h5 stores the sparse cell × peak count
# matrix in HDF5 format — efficient for large matrices. Read10X_h5() parses it
# into a standard dgCMatrix (sparse integer matrix).

log_msg("Reading peak-count matrix …")
counts <- Read10X_h5(cfg$input$peaks_h5)

# Quick sanity check
stopifnot(
  "Peak matrix is empty — check the H5 path in config.yaml" =
    nrow(counts) > 0 && ncol(counts) > 0
)
log_msg(sprintf("  %d peaks × %d barcodes detected.", nrow(counts), ncol(counts)))

# ── 2. Build the ChromatinAssay ──────────────────────────────────────────────
#
# WHY CreateChromatinAssay (not CreateAssayObject)?
# ATAC data is fundamentally *genomic*: each feature is a chromosomal region,
# not a gene name.  ChromatinAssay links the count matrix to:
#   • the fragment file (needed for all QC metrics)
#   • the genome assembly (hg19 here)
#   • genomic ranges objects for downstream peak operations
#
# min.cells  : drops peaks present in < 10 cells (likely noise / artefacts)
# min.features: drops barcodes with < 200 peaks detected (empty droplets)

log_msg("Creating ChromatinAssay …")
chrom_assay <- CreateChromatinAssay(
  counts       = counts,
  sep          = c(":", "-"),           # peak names are "chr:start-end"
  fragments    = cfg$input$fragments,
  min.cells    = 10,
  min.features = 200
)

# ── 3. Create Seurat object + metadata ──────────────────────────────────────
#
# CellRanger generates per-barcode statistics in singlecell.csv:
# total fragments, peak fragments, TSS overlaps, blacklist overlaps, etc.
# We load these as the initial metadata so QC metrics are available immediately.

log_msg("Creating Seurat object …")
metadata <- read.csv(cfg$input$metadata, header = TRUE, row.names = 1)

pbmc <- CreateSeuratObject(
  counts    = chrom_assay,
  meta.data = metadata,
  assay     = "ATAC"
)

log_msg(sprintf("  Seurat object created: %d cells.", ncol(pbmc)))

# ── 4. Add gene annotation ───────────────────────────────────────────────────
#
# The ATAC matrix only contains coordinates (chr1:100-200).
# Adding EnsDb annotations allows:
#   • TSSEnrichment()  — needs TSS positions
#   • GeneActivity()   — needs gene body coordinates (Script 02)
#   • CoveragePlot()   — annotates genes on coverage tracks
#
# IMPORTANT: EnsDb uses Ensembl-style chromosome names ("1", "2" …).
# Seurat/Signac expects UCSC style ("chr1", "chr2" …).
# seqlevelsStyle() conversion is mandatory or coordinate matching silently fails.

log_msg("Adding gene annotations (EnsDb.Hsapiens.v75 / hg19) …")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotations) <- "UCSC"   # convert "1" → "chr1"
Annotation(pbmc) <- annotations

# ── 5. Compute QC metrics ────────────────────────────────────────────────────
#
# TSS Enrichment  : ratio of signal at TSSs vs distal regions.
#                   ATAC-seq from healthy cells shows sharp enrichment at TSS.
#                   Values < 3 indicate poor tagmentation or dead cells.
#
# Nucleosome Signal: ratio of mono-nucleosomal to sub-nucleosomal fragments.
#                   Low values (< 4) confirm the nucleosomal ladder is intact.
#
# Blacklist Ratio : ENCODE blacklist regions produce artifactual signal
#                   regardless of cell type (repetitive elements, centromeres).
#                   Cells with > 5% reads in blacklist are discarded.
#
# Pct reads in peaks: library specificity.  Good libraries: > 15%.

log_msg("Computing nucleosome signal …")
pbmc <- NucleosomeSignal(pbmc)

log_msg("Computing TSS enrichment (fast = FALSE for accuracy) …")
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

pbmc$blacklist_ratio   <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

log_msg("QC metrics computed.")

# ── 6. Visualise QC — inspect BEFORE filtering ──────────────────────────────
#
# TIP: Look at the density scatter plots.  The quantile lines show you where
# most cells fall.  Your filter cutoffs should cleanly separate the main cloud
# from the low-quality tails — adjust config.yaml accordingly.

log_msg("Plotting QC metrics …")

a1 <- DensityScatter(pbmc, x = "nCount_ATAC", y = "TSS.enrichment",
                     log_x = TRUE, quantiles = TRUE) +
      ggtitle("Fragment count vs TSS Enrichment")

a2 <- DensityScatter(pbmc, x = "nucleosome_signal", y = "TSS.enrichment",
                     log_x = TRUE, quantiles = TRUE) +
      ggtitle("Nucleosome Signal vs TSS Enrichment")

vln <- VlnPlot(
  object   = pbmc,
  features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment",
                "nucleosome_signal", "blacklist_ratio", "pct_reads_in_peaks"),
  pt.size  = 0.1,
  ncol     = 6
)

ggsave("results/QC_density_scatter.png", a1 | a2, width = 14, height = 6, dpi = 150)
ggsave("results/QC_violin.png",          vln,      width = 18, height = 5, dpi = 150)

log_msg("  QC plots saved to results/")

# ── 7. Filter low-quality cells ──────────────────────────────────────────────
#
# Thresholds come from config.yaml — do NOT hard-code here.
# After filtering, print a summary so the log records how many cells pass.

n_before <- ncol(pbmc)

pbmc <- subset(
  x      = pbmc,
  subset =
    nCount_ATAC        >  cfg$qc$min_count       &
    nCount_ATAC        <  cfg$qc$max_count        &
    pct_reads_in_peaks >  cfg$qc$min_pct_peaks    &
    blacklist_ratio    <  cfg$qc$max_blacklist     &
    nucleosome_signal  <  cfg$qc$max_nuc_signal   &
    TSS.enrichment     >  cfg$qc$min_tss
)

n_after <- ncol(pbmc)
log_msg(sprintf("QC filtering: %d → %d cells retained (%.1f%% passed).",
                n_before, n_after, 100 * n_after / n_before))

if (n_after < 100) {
  warning("Fewer than 100 cells remain after QC — check your thresholds in config.yaml.")
}

# ── 8. Normalisation: TF-IDF ─────────────────────────────────────────────────
#
# WHY TF-IDF and not log-normalisation?
# ATAC data is binary-like (open = 1, closed = 0) and extremely sparse.
# TF (term frequency) = how often peak i is open in cell j.
# IDF (inverse document frequency) = penalises peaks open in MANY cells
#   (ubiquitous housekeeping regions) and up-weights rare, informative peaks.
# Together, TF-IDF removes library-size bias and highlights discriminative peaks.

log_msg("Running TF-IDF normalisation …")
pbmc <- RunTFIDF(pbmc)

# ── 9. Feature selection ─────────────────────────────────────────────────────
#
# min.cutoff = "q0" keeps ALL peaks (no variance filtering).
# For large datasets consider "q5" (top 95% most variable) to reduce noise.

log_msg("Selecting top features …")
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")

# ── 10. Dimensionality reduction: LSI ────────────────────────────────────────
#
# WHY SVD / LSI instead of PCA?
# PCA assumes data is normally distributed and works on dense matrices.
# ATAC matrices are sparse and binary — LSI (TF-IDF + truncated SVD) is the
# information-retrieval equivalent of PCA and is standard for ATAC-seq.
#
# After RunSVD(), always check DepthCor():
#   If LSI_1 correlates strongly (|r| > 0.75) with sequencing depth, exclude it.
#   We exclude it by default (dims 2:30) following the Signac vignette.

log_msg("Running SVD (LSI) …")
pbmc <- RunSVD(pbmc)

depth_cor_plot <- DepthCor(pbmc)
ggsave("results/LSI_depth_correlation.png", depth_cor_plot, width = 8, height = 4, dpi = 150)
log_msg("  Depth-correlation plot saved — verify LSI_1 is excluded from clustering.")

lsi_dims <- cfg$clustering$lsi_dims[1]:cfg$clustering$lsi_dims[2]

# ── 11. Non-linear reduction: UMAP ───────────────────────────────────────────

log_msg("Running UMAP …")
pbmc <- RunUMAP(object = pbmc, reduction = "lsi", dims = lsi_dims)

# ── 12. Graph-based clustering ───────────────────────────────────────────────
#
# Algorithm 3 = SLM (Smart Local Moving) — robust and fast for large datasets.
# Resolution controls granularity: higher → more clusters.

log_msg("Finding neighbours and clusters …")
pbmc <- FindNeighbors(object = pbmc, reduction = "lsi", dims = lsi_dims)
pbmc <- FindClusters(object  = pbmc,
                     algorithm  = cfg$clustering$algorithm,
                     resolution = cfg$clustering$resolution)

n_clusters <- length(unique(Idents(pbmc)))
log_msg(sprintf("  %d clusters identified at resolution %.2f.",
                n_clusters, cfg$clustering$resolution))

# ── 13. Visualise clusters ───────────────────────────────────────────────────

log_msg("Plotting UMAP …")
umap_plot <- DimPlot(object = pbmc, label = TRUE) + NoLegend() +
             ggtitle(sprintf("scATAC-seq — %d cells, %d clusters", ncol(pbmc), n_clusters))

ggsave("results/UMAP_clusters.png", umap_plot, width = 8, height = 7, dpi = 150)

# ── 14. Save processed object & session info ─────────────────────────────────

log_msg("Saving processed Seurat object …")
saveRDS(pbmc, "results/pbmc_atac_processed.rds")

log_msg("Saving session info …")
writeLines(capture.output(sessionInfo()), "results/sessionInfo_01.txt")

log_msg("=== Script 01 complete.  Run 02_label_transfer_and_DA.R next. ===")
