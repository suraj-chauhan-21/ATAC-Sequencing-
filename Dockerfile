# =============================================================================
# Dockerfile — ATAC-Seq Pipeline (Signac + Seurat)
# =============================================================================
# Produces a fully self-contained container for reproducible scATAC-seq
# analysis.  Based on the official Bioconductor Docker image which bundles
# R, all Bioconductor system dependencies, and a matching R version.
#
# BUILD
# -----
#   docker build -t atac_seq:1.0 .
#
# RUN — interactive R session
#   docker run --rm -it \
#     -v $(pwd)/data:/home/rstudio/data \
#     -v $(pwd)/results:/home/rstudio/results \
#     atac_seq:1.0 R
#
# RUN — execute the full pipeline non-interactively
#   docker run --rm \
#     -v $(pwd)/data:/home/rstudio/data \
#     -v $(pwd)/results:/home/rstudio/results \
#     atac_seq:1.0 Rscript /pipeline/ATAC-Seq.R
#
# RUN — RStudio Server (browser at http://localhost:8787)
#   docker run --rm -p 8787:8787 \
#     -e PASSWORD=yourpassword \
#     -v $(pwd):/home/rstudio/project \
#     atac_seq:1.0
#
# NOTES
# -----
#   • Data files are NOT baked into the image — mount them as volumes.
#   • Image size ~4–5 GB due to Bioconductor + genome annotation packages.
#   • For HPC / Singularity:  singularity pull docker://atac_seq:1.0
# =============================================================================

# ── Base image ────────────────────────────────────────────────────────────────
# bioconductor/bioconductor_docker bundles:
#   • R 4.3.x  (matched to Bioconductor 3.18)
#   • RStudio Server
#   • All system libraries needed by Bioconductor (libcurl, libxml2, etc.)
FROM bioconductor/bioconductor_docker:RELEASE_3_18

LABEL maintainer="suraj-chauhan-21"
LABEL description="Reproducible scATAC-seq pipeline: Signac + Seurat + Bioconductor 3.18"
LABEL version="1.0"

# ── System dependencies ───────────────────────────────────────────────────────
# libhdf5-dev   : required by hdf5r (Read10X_h5)
# libbz2-dev    : required by Rsamtools / htslib
# liblzma-dev   : required by htslib
# libcurl4      : required by BiocManager download
# tabix/samtools: command-line genomics tools for fragment-file operations
RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    samtools \
    tabix \
    bedtools \
    wget \
    curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ── R package installation ────────────────────────────────────────────────────
# Install in a single RUN layer to keep the image layer count low.
# Use BiocManager for all Bioconductor packages to ensure version compatibility.
RUN Rscript -e " \
  options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/jammy/latest')); \
  \
  # CRAN packages \
  install.packages(c( \
    'Seurat', \
    'tidyverse', \
    'patchwork', \
    'cowplot', \
    'yaml', \
    'hdf5r', \
    'remotes', \
    'renv' \
  ), dependencies = TRUE); \
  \
  # Bioconductor packages \
  BiocManager::install(c( \
    'Signac', \
    'EnsDb.Hsapiens.v75', \
    'BSgenome.Hsapiens.UCSC.hg19', \
    'JASPAR2020', \
    'TFBSTools', \
    'chromVAR', \
    'GenomicRanges', \
    'IRanges', \
    'BiocGenerics' \
  ), ask = FALSE, update = FALSE); \
"

# ── Copy pipeline scripts ─────────────────────────────────────────────────────
# The /pipeline directory holds all R scripts and config.yaml.
# Data is NOT copied — it must be mounted at runtime.
WORKDIR /pipeline
COPY ATAC-Seq.R                  ./ATAC-Seq.R
COPY 02_label_transfer_and_DA.R  ./02_label_transfer_and_DA.R
COPY config.yaml                 ./config.yaml

# ── Create expected directory structure ───────────────────────────────────────
RUN mkdir -p /pipeline/data /pipeline/results

# ── Default command: launch R ─────────────────────────────────────────────────
# Override with `docker run atac_seq:1.0 Rscript ATAC-Seq.R` to run headlessly.
CMD ["R"]

# =============================================================================
# SINGULARITY (HPC) USAGE
# =============================================================================
# Build image locally, push to Docker Hub, then on the HPC:
#
#   singularity pull atac_seq.sif docker://yourhubuser/atac_seq:1.0
#
#   singularity exec \
#     --bind /scratch/mydata:/pipeline/data \
#     --bind /scratch/results:/pipeline/results \
#     atac_seq.sif Rscript /pipeline/ATAC-Seq.R
# =============================================================================
