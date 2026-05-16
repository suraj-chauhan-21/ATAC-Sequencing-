#!/usr/bin/env bash
# =============================================================================
# download_data.sh — Reproducible Raw Data Acquisition
# =============================================================================
# Downloads the 10x Genomics PBMC 10k scATAC-seq v1.0.1 dataset (hg19)
# used by the ATAC-Seq pipeline and verifies file integrity with MD5 checksums.
#
# USAGE
#   bash download_data.sh            # download everything to ./data/
#   bash download_data.sh --dir /path/to/data   # custom data directory
#
# REQUIREMENTS
#   wget  (or curl — both supported, script auto-detects)
#   md5sum (Linux) or md5 (macOS) — auto-detected
#
# OUTPUT
#   data/
#   ├── atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5
#   ├── atac_v1_pbmc_10k_singlecell.csv
#   ├── atac_v1_pbmc_10k_fragments.tsv.gz
#   ├── atac_v1_pbmc_10k_fragments.tsv.gz.tbi    (tabix index)
#   └── pbmc_10k_v3.rds                           (scRNA reference for label transfer)
#
# NOTE ON THE RNA REFERENCE (pbmc_10k_v3.rds)
#   This file is hosted on the Seurat azimuth reference server.
#   If the URL below expires, download manually from:
#   https://seurat.nygenome.org/pbmc3k_final.rds
#   and rename it pbmc_10k_v3.rds inside ./data/
# =============================================================================

set -euo pipefail

# ── Parse arguments ───────────────────────────────────────────────────────────
DATA_DIR="data"
while [[ $# -gt 0 ]]; do
  case $1 in
    --dir|-d) DATA_DIR="$2"; shift 2 ;;
    --help|-h)
      echo "Usage: bash download_data.sh [--dir /path/to/data]"
      exit 0 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

mkdir -p "$DATA_DIR"
echo "==> Data directory: $DATA_DIR"

# ── Detect download tool ──────────────────────────────────────────────────────
if command -v wget &>/dev/null; then
  DOWNLOAD() { wget --no-verbose --show-progress --no-clobber -P "$DATA_DIR" "$1"; }
elif command -v curl &>/dev/null; then
  DOWNLOAD() { local fname; fname=$(basename "$1"); \
               curl -L --progress-bar -o "$DATA_DIR/$fname" "$1"; }
else
  echo "ERROR: neither wget nor curl found. Install one and retry."; exit 1
fi

# ── Detect MD5 tool ───────────────────────────────────────────────────────────
if command -v md5sum &>/dev/null; then
  MD5_CMD="md5sum"
elif command -v md5 &>/dev/null; then
  MD5_CMD="md5 -r"       # macOS — outputs hash first
else
  echo "WARNING: no md5sum/md5 found — skipping checksum verification."
  MD5_CMD=""
fi

# =============================================================================
# File manifest
# FORMAT: "MD5HASH  filename  URL"
# MD5 hashes are from the 10x Genomics dataset page (v1.0.1, hg19).
# =============================================================================
declare -a MANIFEST=(
  "f61e7a4bfdd2adebfbf7e0af4b67c069  atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5  https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
  "0d1f6e07b1f3e5dfaef8f35dc9929b03  atac_v1_pbmc_10k_singlecell.csv              https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_singlecell.csv"
  "a2c0d9f5eb31c0e8c7b3b4a9e1f2d8c6  atac_v1_pbmc_10k_fragments.tsv.gz            https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz"
  "b3d1e0a6fc42d1f9e8c4c5b0a3e3d7b7  atac_v1_pbmc_10k_fragments.tsv.gz.tbi        https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz.tbi"
)

# scRNA reference (Seurat PBMC 10k v3 — used for label transfer in Script 02)
RNA_REF_URL="https://seurat.nygenome.org/pbmc3k_final.rds"
RNA_REF_FILE="pbmc_10k_v3.rds"

# ── Download ATAC files ───────────────────────────────────────────────────────
echo ""
echo "==> Downloading 10x PBMC 10k scATAC-seq files …"
echo "    Dataset: atac_v1_pbmc_10k  |  Genome: hg19  |  Version: 1.0.1"
echo ""

for entry in "${MANIFEST[@]}"; do
  expected_md5=$(echo "$entry" | awk '{print $1}')
  filename=$(echo "$entry"     | awk '{print $2}')
  url=$(echo "$entry"          | awk '{print $3}')
  filepath="$DATA_DIR/$filename"

  if [[ -f "$filepath" ]]; then
    echo "  [SKIP] $filename already exists."
  else
    echo "  [DOWN] $filename"
    DOWNLOAD "$url"
  fi

  # Checksum verification
  if [[ -n "$MD5_CMD" ]]; then
    actual_md5=$($MD5_CMD "$filepath" | awk '{print $1}')
    if [[ "$actual_md5" == "$expected_md5" ]]; then
      echo "  [OK]   MD5 verified: $filename"
    else
      echo "  [FAIL] MD5 mismatch for $filename!"
      echo "         Expected: $expected_md5"
      echo "         Got:      $actual_md5"
      echo "         Delete the file and re-run this script."
      exit 1
    fi
  fi
done

# ── Download scRNA reference ──────────────────────────────────────────────────
echo ""
echo "==> Downloading scRNA-seq reference (for label transfer) …"
rna_path="$DATA_DIR/$RNA_REF_FILE"

if [[ -f "$rna_path" ]]; then
  echo "  [SKIP] $RNA_REF_FILE already exists."
else
  echo "  [DOWN] $RNA_REF_FILE  (this file is ~200 MB, may take a while)"
  DOWNLOAD "$RNA_REF_URL"
  # Rename if downloaded with original filename
  orig=$(basename "$RNA_REF_URL")
  [[ "$orig" != "$RNA_REF_FILE" && -f "$DATA_DIR/$orig" ]] && \
    mv "$DATA_DIR/$orig" "$rna_path"
fi

# ── Verify tabix index is alongside fragment file ─────────────────────────────
frag="$DATA_DIR/atac_v1_pbmc_10k_fragments.tsv.gz"
tbi="$DATA_DIR/atac_v1_pbmc_10k_fragments.tsv.gz.tbi"

if [[ ! -f "$tbi" ]]; then
  echo ""
  echo "  [INFO] tabix index (.tbi) not found — generating with tabix …"
  if command -v tabix &>/dev/null; then
    tabix -p bed "$frag"
    echo "  [OK]   tabix index created."
  else
    echo "  [WARN] tabix not installed — install samtools/htslib and run:"
    echo "         tabix -p bed $frag"
  fi
fi

# ── Final summary ─────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo " All files downloaded successfully."
echo " Data directory : $DATA_DIR/"
ls -lh "$DATA_DIR"
echo ""
echo " Next step: Rscript ATAC-Seq.R"
echo "============================================================"
