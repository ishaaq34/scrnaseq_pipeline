#!/bin/bash

# Run complete PBMC 3k analysis
# This script runs the entire pipeline from QC to final report

set -e

echo "=========================================="
echo "PBMC 3k scRNA-seq Analysis Pipeline"
echo "=========================================="
echo ""

# Check if data exists
if [ ! -d "data/cellranger_output/hg19" ]; then
    echo "Error: Data not found. Please run download_data.sh first."
    echo "Run: bash download_data.sh"
    exit 1
fi

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: Conda not found. Please install Miniconda or Anaconda."
    exit 1
fi

# Check if snakemake is available
if ! command -v snakemake &> /dev/null; then
    echo "Error: Snakemake not found. Please install it in your conda environment."
    echo "Run: conda install -c conda-forge -c bioconda snakemake"
    exit 1
fi

# Set number of cores (adjust based on your system)
CORES=${1:-4}

echo "Configuration:"
echo "  - Cores: $CORES"
echo "  - Config: config_pbmc3k.yaml"
echo "  - Framework: $(grep '^framework:' config_pbmc3k.yaml | awk '{print $2}')"
echo ""

# Run dry run first
echo "Running dry run to check workflow..."
snakemake --dry-run --cores 1 --configfile config_pbmc3k.yaml

echo ""
echo "Dry run successful! Starting analysis..."
echo ""

# Run the pipeline
snakemake \
    --cores $CORES \
    --use-conda \
    --configfile config_pbmc3k.yaml \
    --printshellcmds \
    --keep-going

echo ""
echo "=========================================="
echo "Analysis Complete!"
echo "=========================================="
echo ""
echo "Results are available in:"
echo "  - QC metrics: results/qc/"
echo "  - Clustering: results/clustering/"
echo "  - Cell types: results/annotation/"
echo "  - Final report: results/reports/analysis_report.html"
echo ""
echo "To view the report:"
echo "  open results/reports/analysis_report.html"
echo ""
