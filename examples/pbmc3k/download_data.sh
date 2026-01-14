#!/bin/bash

# Download PBMC 3k dataset from 10x Genomics
# Dataset: 2,700 Peripheral Blood Mononuclear Cells from a Healthy Donor
# Source: https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k

set -e

echo "Downloading PBMC 3k dataset from 10x Genomics..."

# Create data directory
mkdir -p data
cd data

# Download filtered gene-barcode matrices
echo "Downloading filtered gene-barcode matrices..."
wget -O pbmc3k_filtered_gene_bc_matrices.tar.gz \
    https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

# Extract
echo "Extracting data..."
tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz

# Organize files
echo "Organizing files..."
mkdir -p cellranger_output
mv filtered_gene_bc_matrices/hg19 cellranger_output/

# Convert to H5 format if needed (requires cellranger or scanpy)
if command -v cellranger &> /dev/null; then
    echo "Converting to H5 format..."
    # This would require cellranger
    echo "Note: H5 conversion requires Cell Ranger or manual conversion"
fi

# Clean up
rm pbmc3k_filtered_gene_bc_matrices.tar.gz
rm -rf filtered_gene_bc_matrices

echo ""
echo "Download complete!"
echo "Data location: data/cellranger_output/hg19/"
echo ""
echo "Files downloaded:"
ls -lh cellranger_output/hg19/

echo ""
echo "Next steps:"
echo "1. Run the pipeline: snakemake --cores 4 --use-conda"
echo "2. Or use the pre-configured workflow: bash run_analysis.sh"
