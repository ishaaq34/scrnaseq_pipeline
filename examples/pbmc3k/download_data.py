#!/usr/bin/env python3

"""
Download and prepare PBMC 3k dataset for analysis
This script downloads the dataset from 10x Genomics and converts it to H5AD format
"""

import os
import sys
import urllib.request
import tarfile
from pathlib import Path

def download_file(url, output_path):
    """Download file with progress bar"""
    print(f"Downloading from {url}...")
    
    def progress_hook(count, block_size, total_size):
        percent = int(count * block_size * 100 / total_size)
        sys.stdout.write(f"\rProgress: {percent}% ")
        sys.stdout.flush()
    
    urllib.request.urlretrieve(url, output_path, progress_hook)
    print("\nDownload complete!")

def extract_tarfile(tar_path, extract_path):
    """Extract tar.gz file"""
    print(f"Extracting {tar_path}...")
    with tarfile.open(tar_path, 'r:gz') as tar:
        tar.extractall(extract_path)
    print("Extraction complete!")

def convert_to_h5ad():
    """Convert to H5AD format using scanpy"""
    try:
        import scanpy as sc
        import numpy as np
        
        print("Converting to H5AD format...")
        
        # Read 10x data
        data_dir = "data/cellranger_output/hg19"
        adata = sc.read_10x_mtx(
            data_dir,
            var_names='gene_symbols',
            cache=True
        )
        
        # Make variable names unique
        adata.var_names_make_unique()
        
        # Add basic metadata
        adata.obs['sample'] = 'pbmc3k'
        adata.obs['condition'] = 'healthy'
        adata.obs['batch'] = 'batch1'
        
        # Save as H5AD
        output_file = "data/pbmc3k.h5ad"
        adata.write(output_file)
        print(f"Saved to {output_file}")
        
        # Print summary
        print(f"\nDataset summary:")
        print(f"  Cells: {adata.n_obs}")
        print(f"  Genes: {adata.n_vars}")
        
        return True
        
    except ImportError:
        print("Warning: scanpy not installed. Skipping H5AD conversion.")
        print("You can still use the MTX format with the pipeline.")
        return False

def main():
    """Main function"""
    print("=" * 60)
    print("PBMC 3k Dataset Download and Preparation")
    print("=" * 60)
    print()
    
    # Create data directory
    data_dir = Path("data")
    data_dir.mkdir(exist_ok=True)
    
    # Download URL
    url = "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
    tar_file = data_dir / "pbmc3k_filtered_gene_bc_matrices.tar.gz"
    
    # Download if not exists
    if not tar_file.exists():
        download_file(url, tar_file)
    else:
        print(f"File already exists: {tar_file}")
    
    # Extract
    print()
    extract_tarfile(tar_file, data_dir)
    
    # Organize files
    print("\nOrganizing files...")
    cellranger_dir = data_dir / "cellranger_output"
    cellranger_dir.mkdir(exist_ok=True)
    
    # Move extracted data
    source = data_dir / "filtered_gene_bc_matrices" / "hg19"
    dest = cellranger_dir / "hg19"
    
    if source.exists() and not dest.exists():
        source.rename(dest)
        print(f"Moved data to {dest}")
    
    # Clean up
    print("\nCleaning up...")
    if tar_file.exists():
        tar_file.unlink()
        print(f"Removed {tar_file}")
    
    temp_dir = data_dir / "filtered_gene_bc_matrices"
    if temp_dir.exists():
        import shutil
        shutil.rmtree(temp_dir)
        print(f"Removed {temp_dir}")
    
    # Convert to H5AD
    print()
    convert_to_h5ad()
    
    # Print summary
    print()
    print("=" * 60)
    print("Setup Complete!")
    print("=" * 60)
    print()
    print("Data location: data/cellranger_output/hg19/")
    print()
    print("Files:")
    if dest.exists():
        for f in dest.iterdir():
            size = f.stat().st_size / (1024 * 1024)  # MB
            print(f"  - {f.name} ({size:.2f} MB)")
    
    print()
    print("Next steps:")
    print("  1. Run the pipeline: bash run_analysis.sh")
    print("  2. Or manually: snakemake --cores 4 --use-conda --configfile config_pbmc3k.yaml")
    print()

if __name__ == "__main__":
    main()
