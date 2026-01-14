#!/usr/bin/env python3

"""
Quality Control Analysis with Scanpy
=====================================
"""

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set plotting parameters
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white', frameon=False)

# Parse snakemake inputs
samples_file = snakemake.params["samples_file"]
min_genes = snakemake.params["min_genes"]
max_genes = snakemake.params["max_genes"]
max_mito = snakemake.params["max_mito"]
min_cells = snakemake.params["min_cells"]

# Output files
qc_metrics_file = snakemake.output["metrics"]
qc_plots_file = snakemake.output["plots"]
adata_file = snakemake.output["adata"]

# Load sample metadata
samples_df = pd.read_csv(samples_file, sep='\t')

# Initialize list to store AnnData objects
adata_list = []

# Load data for each sample
for idx, row in samples_df.iterrows():
    sample_id = row['sample_id']
    condition = row['condition']
    batch = row['batch']
    
    print(f"Loading sample: {sample_id}")
    
    # Load 10x data
    h5_file = f"results/cellranger/{sample_id}/outs/filtered_feature_bc_matrix.h5"
    
    if Path(h5_file).exists():
        adata = sc.read_10x_h5(h5_file)
    else:
        # Alternative: load from matrix directory
        data_dir = f"results/cellranger/{sample_id}/outs/filtered_feature_bc_matrix"
        adata = sc.read_10x_mtx(data_dir, var_names='gene_symbols', cache=True)
    
    # Make variable names unique
    adata.var_names_make_unique()
    
    # Add metadata
    adata.obs['sample'] = sample_id
    adata.obs['condition'] = condition
    adata.obs['batch'] = batch
    
    adata_list.append(adata)

# Concatenate all samples
adata = sc.concat(adata_list, label="sample", keys=[row['sample_id'] for _, row in samples_df.iterrows()])

# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.match('^RP[SL]')

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=['mt', 'ribo'],
    percent_top=None,
    log1p=False,
    inplace=True
)

# Extract QC metrics before filtering
qc_metrics = pd.DataFrame({
    'cell_id': adata.obs_names,
    'sample': adata.obs['sample'],
    'condition': adata.obs['condition'],
    'batch': adata.obs['batch'],
    'n_genes_by_counts': adata.obs['n_genes_by_counts'],
    'total_counts': adata.obs['total_counts'],
    'pct_counts_mt': adata.obs['pct_counts_mt'],
    'pct_counts_ribo': adata.obs['pct_counts_ribo']
})

# Save pre-filtering metrics
qc_metrics.to_csv(qc_metrics_file, index=False)

# Generate QC plots
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Violin plots
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, ax=axes[0], show=False)

# Scatter plots
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=axes[1, 0], show=False)
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=axes[1, 1], show=False)

# Histogram
axes[1, 2].hist(adata.obs['n_genes_by_counts'], bins=50, color='steelblue', edgecolor='black')
axes[1, 2].axvline(min_genes, color='red', linestyle='--', label=f'Min: {min_genes}')
axes[1, 2].axvline(max_genes, color='red', linestyle='--', label=f'Max: {max_genes}')
axes[1, 2].set_xlabel('Number of Genes')
axes[1, 2].set_ylabel('Count')
axes[1, 2].set_title('Distribution of Genes per Cell')
axes[1, 2].legend()

plt.tight_layout()
plt.savefig(qc_plots_file)
plt.close()

# Apply QC filters
n_cells_before = adata.n_obs

sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

# Additional filters
adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
adata = adata[adata.obs.pct_counts_mt < max_mito, :]

n_cells_after = adata.n_obs

# Print filtering summary
print("\n=== QC Filtering Summary ===")
print(f"Cells before filtering: {n_cells_before}")
print(f"Cells after filtering: {n_cells_after}")
print(f"Cells removed: {n_cells_before - n_cells_after} ({100 * (n_cells_before - n_cells_after) / n_cells_before:.1f}%)")

# Save filtered AnnData object
adata.write(adata_file)

print("\nQC analysis complete!")
