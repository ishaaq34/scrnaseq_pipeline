# Quick Start Guide

This guide will walk you through running your first scRNA-seq analysis with the pipeline.

## Prerequisites

- Pipeline installed (see [Installation Guide](installation.md))
- Conda environment activated: `conda activate scrnaseq`
- Input data prepared (FASTQ files or Cell Ranger outputs)

## Step 1: Prepare Your Data

### Option A: Starting from Cell Ranger Outputs

If you already have Cell Ranger outputs:

```bash
# Create data directory structure
mkdir -p results/cellranger/SAMPLE_NAME/outs/

# Copy your Cell Ranger outputs
cp /path/to/cellranger/SAMPLE_NAME/outs/filtered_feature_bc_matrix.h5 \
   results/cellranger/SAMPLE_NAME/outs/
```

### Option B: Starting from FASTQ Files

If you have raw FASTQ files:

```bash
# Organize FASTQ files
mkdir -p data/fastqs/SAMPLE_NAME/
cp /path/to/fastqs/*.fastq.gz data/fastqs/SAMPLE_NAME/

# Enable Cell Ranger in config
# Set run_cellranger: true in config/config.yaml
```

## Step 2: Configure the Pipeline

### Edit Sample Metadata

Edit `config/samples.tsv`:

```tsv
sample_id       condition       batch   fastq_dir       notes
PBMC_1          control         batch1  data/fastqs/PBMC_1      Control sample
PBMC_2          treated         batch1  data/fastqs/PBMC_2      Treatment sample
```

### Edit Configuration

Edit `config/config.yaml` to set your analysis parameters:

```yaml
# Choose analysis framework
framework: "seurat"  # or "scanpy"

# Set QC thresholds
qc:
  min_genes: 200
  max_genes: 5000
  max_mito_percent: 20

# Configure clustering
clustering:
  resolution: [0.4, 0.6, 0.8]
  algorithm: "leiden"
```

## Step 3: Validate Configuration

Run a dry run to check the workflow:

```bash
snakemake --dry-run --cores 1
```

This will show you all the jobs that will be executed without actually running them.

## Step 4: Run the Pipeline

### Local Execution

Run on your local machine:

```bash
# Run with 8 cores
snakemake --cores 8 --use-conda

# Monitor progress
tail -f logs/qc/qc_seurat.log
```

### HPC Cluster Execution

Run on a SLURM cluster:

```bash
snakemake --profile config/slurm --jobs 100
```

## Step 5: Monitor Progress

The pipeline will create several output directories:

```bash
results/
├── qc/                 # Quality control results
├── processed/          # Normalized data
├── clustering/         # Clustering results
├── annotation/         # Cell type annotations
├── differential/       # DE analysis
└── reports/           # Final HTML report
```

Check the main log files:

```bash
# View QC log
cat logs/qc/qc_seurat.log

# View clustering log
cat logs/clustering/clustering_seurat.log
```

## Step 6: View Results

### Quality Control Results

```bash
# View QC metrics
head results/qc/qc_metrics.csv

# Open QC plots
open results/qc/qc_plots.pdf
```

### Clustering Results

```bash
# View cluster markers
head results/clustering/cluster_markers.csv

# Open UMAP plot
open results/clustering/umap_clusters.pdf
```

### Final Report

```bash
# Open HTML report
open results/reports/analysis_report.html
```

## Step 7: Interactive Exploration

Launch the Streamlit dashboard for interactive exploration:

```bash
streamlit run dashboard/app.py
```

Then open your browser to `http://localhost:8501`

## Common Workflows

### QC Only

To run only quality control:

```bash
snakemake qc_only --cores 4 --use-conda
```

### Up to Clustering

To run up to clustering analysis:

```bash
snakemake clustering_only --cores 8 --use-conda
```

### Generate Report

To generate only the HTML report:

```bash
snakemake results/reports/analysis_report.html --cores 1
```

## Troubleshooting

### Pipeline Fails

Check the log files in `logs/` directory:

```bash
ls -lh logs/
cat logs/qc/qc_seurat.log
```

### Out of Memory

Reduce the number of parallel jobs:

```bash
snakemake --cores 4 --use-conda  # Instead of 8
```

### Conda Environment Issues

Recreate the environment:

```bash
conda env remove -n scrnaseq
conda env create -f environment.yaml
```

## Next Steps

- Explore [Advanced Usage](../README.md#advanced-usage) for customization
- Read the [User Guide](user_guide/qc.md) for detailed explanations
- Try the [PBMC 3k Tutorial](tutorials/pbmc3k.md) with example data

## Getting Help

- Check the [FAQ](faq.md)
- Open an issue on [GitHub](https://github.com/ishaq7334/scrnaseq_pipeline/issues)
- Join the discussion on [GitHub Discussions](https://github.com/ishaq7334/scrnaseq_pipeline/discussions)
