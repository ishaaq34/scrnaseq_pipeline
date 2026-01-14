# Single-Cell RNA-seq Analysis Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0.0-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive, production-ready single-cell RNA-seq analysis pipeline built with Snakemake, featuring automated quality control, cell type annotation, differential expression analysis, and interactive visualizations.

## Features

### Complete Analysis Workflow

- **Pre-processing**: Cell Ranger alignment, UMI counting, and matrix generation
- **Quality Control**: Automated filtering, doublet detection, and QC metrics
- **Normalization**: SCTransform, log-normalization, and batch correction
- **Dimensionality Reduction**: PCA, UMAP, t-SNE
- **Clustering**: Leiden, Louvain algorithms with resolution optimization
- **Cell Type Annotation**: Automated annotation using reference datasets
- **Differential Expression**: DESeq2, MAST, Wilcoxon rank-sum tests
- **Pathway Analysis**: GO, KEGG, Reactome enrichment
- **Trajectory Analysis**: RNA velocity, pseudotime inference

### Visualization

- Publication-ready plots (PDF/PNG)
- QC metrics and statistics (CSV)
- Standard scRNA-seq visualizations (UMAP, heatmaps, violin plots)

### Dual Framework Support

- **Seurat (R)**: Industry-standard scRNA-seq analysis
- **Scanpy (Python)**: Fast, scalable analysis for large datasets

### Reproducibility and Scalability

- **Conda/Mamba**: Isolated environments for each tool
- **Docker Support**: Containerized execution
- **HPC Ready**: SLURM, SGE, PBS job submission
- **Cloud Compatible**: AWS, GCP, Azure execution

---

## Repository Structure

```text
scrnaseq_pipeline/
├── workflow/
│   ├── Snakefile                    # Main pipeline orchestration
│   ├── rules/                       # Modular Snakemake rules
│   │   ├── cellranger.smk          # Alignment & counting
│   │   ├── qc.smk                  # Quality control
│   │   ├── normalization.smk       # Data normalization
│   │   ├── clustering.smk          # Clustering & UMAP
│   │   ├── annotation.smk          # Cell type annotation
│   │   ├── differential.smk        # DE analysis
│   │   └── visualization.smk       # Plotting & reports
│   ├── scripts/                     # Analysis scripts
│   │   ├── seurat/                 # R/Seurat scripts
│   │   └── scanpy/                 # Python/Scanpy scripts
│   └── envs/                       # Conda environments
│       ├── cellranger.yaml
│       ├── seurat.yaml
│       └── scanpy.yaml
├── config/
│   ├── config.yaml                 # Main configuration
│   ├── samples.tsv                 # Sample metadata
│   └── references.yaml             # Reference datasets
├── resources/                      # Reference data
│   ├── genomes/                   # Genome references
│   └── annotations/               # Cell type markers
├── results/                       # Pipeline outputs
│   ├── cellranger/               # Alignment results
│   ├── qc/                       # QC metrics & plots
│   ├── processed/                # Normalized data
│   ├── clustering/               # Clusters & UMAPs
│   ├── annotation/               # Cell type labels
│   ├── differential/             # DE results
│   └── reports/                  # HTML reports
├── dashboard/                     # Streamlit dashboard
│   ├── app.py                    # Main dashboard
│   └── components/               # UI components
├── docs/                         # MkDocs documentation
│   ├── index.md
│   ├── installation.md
│   ├── quickstart.md
│   ├── tutorials/
│   └── api/
├── tests/                        # Unit & integration tests
│   ├── test_qc.py
│   └── test_clustering.py
├── examples/                     # Example datasets
│   └── pbmc3k/                  # 3k PBMCs tutorial
├── .github/
│   └── workflows/
│       └── ci.yml               # Continuous integration
├── environment.yaml             # Main conda environment
├── Dockerfile                   # Docker container
├── LICENSE
└── README.md
```

---

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/ishaq7334/scrnaseq_pipeline.git
cd scrnaseq_pipeline

# Create conda environment
conda env create -f environment.yaml
conda activate scrnaseq

# Install Cell Ranger (if using 10x data)
# Download from 10x Genomics website
```

### Configure Your Analysis

Edit `config/config.yaml`:

```yaml
# Analysis framework
framework: "seurat"  # or "scanpy"

# Sample information
samples: "config/samples.tsv"

# Reference genome
genome: "GRCh38"
transcriptome: "resources/genomes/GRCh38"

# QC thresholds
qc:
  min_genes: 200
  max_genes: 5000
  max_mito_percent: 20
  min_cells: 3

# Clustering parameters
clustering:
  resolution: [0.4, 0.6, 0.8, 1.0]
  algorithm: "leiden"

# Differential expression
de:
  method: "wilcox"  # wilcox, MAST, DESeq2
  logfc_threshold: 0.25
  min_pct: 0.1
```

Edit `config/samples.tsv`:

```tsv
sample_id condition batch fastq_dir
PBMC_1 control batch1 data/fastqs/PBMC_1
PBMC_2 treated batch1 data/fastqs/PBMC_2
PBMC_3 control batch2 data/fastqs/PBMC_3
```

### Run the Pipeline

```bash
# Dry run to check workflow
snakemake --dry-run --cores 1

# Run locally (8 cores)
snakemake --cores 8 --use-conda

# Run on HPC cluster (SLURM)
snakemake --profile config/slurm --jobs 100
```

### View Results

```bash
# View QC plots
open results/qc/qc_plots.pdf

# View UMAP clustering
open results/clustering/umap_clusters.pdf
```

---

## Analysis Outputs

### Quality Control

- `results/qc/qc_metrics.csv` - Per-cell QC statistics
- `results/qc/violin_plots.pdf` - QC metric distributions
- `results/qc/scatter_plots.pdf` - Gene vs UMI counts

### Clustering and Visualization

- `results/clustering/umap_clusters.pdf` - UMAP colored by cluster
- `results/clustering/umap_conditions.pdf` - UMAP colored by condition
- `results/clustering/cluster_markers.csv` - Top markers per cluster

### Cell Type Annotation

- `results/annotation/cell_types.csv` - Predicted cell types
- `results/annotation/annotation_umap.pdf` - UMAP with cell type labels
- `results/annotation/marker_heatmap.pdf` - Canonical marker expression

### Differential Expression

- `results/differential/{comparison}/de_genes.csv` - DE results
- `results/differential/{comparison}/volcano_plot.pdf` - Volcano plots
- `results/differential/{comparison}/pathway_enrichment.csv` - GO/KEGG

---

## Tutorials

### Tutorial 1: PBMC 3k Analysis

Analyze 3,000 peripheral blood mononuclear cells from 10x Genomics.

```bash
cd examples/pbmc3k
snakemake --cores 8 --use-conda
```

Expected Results:

- ~2,700 cells after QC
- 8-10 cell type clusters (T cells, B cells, monocytes, NK cells, etc.)
- ~15 minutes runtime on 8 cores

### Tutorial 2: Batch Correction

Integrate multiple samples with Harmony or Seurat integration.

```bash
snakemake --cores 8 --use-conda --config batch_correction=harmony
```

### Tutorial 3: Trajectory Analysis

Perform pseudotime analysis with Monocle3 or scVelo.

```bash
snakemake trajectory --cores 8 --use-conda
```

---

## Example Datasets

The repository includes ready-to-use example datasets for testing and learning:

### PBMC 3k Dataset

**Source**: 10x Genomics Public Datasets  
**URL**: <https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k>  
**Description**: 2,700 Peripheral Blood Mononuclear Cells from a Healthy Donor  
**Size**: ~20 MB (compressed)

**Quick Start**:

```bash
cd examples/pbmc3k
python download_data.py  # Downloads from 10x Genomics
bash run_analysis.sh     # Runs complete pipeline
```

**What's Included**:

- Automated download script
- Pre-configured analysis parameters
- Expected results documentation
- Complete tutorial walkthrough

See `examples/README.md` for more datasets and detailed instructions.

---

---

## Advanced Usage

### Custom Cell Type Markers

Create `resources/annotations/custom_markers.yaml`:

```yaml
T_cells:
  - CD3D
  - CD3E
  - CD8A
B_cells:
  - CD19
  - MS4A1
Monocytes:
  - CD14
  - FCGR3A
NK_cells:
  - GNLY
  - NKG7
```

### Multi-Sample Integration

```yaml
# config/config.yaml
integration:
  enabled: true
  method: "harmony"  # harmony, seurat, scvi
  batch_key: "batch"
```

### Custom Differential Expression

```yaml
# config/config.yaml
comparisons:
  - name: "treated_vs_control"
    group1: ["PBMC_2"]
    group2: ["PBMC_1", "PBMC_3"]
    cell_type: "T_cells"  # Optional: restrict to cell type
```

---

## Documentation

Full documentation available at: **[https://ishaq7334.github.io/scrnaseq_pipeline](https://ishaq7334.github.io/scrnaseq_pipeline)**

- **[Installation Guide](docs/installation.md)** - Detailed setup instructions
- **[User Manual](docs/manual.md)** - Complete pipeline reference
- **[Tutorials](docs/tutorials/)** - Step-by-step walkthroughs
- **[API Reference](docs/api/)** - Script documentation
- **[FAQ](docs/faq.md)** - Common questions & troubleshooting

---

## Testing

```bash
# Run unit tests
pytest tests/

# Test with example data
snakemake --cores 2 --use-conda --directory examples/pbmc3k
```

---

### Development Setup

```bash
# Install development dependencies
conda env create -f environment.dev.yaml
conda activate scrnaseq-dev

# Install pre-commit hooks
pre-commit install

# Run linters
snakemake --lint
black scripts/
flake8 scripts/
```

---

Key Dependencies:

- Cell Ranger (10x Genomics)
- Seurat (Hao et al., Cell 2021)
- Scanpy (Wolf et al., Genome Biology 2018)
- Snakemake (Mölder et al., F1000Research 2021)

---

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

---

---

## Acknowledgments

- **10x Genomics** - Cell Ranger software and example datasets
- **Satija Lab** - Seurat R package
- **Theis Lab** - Scanpy Python package
- **Snakemake Community** - Workflow management framework

---
