# scRNA-seq Analysis Pipeline

Welcome to the comprehensive single-cell RNA-seq analysis pipeline documentation.

## Overview

This pipeline provides an end-to-end solution for analyzing single-cell RNA-sequencing data, from raw FASTQ files to publication-ready figures and differential expression results.

## Key Features

- **Automated Workflow**: Snakemake-based pipeline handles all analysis steps
- **Dual Framework Support**: Choose between Seurat (R) or Scanpy (Python)
- **Quality Control**: Comprehensive filtering and doublet detection
- **Cell Type Annotation**: Automated annotation using reference datasets
- **Differential Expression**: Multiple statistical methods available
- **Visualization**: Publication-ready plots and interactive dashboards
- **Reproducibility**: Conda environments and Docker support
- **Scalability**: HPC cluster and cloud computing ready

## Workflow Overview

```text
FASTQ Files
    ↓
Cell Ranger Alignment
    ↓
Quality Control
    ↓
Normalization & Scaling
    ↓
Dimensionality Reduction (PCA, UMAP)
    ↓
Clustering
    ↓
Cell Type Annotation
    ↓
Differential Expression
    ↓
Pathway Enrichment
    ↓
Reports & Visualizations
```

## Quick Links

- [Installation Guide](installation.md) - Get started with installation
- [Quick Start](quickstart.md) - Run your first analysis
- [Configuration](configuration.md) - Customize pipeline parameters
- [Tutorials](tutorials/pbmc3k.md) - Step-by-step examples
- [API Reference](api/rules.md) - Detailed documentation

## Supported Data Types

- 10x Genomics Chromium (3' and 5')
- Drop-seq
- inDrop
- Smart-seq2
- Any UMI-based scRNA-seq protocol

## Supported Organisms

- Human (Homo sapiens)
- Mouse (Mus musculus)
- Custom genomes (with appropriate reference)

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{scrnaseq_pipeline,
  author = {RAJ ISHAQ NABI KHAN},
  title = {scRNA-seq Analysis Pipeline},
  year = {2026},
  url = {https://github.com/username/scrnaseq_pipeline}
}
```

## Getting Help

- **Issues**: Report bugs on [GitHub Issues](https://github.com/username/scrnaseq_pipeline/issues)
- **Discussions**: Ask questions on [GitHub Discussions](https://github.com/username/scrnaseq_pipeline/discussions)
- **Email**: Contact the maintainers at <contact@example.com>

## License

This project is licensed under the MIT License - see the [LICENSE](../LICENSE) file for details.
