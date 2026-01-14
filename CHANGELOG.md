# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Initial release of scRNA-seq analysis pipeline
- Snakemake workflow with modular rules
- Support for both Seurat (R) and Scanpy (Python) frameworks
- Cell Ranger integration for 10x Genomics data
- Automated quality control and filtering
- Normalization with multiple methods (log-norm, SCTransform)
- Batch correction (Harmony, Seurat integration)
- Clustering with Leiden and Louvain algorithms
- Cell type annotation (SingleR, CellTypist)
- Differential expression analysis
- Pathway enrichment analysis
- Comprehensive visualization suite
- Interactive Streamlit dashboard
- MkDocs documentation
- Example PBMC 3k tutorial
- Conda environments for reproducibility
- Docker support
- HPC cluster profiles (SLURM, SGE, PBS)
- Continuous integration with GitHub Actions
- Unit and integration tests

### Changed

- N/A (initial release)

### Deprecated

- N/A (initial release)

### Removed

- N/A (initial release)

### Fixed

- N/A (initial release)

### Security

- N/A (initial release)

## [1.0.0] - 2026-01-14

### Added

- First stable release

---

## Release Notes

### Version 1.0.0

This is the first stable release of the scRNA-seq Analysis Pipeline. Key features include:

**Core Functionality:**

- Complete end-to-end scRNA-seq analysis workflow
- Dual framework support (Seurat and Scanpy)
- Automated quality control and cell filtering
- Multiple normalization strategies
- Advanced clustering algorithms
- Reference-based cell type annotation

**Reproducibility:**

- Conda/Mamba environment management
- Docker containerization
- Snakemake workflow management
- Version-controlled analysis parameters

**Scalability:**

- HPC cluster support
- Cloud computing compatibility
- Efficient resource management
- Parallel processing capabilities

**Documentation:**

- Comprehensive user manual
- Step-by-step tutorials
- API reference
- Troubleshooting guide

**Quality Assurance:**

- Unit tests for core functions
- Integration tests with example data
- Continuous integration pipeline
- Code quality checks

For detailed usage instructions, see the [documentation](docs/index.md).
