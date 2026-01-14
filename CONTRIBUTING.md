# Contributing Guidelines

Thank you for your interest in contributing to the scRNA-seq Analysis Pipeline! This document provides guidelines for contributing to the project.

## Code of Conduct

We are committed to providing a welcoming and inclusive environment. Please be respectful and professional in all interactions.

## How to Contribute

### Reporting Bugs

If you find a bug, please create an issue with:

- A clear, descriptive title
- Steps to reproduce the problem
- Expected vs actual behavior
- Your environment (OS, Conda version, etc.)
- Relevant log files or error messages

### Suggesting Enhancements

Enhancement suggestions are welcome! Please create an issue describing:

- The motivation for the enhancement
- Detailed description of the proposed feature
- Potential implementation approach
- Any alternative solutions considered

### Pull Requests

1. **Fork the repository** and create a new branch:

   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes** following our coding standards

3. **Test your changes**:

   ```bash
   snakemake --dry-run --cores 1
   pytest tests/
   ```

4. **Commit your changes** with clear, descriptive messages:

   ```bash
   git commit -m "Add feature: description of changes"
   ```

5. **Push to your fork** and submit a pull request

## Development Setup

```bash
# Clone your fork
git clone https://github.com/your-ishaq7334/scrnaseq_pipeline.git
cd scrnaseq_pipeline

# Add upstream remote
git remote add upstream https://github.com/original-owner/scrnaseq_pipeline.git

# Create development environment
conda env create -f environment.dev.yaml
conda activate scrnaseq-dev

# Install pre-commit hooks
pre-commit install
```

## Coding Standards

### Python Code

- Follow PEP 8 style guidelines
- Use type hints where appropriate
- Write docstrings for all functions and classes
- Maximum line length: 100 characters

```python
def process_data(adata: sc.AnnData, min_genes: int = 200) -> sc.AnnData:
    """
    Process single-cell data with quality filters.
    
    Parameters
    ----------
    adata : sc.AnnData
        Annotated data matrix
    min_genes : int, default=200
        Minimum number of genes per cell
        
    Returns
    -------
    sc.AnnData
        Filtered annotated data matrix
    """
    return adata
```

### R Code

- Follow tidyverse style guide
- Use roxygen2 documentation
- Consistent naming conventions

```r
#' Process Seurat Object
#'
#' @param seurat_obj Seurat object
#' @param min_genes Minimum genes per cell
#' @return Filtered Seurat object
process_seurat <- function(seurat_obj, min_genes = 200) {
  # Implementation
  return(seurat_obj)
}
```

### Snakemake Rules

- Clear rule names describing the operation
- Comprehensive docstrings
- Proper resource specifications

```python
rule example_rule:
    """
    Brief description of what this rule does
    """
    input:
        "path/to/input.txt"
    output:
        "path/to/output.txt"
    params:
        param1 = config["param1"]
    threads: 4
    resources:
        mem_gb = 16
    conda:
        "../envs/environment.yaml"
    log:
        "logs/example_rule.log"
    script:
        "../scripts/example_script.py"
```

## Testing

### Unit Tests

Add tests for new functionality:

```python
# tests/test_qc.py
import pytest
import scanpy as sc

def test_qc_filtering():
    """Test QC filtering removes low-quality cells"""
    # Test implementation
    pass
```

### Integration Tests

Test complete workflows:

```bash
# Test with example data
snakemake --cores 2 --use-conda --directory examples/pbmc3k
```

## Documentation

- Update relevant documentation files
- Add docstrings to new functions
- Include usage examples
- Update CHANGELOG.md

## Review Process

All pull requests will be reviewed for:

- Code quality and style
- Test coverage
- Documentation completeness
- Compatibility with existing features

## Questions?

Feel free to open an issue for any questions about contributing.

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
