# PBMC 3k Tutorial

This tutorial demonstrates the complete scRNA-seq analysis pipeline using the publicly available PBMC 3k dataset from 10x Genomics.

## Dataset Information

**Dataset**: Peripheral Blood Mononuclear Cells (PBMCs) from a Healthy Donor  
**Source**: 10x Genomics  
**Technology**: Chromium Single Cell 3' v1  
**Cells**: ~2,700 cells after QC  
**Sequencing**: Illumina NextSeq 500  
**URL**: <https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz>

## Expected Cell Types

This dataset contains several major immune cell populations:

- CD4+ T cells
- CD8+ T cells
- B cells
- NK cells
- Monocytes (CD14+ and FCGR3A+)
- Dendritic cells
- Megakaryocytes

## Tutorial Steps

### Step 1: Download the Dataset

The dataset is automatically downloaded when you run the setup script:

```bash
cd examples/pbmc3k
bash download_data.sh
```

This will download and extract:

- Filtered gene-barcode matrix (H5 format)
- Feature and barcode files
- Pre-computed analysis from 10x

### Step 2: Configure the Analysis

The configuration is already set up in `config_pbmc3k.yaml`:

```yaml
framework: "seurat"  # or "scanpy"

qc:
  min_genes: 200
  max_genes: 2500
  max_mito_percent: 5

clustering:
  resolution: [0.5, 0.8, 1.0]
  algorithm: "leiden"

annotation:
  method: "singleR"
  reference: "celldex::HumanPrimaryCellAtlasData"
```

### Step 3: Run the Analysis

Execute the complete pipeline:

```bash
# From the examples/pbmc3k directory
snakemake --cores 4 --use-conda --configfile config_pbmc3k.yaml

# Or from the main repository directory
snakemake --cores 4 --use-conda --directory examples/pbmc3k
```

### Step 4: Expected Results

**Quality Control**:

- Starting cells: ~2,700
- Cells after filtering: ~2,700 (high-quality dataset)
- Median genes per cell: ~900
- Median UMI counts: ~2,500

**Clustering**:

- Number of clusters: 8-9 (at resolution 0.8)
- Clear separation of major cell types
- Well-defined UMAP structure

**Cell Type Annotation**:

- CD4+ T cells: ~40%
- CD8+ T cells: ~15%
- B cells: ~15%
- NK cells: ~10%
- Monocytes: ~15%
- Dendritic cells: ~3%
- Megakaryocytes: ~2%

### Step 5: Explore Results

Results will be saved in `examples/pbmc3k/results/`:

```bash
# View QC metrics
head results/qc/qc_metrics.csv

# Open UMAP plot
open results/clustering/umap_clusters.pdf

# View cell type annotations
head results/annotation/cell_types.csv

# Open final report
open results/reports/analysis_report.html
```

## Key Marker Genes

**T cells**: CD3D, CD3E, IL7R  
**CD8+ T cells**: CD8A, CD8B  
**CD4+ T cells**: CD4, IL7R  
**B cells**: CD19, MS4A1 (CD20)  
**NK cells**: GNLY, NKG7  
**Monocytes**: CD14, LYZ  
**FCGR3A+ Monocytes**: FCGR3A, MS4A7  
**Dendritic cells**: FCER1A, CST3  
**Megakaryocytes**: PPBP, PF4

## Runtime Expectations

On a standard laptop (8 cores, 16 GB RAM):

- QC: ~2 minutes
- Normalization: ~3 minutes
- Clustering: ~5 minutes
- Annotation: ~3 minutes
- Differential expression: ~5 minutes
- **Total**: ~20 minutes

## Troubleshooting

### Issue: Download fails

**Solution**: Manually download from 10x Genomics website:

```bash
wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
```

### Issue: Out of memory

**Solution**: Reduce number of cores:

```bash
snakemake --cores 2 --use-conda
```

### Issue: Conda environment creation fails

**Solution**: Use mamba for faster resolution:

```bash
conda install -n base -c conda-forge mamba
snakemake --cores 4 --use-conda --conda-frontend mamba
```

## Validation

Compare your results with the published 10x analysis:

- <https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k>

Your results should show similar:

- Number of cells
- Median genes per cell
- Clustering structure
- Cell type proportions

## Next Steps

After completing this tutorial:

1. Try different clustering resolutions
2. Test both Seurat and Scanpy frameworks
3. Modify QC thresholds and observe effects
4. Add custom marker genes for annotation
5. Perform differential expression between cell types

## References

- 10x Genomics PBMC 3k dataset: <https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k>
- Seurat PBMC tutorial: <https://satijalab.org/seurat/articles/pbmc3k_tutorial.html>
- Scanpy PBMC tutorial: <https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html>
