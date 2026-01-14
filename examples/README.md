# Example Datasets

This directory contains example datasets and tutorials for the scRNA-seq analysis pipeline.

## Available Examples

### 1. PBMC 3k (Recommended for First-Time Users)

**Directory**: `pbmc3k/`  
**Dataset**: 2,700 Peripheral Blood Mononuclear Cells from a Healthy Donor  
**Source**: 10x Genomics  
**Size**: ~20 MB (compressed)  
**Runtime**: ~20 minutes (8 cores)  
**Difficulty**: Beginner

**Description**: Classic scRNA-seq dataset with well-defined cell populations. Perfect for learning the pipeline and validating installation.

**Quick Start**:

```bash
cd pbmc3k
python download_data.py  # or bash download_data.sh
bash run_analysis.sh
```

**Expected Results**:

- 8-9 distinct cell clusters
- Major immune cell types (T cells, B cells, NK cells, monocytes)
- Clear UMAP separation

**Tutorial**: See `pbmc3k/README.md` for detailed walkthrough

---

### 2. PBMC 10k (Coming Soon)

**Dataset**: 10,000 PBMCs from a Healthy Donor  
**Source**: 10x Genomics  
**Size**: ~50 MB  
**Difficulty**: Intermediate

Larger dataset for testing scalability and performance.

---

### 3. Mouse Brain (Coming Soon)

**Dataset**: 1,300 cells from mouse cortex and hippocampus  
**Source**: 10x Genomics  
**Difficulty**: Intermediate

Example of mouse data analysis with different cell types.

---

### 4. Batch Correction Example (Coming Soon)

**Dataset**: Multiple PBMC samples across batches  
**Source**: Public repositories  
**Difficulty**: Advanced

Demonstrates batch effect correction with Harmony and Seurat integration.

---

## Dataset Sources

All example datasets are from publicly available sources:

1. **10x Genomics**: <https://support.10xgenomics.com/single-cell-gene-expression/datasets>
   - PBMC datasets (3k, 5k, 10k)
   - Mouse brain datasets
   - Cancer cell lines

2. **Human Cell Atlas**: <https://data.humancellatlas.org/>
   - Diverse tissue types
   - Multiple organisms
   - Large-scale datasets

3. **GEO (Gene Expression Omnibus)**: <https://www.ncbi.nlm.nih.gov/geo/>
   - Published datasets
   - Disease vs control comparisons
   - Time-course experiments

4. **Single Cell Portal**: <https://singlecell.broadinstitute.org/>
   - Curated datasets
   - Pre-processed data
   - Interactive exploration

---

## Using Your Own Data

To analyze your own data:

1. **Organize your data**:

   ```
   my_project/
   ├── data/
   │   └── cellranger_output/
   │       └── sample_name/
   │           └── outs/
   │               └── filtered_feature_bc_matrix.h5
   ```

2. **Create configuration**:
   - Copy `pbmc3k/config_pbmc3k.yaml` as template
   - Edit `samples.tsv` with your sample information
   - Adjust QC thresholds based on your data

3. **Run the pipeline**:

   ```bash
   snakemake --cores 8 --use-conda --configfile your_config.yaml
   ```

---

## Data Format Requirements

The pipeline accepts:

### Input Formats

- **10x Genomics outputs**: H5, MEX (MTX + genes + barcodes)
- **H5AD files**: AnnData format (Scanpy)
- **RDS files**: Seurat objects
- **Loom files**: Common scRNA-seq format

### Output Formats

- **H5AD**: For Scanpy and Python tools
- **RDS**: For Seurat and R tools
- **CSV/TSV**: For spreadsheet analysis
- **PDF/PNG**: For publication figures

---

## Dataset Size Guidelines

**Small datasets** (< 5,000 cells):

- RAM: 8-16 GB
- Cores: 2-4
- Runtime: 10-30 minutes

**Medium datasets** (5,000-50,000 cells):

- RAM: 16-32 GB
- Cores: 4-8
- Runtime: 30-120 minutes

**Large datasets** (> 50,000 cells):

- RAM: 32-64 GB
- Cores: 8-16
- Runtime: 2-6 hours
- Consider using HPC cluster

---

## Downloading Example Data

### Method 1: Python Script (Recommended)

```bash
cd pbmc3k
python download_data.py
```

### Method 2: Shell Script

```bash
cd pbmc3k
bash download_data.sh
```

### Method 3: Manual Download

1. Visit: <https://cf.10xgenomics.com/samples/cell/pbmc3k/>
2. Download: `pbmc3k_filtered_gene_bc_matrices.tar.gz`
3. Extract to `data/` directory

---

## Validation

After running an example:

1. **Check cell counts**: Should match expected numbers
2. **Verify clusters**: Compare with published results
3. **Inspect QC metrics**: Ensure reasonable values
4. **Review cell types**: Match known biology

---

## Troubleshooting

### Download fails

- Check internet connection
- Try manual download from 10x website
- Use alternative mirror if available

### Out of memory

- Reduce number of cores
- Close other applications
- Use HPC cluster for large datasets

### Different results

- Random seed variation is normal
- Check parameter settings
- Verify data integrity

---

## Contributing Examples

To add a new example dataset:

1. Create directory: `examples/dataset_name/`
2. Add README with dataset description
3. Include download script
4. Provide optimized configuration
5. Document expected results
6. Submit pull request

---

## Citation

When using these example datasets in publications, please cite:

**PBMC 3k**:

- 10x Genomics (2016). 3k PBMCs from a Healthy Donor. <https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k>

**Pipeline**:

- Your pipeline citation here

---

## Support

For questions about example datasets:

- Open an issue on GitHub
- Check the FAQ
- Join the discussion forum

---

**Note**: All example datasets are provided for educational and research purposes. Please review and comply with the original data licenses and usage terms.
