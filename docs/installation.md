# Installation Guide

## Prerequisites

### System Requirements

- **Operating System**: Linux, macOS, or Windows (with WSL2)
- **Memory**: Minimum 16 GB RAM (32 GB recommended for large datasets)
- **Storage**: At least 100 GB free disk space
- **CPU**: Multi-core processor (8+ cores recommended)

### Required Software

1. **Conda or Mamba** (package manager)
2. **Git** (version control)
3. **Cell Ranger** (for 10x Genomics data, optional)

---

## Installation Steps

### Step 1: Install Conda/Mamba

If you don't have Conda installed:

```bash
# Download Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install
bash Miniconda3-latest-Linux-x86_64.sh

# Initialize
conda init bash
source ~/.bashrc
```

For faster package installation, install Mamba:

```bash
conda install -n base -c conda-forge mamba
```

### Step 2: Clone the Repository

```bash
git clone https://github.com/username/scrnaseq_pipeline.git
cd scrnaseq_pipeline
```

### Step 3: Create Conda Environment

```bash
# Create main environment
conda env create -f environment.yaml

# Activate environment
conda activate scrnaseq
```

### Step 4: Install Cell Ranger (Optional)

If you're starting from FASTQ files:

1. Download Cell Ranger from [10x Genomics](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)
2. Extract and add to PATH:

```bash
tar -xzvf cellranger-7.2.0.tar.gz
export PATH=/path/to/cellranger-7.2.0:$PATH
```

1. Download reference genome:

```bash
# Human GRCh38
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz
mv refdata-gex-GRCh38-2024-A resources/genomes/
```

### Step 5: Verify Installation

```bash
# Check Snakemake
snakemake --version

# Test pipeline
snakemake --dry-run --cores 1
```

---

## Framework-Specific Setup

### For Seurat (R) Analysis

```bash
conda env create -f workflow/envs/seurat.yaml
conda activate seurat

# Test R installation
R --version
```

### For Scanpy (Python) Analysis

```bash
conda env create -f workflow/envs/scanpy.yaml
conda activate scanpy

# Test Python installation
python -c "import scanpy; print(scanpy.__version__)"
```

---

## HPC Cluster Setup

### SLURM Configuration

Create `config/slurm/config.yaml`:

```yaml
cluster: "sbatch --partition={resources.partition} --cpus-per-task={threads} --mem={resources.mem_gb}G --time={resources.time} --job-name={rule} --output=logs/slurm/{rule}_%j.out"
jobs: 100
```

### Running on SLURM

```bash
snakemake --profile config/slurm
```

---

## Troubleshooting

### Common Issues

**Issue**: Conda environment creation fails

**Solution**: Try using mamba instead:

```bash
mamba env create -f environment.yaml
```

**Issue**: Cell Ranger not found

**Solution**: Ensure Cell Ranger is in your PATH:

```bash
which cellranger
export PATH=/path/to/cellranger:$PATH
```

**Issue**: Out of memory errors

**Solution**: Reduce the number of parallel jobs or increase memory allocation in `config/config.yaml`

---

## Next Steps

After installation, proceed to the [Quick Start Guide](quickstart.md) to run your first analysis.
