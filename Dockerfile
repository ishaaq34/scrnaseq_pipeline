FROM continuumio/miniconda3:latest

LABEL maintainer="your.email@example.com"
LABEL description="Docker container for scRNA-seq analysis pipeline"
LABEL version="1.0.0"

# Set working directory
WORKDIR /pipeline

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    wget \
    curl \
    vim \
    && rm -rf /var/lib/apt/lists/*

# Copy environment files
COPY environment.yaml /pipeline/
COPY workflow/envs/ /pipeline/workflow/envs/

# Create conda environments
RUN conda env create -f environment.yaml && \
    conda clean -afy

# Install mamba for faster package resolution
RUN conda install -n base -c conda-forge mamba && \
    conda clean -afy

# Copy pipeline files
COPY . /pipeline/

# Activate conda environment
SHELL ["conda", "run", "-n", "scrnaseq", "/bin/bash", "-c"]

# Set PATH
ENV PATH /opt/conda/envs/scrnaseq/bin:$PATH

# Create necessary directories
RUN mkdir -p /pipeline/results /pipeline/logs /pipeline/resources

# Set entrypoint
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "scrnaseq"]
CMD ["snakemake", "--help"]

# Usage:
# Build: docker build -t scrnaseq-pipeline .
# Run: docker run -v $(pwd)/data:/pipeline/data -v $(pwd)/results:/pipeline/results scrnaseq-pipeline snakemake --cores 8
