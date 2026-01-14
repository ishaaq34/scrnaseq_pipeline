"""
Quality Control Rules
======================

Rules for filtering cells and genes based on QC metrics.
"""

if FRAMEWORK == "seurat":
    
    rule qc_seurat:
        """
        Perform quality control with Seurat
        """
        input:
            h5_files = expand("results/cellranger/{sample}/outs/filtered_feature_bc_matrix.h5",
                            sample=SAMPLES) if config.get("run_cellranger", False) else []
        output:
            metrics = "results/qc/qc_metrics.csv",
            plots = "results/qc/qc_plots.pdf",
            seurat_obj = "results/qc/seurat_qc.rds"
        params:
            min_genes = config["qc"]["min_genes"],
            max_genes = config["qc"]["max_genes"],
            max_mito = config["qc"]["max_mito_percent"],
            min_cells = config["qc"]["min_cells"],
            samples_file = config["samples"]
        threads: config["resources"]["qc"]["threads"]
        resources:
            mem_gb = config["resources"]["qc"]["memory_gb"]
        conda:
            "../envs/seurat.yaml"
        log:
            "logs/qc/qc_seurat.log"
        script:
            "../scripts/seurat/01_qc.R"

elif FRAMEWORK == "scanpy":
    
    rule qc_scanpy:
        """
        Perform quality control with Scanpy
        """
        input:
            h5_files = expand("results/cellranger/{sample}/outs/filtered_feature_bc_matrix.h5",
                            sample=SAMPLES) if config.get("run_cellranger", False) else []
        output:
            metrics = "results/qc/qc_metrics.csv",
            plots = "results/qc/qc_plots.pdf",
            adata = "results/qc/adata_qc.h5ad"
        params:
            min_genes = config["qc"]["min_genes"],
            max_genes = config["qc"]["max_genes"],
            max_mito = config["qc"]["max_mito_percent"],
            min_cells = config["qc"]["min_cells"],
            samples_file = config["samples"]
        threads: config["resources"]["qc"]["threads"]
        resources:
            mem_gb = config["resources"]["qc"]["memory_gb"]
        conda:
            "../envs/scanpy.yaml"
        log:
            "logs/qc/qc_scanpy.log"
        script:
            "../scripts/scanpy/01_qc.py"


rule detect_doublets:
    """
    Detect doublets using DoubletFinder (Seurat) or Scrublet (Scanpy)
    """
    input:
        seurat_obj = "results/qc/seurat_qc.rds" if FRAMEWORK == "seurat" else "results/qc/adata_qc.h5ad"
    output:
        doublet_scores = "results/qc/doublet_scores.csv",
        doublet_plot = "results/qc/doublet_plot.pdf"
    params:
        expected_rate = config["qc"]["expected_doublet_rate"]
    conda:
        "../envs/seurat.yaml" if FRAMEWORK == "seurat" else "../envs/scanpy.yaml"
    log:
        "logs/qc/detect_doublets.log"
    script:
        "../scripts/seurat/02_doublets.R" if FRAMEWORK == "seurat" else "../scripts/scanpy/02_doublets.py"
