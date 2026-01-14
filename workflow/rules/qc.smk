"""
Quality Control Rules

Rules for filtering cells and genes based on QC metrics using Seurat.
"""

rule qc_seurat:
    """
    Perform quality control with Seurat
    """
    input:
        samples_file = config["samples"]
    output:
        metrics = "results/qc/qc_metrics.csv",
        plots = "results/qc/qc_plots.pdf",
        seurat_obj = "results/qc/seurat_qc.rds"
    params:
        min_genes = config["qc"]["min_genes"],
        max_genes = config["qc"]["max_genes"],
        max_mito = config["qc"]["max_mito_percent"],
        min_cells = config["qc"]["min_cells"]
    threads: 4
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/qc/qc_seurat.log"
    script:
        "../scripts/seurat/01_qc.R"
