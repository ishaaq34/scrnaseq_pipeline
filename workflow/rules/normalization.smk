"""
Normalization and Scaling Rules

Rules for data normalization, scaling, and feature selection using Seurat.
"""

rule normalize_seurat:
    """
    Normalize and scale data with Seurat
    """
    input:
        seurat_obj = "results/qc/seurat_qc.rds"
    output:
        seurat_obj = "results/processed/seurat_object.rds",
        hvg_plot = "results/processed/highly_variable_genes.pdf"
    params:
        method = config["normalization"]["method"],
        n_top_genes = config["normalization"]["n_top_genes"]
    threads: 4
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/normalization/normalize_seurat.log"
    script:
        "../scripts/seurat/03_normalize.R"
