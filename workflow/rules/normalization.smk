"""
Normalization and Scaling Rules
================================

Rules for data normalization, scaling, and feature selection.
"""

if FRAMEWORK == "seurat":
    
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
            n_top_genes = config["normalization"]["n_top_genes"],
            vars_to_regress = config["normalization"].get("sct_vars_to_regress", [])
        threads: 4
        conda:
            "../envs/seurat.yaml"
        log:
            "logs/normalization/normalize_seurat.log"
        script:
            "../scripts/seurat/03_normalize.R"

elif FRAMEWORK == "scanpy":
    
    rule normalize_scanpy:
        """
        Normalize and scale data with Scanpy
        """
        input:
            adata = "results/qc/adata_qc.h5ad"
        output:
            adata = "results/processed/adata.h5ad",
            hvg_plot = "results/processed/highly_variable_genes.pdf"
        params:
            method = config["normalization"]["method"],
            target_sum = config["normalization"]["target_sum"],
            n_top_genes = config["normalization"]["n_top_genes"]
        threads: 4
        conda:
            "../envs/scanpy.yaml"
        log:
            "logs/normalization/normalize_scanpy.log"
        script:
            "../scripts/scanpy/03_normalize.py"


rule batch_correction:
    """
    Perform batch correction if enabled
    """
    input:
        obj = "results/processed/seurat_object.rds" if FRAMEWORK == "seurat" else "results/processed/adata.h5ad"
    output:
        obj = "results/processed/integrated_object.rds" if FRAMEWORK == "seurat" else "results/processed/adata_integrated.h5ad"
    params:
        method = config["integration"]["method"],
        batch_key = config["integration"]["batch_key"]
    conda:
        "../envs/seurat.yaml" if FRAMEWORK == "seurat" else "../envs/scanpy.yaml"
    log:
        "logs/normalization/batch_correction.log"
    script:
        "../scripts/seurat/04_integration.R" if FRAMEWORK == "seurat" else "../scripts/scanpy/04_integration.py"
