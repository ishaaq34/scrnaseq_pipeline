"""
Differential Expression Analysis Rules
=======================================

Rules for comparing gene expression between conditions.
"""

if FRAMEWORK == "seurat":
    
    rule differential_expression_seurat:
        """
        Perform differential expression analysis with Seurat
        """
        input:
            seurat_obj = "results/annotation/seurat_annotated.rds"
        output:
            de_results = "results/differential/{comparison}/de_genes.csv",
            volcano_plot = "results/differential/{comparison}/volcano_plot.pdf",
            heatmap = "results/differential/{comparison}/de_heatmap.pdf",
            top_genes = "results/differential/{comparison}/top_de_genes.csv"
        params:
            method = config["differential"]["method"],
            logfc_threshold = config["differential"]["logfc_threshold"],
            min_pct = config["differential"]["min_pct"],
            pval_cutoff = config["differential"]["pval_cutoff"],
            comparison = lambda wildcards: wildcards.comparison
        conda:
            "../envs/seurat.yaml"
        log:
            "logs/differential/{comparison}_seurat.log"
        script:
            "../scripts/seurat/09_differential.R"

elif FRAMEWORK == "scanpy":
    
    rule differential_expression_scanpy:
        """
        Perform differential expression analysis with Scanpy
        """
        input:
            adata = "results/annotation/adata_annotated.h5ad"
        output:
            de_results = "results/differential/{comparison}/de_genes.csv",
            volcano_plot = "results/differential/{comparison}/volcano_plot.pdf",
            heatmap = "results/differential/{comparison}/de_heatmap.pdf",
            top_genes = "results/differential/{comparison}/top_de_genes.csv"
        params:
            method = config["differential"]["method"],
            logfc_threshold = config["differential"]["logfc_threshold"],
            pval_cutoff = config["differential"]["pval_cutoff"],
            comparison = lambda wildcards: wildcards.comparison
        conda:
            "../envs/scanpy.yaml"
        log:
            "logs/differential/{comparison}_scanpy.log"
        script:
            "../scripts/scanpy/09_differential.py"


rule pathway_enrichment:
    """
    Perform pathway enrichment analysis on DE genes
    """
    input:
        de_results = "results/differential/{comparison}/de_genes.csv"
    output:
        enrichment = "results/differential/{comparison}/pathway_enrichment.csv",
        dotplot = "results/differential/{comparison}/enrichment_dotplot.pdf",
        network = "results/differential/{comparison}/enrichment_network.pdf"
    params:
        databases = config["enrichment"]["databases"],
        organism = config["enrichment"]["organism"],
        pval_cutoff = config["enrichment"]["pval_cutoff"]
    conda:
        "../envs/enrichment.yaml"
    log:
        "logs/differential/{comparison}_enrichment.log"
    script:
        "../scripts/10_enrichment.R"
