"""
Differential Expression Analysis Rules

Rules for comparing gene expression between conditions using Seurat.
"""

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
