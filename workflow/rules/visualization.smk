"""
Visualization Rules

Rules for generating plots using Seurat.
"""

rule plot_umap_features:
    """
    Generate UMAP plots colored by various features
    """
    input:
        obj = "results/annotation/seurat_annotated.rds"
    output:
        plots = expand("results/visualization/umap_{feature}.pdf",
                      feature=config["visualization"]["umap_features"])
    params:
        features = config["visualization"]["umap_features"]
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/visualization/umap_features.log"
    script:
        "../scripts/seurat/12_plot_features.R"


rule plot_gene_expression:
    """
    Generate feature plots for specific genes
    """
    input:
        obj = "results/annotation/seurat_annotated.rds"
    output:
        feature_plots = "results/visualization/feature_genes.pdf",
        violin_plots = "results/visualization/violin_genes.pdf",
        dotplot = "results/visualization/dotplot_genes.pdf"
    params:
        genes = config["visualization"]["feature_genes"]
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/visualization/gene_expression.log"
    script:
        "../scripts/seurat/13_plot_genes.R"
