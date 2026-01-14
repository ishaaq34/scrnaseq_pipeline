"""
Visualization and Reporting Rules
==================================

Rules for generating plots and HTML reports.
"""

rule generate_report:
    """
    Generate comprehensive HTML analysis report
    """
    input:
        qc_metrics = "results/qc/qc_metrics.csv",
        cluster_markers = "results/clustering/cluster_markers.csv",
        cell_types = "results/annotation/cell_types.csv",
        de_results = expand("results/differential/{comparison}/de_genes.csv",
                          comparison=[c["name"] for c in config.get("comparisons", [])])
    output:
        report = "results/reports/analysis_report.html"
    params:
        title = config["report"]["title"],
        author = config["report"]["author"]
    conda:
        "../envs/reporting.yaml"
    log:
        "logs/reporting/generate_report.log"
    script:
        "../scripts/11_generate_report.Rmd"


rule plot_umap_features:
    """
    Generate UMAP plots colored by various features
    """
    input:
        obj = "results/annotation/seurat_annotated.rds" if FRAMEWORK == "seurat" else "results/annotation/adata_annotated.h5ad"
    output:
        plots = expand("results/visualization/umap_{feature}.pdf",
                      feature=config["visualization"]["umap_features"])
    params:
        features = config["visualization"]["umap_features"]
    conda:
        "../envs/seurat.yaml" if FRAMEWORK == "seurat" else "../envs/scanpy.yaml"
    log:
        "logs/visualization/umap_features.log"
    script:
        "../scripts/seurat/12_plot_features.R" if FRAMEWORK == "seurat" else "../scripts/scanpy/12_plot_features.py"


rule plot_gene_expression:
    """
    Generate feature plots for specific genes
    """
    input:
        obj = "results/annotation/seurat_annotated.rds" if FRAMEWORK == "seurat" else "results/annotation/adata_annotated.h5ad"
    output:
        feature_plots = "results/visualization/feature_genes.pdf",
        violin_plots = "results/visualization/violin_genes.pdf",
        dotplot = "results/visualization/dotplot_genes.pdf"
    params:
        genes = config["visualization"]["feature_genes"]
    conda:
        "../envs/seurat.yaml" if FRAMEWORK == "seurat" else "../envs/scanpy.yaml"
    log:
        "logs/visualization/gene_expression.log"
    script:
        "../scripts/seurat/13_plot_genes.R" if FRAMEWORK == "seurat" else "../scripts/scanpy/13_plot_genes.py"


rule create_interactive_viewer:
    """
    Create interactive cellxgene viewer
    """
    input:
        adata = "results/annotation/adata_annotated.h5ad"
    output:
        viewer = directory("results/interactive/cellxgene")
    conda:
        "../envs/scanpy.yaml"
    log:
        "logs/visualization/cellxgene.log"
    shell:
        """
        mkdir -p {output.viewer}
        cp {input.adata} {output.viewer}/data.h5ad
        echo "Run: cellxgene launch {output.viewer}/data.h5ad" > {output.viewer}/README.txt
        """
