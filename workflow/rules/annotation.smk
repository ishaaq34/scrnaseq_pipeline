"""
Cell Type Annotation Rules

Rules for automated and manual cell type annotation using Seurat.
"""

rule annotate_celltypes_seurat:
    """
    Annotate cell types using SingleR or manual markers
    """
    input:
        seurat_obj = "results/clustering/seurat_clustered.rds",
        markers = config["annotation"].get("custom_markers", "")
    output:
        seurat_obj = "results/annotation/seurat_annotated.rds",
        cell_types = "results/annotation/cell_types.csv",
        umap_plot = "results/annotation/annotation_umap.pdf",
        marker_plot = "results/annotation/marker_dotplot.pdf"
    params:
        method = config["annotation"]["method"],
        reference = config["annotation"].get("reference", ""),
        min_confidence = config["annotation"]["min_confidence"]
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/annotation/annotate_seurat.log"
    script:
        "../scripts/seurat/07_annotation.R"
