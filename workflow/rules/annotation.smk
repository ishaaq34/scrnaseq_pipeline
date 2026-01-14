"""
Cell Type Annotation Rules
===========================

Rules for automated and manual cell type annotation.
"""

if FRAMEWORK == "seurat":
    
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

elif FRAMEWORK == "scanpy":
    
    rule annotate_celltypes_scanpy:
        """
        Annotate cell types using CellTypist or manual markers
        """
        input:
            adata = "results/clustering/adata_clustered.h5ad",
            markers = config["annotation"].get("custom_markers", "")
        output:
            adata = "results/annotation/adata_annotated.h5ad",
            cell_types = "results/annotation/cell_types.csv",
            umap_plot = "results/annotation/annotation_umap.pdf",
            marker_plot = "results/annotation/marker_dotplot.pdf"
        params:
            method = config["annotation"]["method"],
            min_confidence = config["annotation"]["min_confidence"]
        conda:
            "../envs/scanpy.yaml"
        log:
            "logs/annotation/annotate_scanpy.log"
        script:
            "../scripts/scanpy/07_annotation.py"


rule validate_annotations:
    """
    Validate cell type annotations with canonical markers
    """
    input:
        obj = "results/annotation/seurat_annotated.rds" if FRAMEWORK == "seurat" else "results/annotation/adata_annotated.h5ad"
    output:
        validation_report = "results/annotation/validation_report.html",
        marker_heatmap = "results/annotation/canonical_markers_heatmap.pdf"
    conda:
        "../envs/seurat.yaml" if FRAMEWORK == "seurat" else "../envs/scanpy.yaml"
    log:
        "logs/annotation/validate.log"
    script:
        "../scripts/seurat/08_validate_annotation.R" if FRAMEWORK == "seurat" else "../scripts/scanpy/08_validate_annotation.py"
