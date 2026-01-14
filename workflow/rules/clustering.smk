"""
Clustering and Dimensionality Reduction Rules
==============================================

Rules for PCA, UMAP, t-SNE, and clustering.
"""

if FRAMEWORK == "seurat":
    
    rule clustering_seurat:
        """
        Perform dimensionality reduction and clustering with Seurat
        """
        input:
            seurat_obj = "results/processed/seurat_object.rds"
        output:
            seurat_obj = "results/clustering/seurat_clustered.rds",
            umap_plot = "results/clustering/umap_clusters.pdf",
            pca_plot = "results/clustering/pca_plot.pdf",
            elbow_plot = "results/clustering/elbow_plot.pdf",
            markers = "results/clustering/cluster_markers.csv"
        params:
            n_pcs = config["clustering"]["n_pcs"],
            resolution = config["clustering"]["resolution"],
            algorithm = config["clustering"]["algorithm"]
        threads: config["resources"]["clustering"]["threads"]
        resources:
            mem_gb = config["resources"]["clustering"]["memory_gb"]
        conda:
            "../envs/seurat.yaml"
        log:
            "logs/clustering/clustering_seurat.log"
        script:
            "../scripts/seurat/05_clustering.R"

elif FRAMEWORK == "scanpy":
    
    rule clustering_scanpy:
        """
        Perform dimensionality reduction and clustering with Scanpy
        """
        input:
            adata = "results/processed/adata.h5ad"
        output:
            adata = "results/clustering/adata_clustered.h5ad",
            umap_plot = "results/clustering/umap_clusters.pdf",
            pca_plot = "results/clustering/pca_plot.pdf",
            markers = "results/clustering/cluster_markers.csv"
        params:
            n_pcs = config["clustering"]["n_pcs"],
            n_neighbors = config["clustering"]["n_neighbors"],
            resolution = config["clustering"]["resolution"]
        threads: config["resources"]["clustering"]["threads"]
        resources:
            mem_gb = config["resources"]["clustering"]["memory_gb"]
        conda:
            "../envs/scanpy.yaml"
        log:
            "logs/clustering/clustering_scanpy.log"
        script:
            "../scripts/scanpy/05_clustering.py"


rule find_cluster_markers:
    """
    Find marker genes for each cluster
    """
    input:
        obj = "results/clustering/seurat_clustered.rds" if FRAMEWORK == "seurat" else "results/clustering/adata_clustered.h5ad"
    output:
        markers = "results/clustering/all_markers.csv",
        top_markers = "results/clustering/top_markers.csv",
        heatmap = "results/clustering/marker_heatmap.pdf"
    params:
        method = config["markers"]["method"],
        min_pct = config["markers"]["min_pct"],
        logfc_threshold = config["markers"]["logfc_threshold"],
        top_n = config["markers"]["top_n"]
    conda:
        "../envs/seurat.yaml" if FRAMEWORK == "seurat" else "../envs/scanpy.yaml"
    log:
        "logs/clustering/find_markers.log"
    script:
        "../scripts/seurat/06_markers.R" if FRAMEWORK == "seurat" else "../scripts/scanpy/06_markers.py"
