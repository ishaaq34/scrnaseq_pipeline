"""
Clustering and Dimensionality Reduction Rules

Rules for PCA, UMAP, and clustering using Seurat.
"""

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
    threads: 8
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/clustering/clustering_seurat.log"
    script:
        "../scripts/seurat/05_clustering.R"
