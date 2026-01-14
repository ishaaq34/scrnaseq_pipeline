#!/usr/bin/env Rscript

# Quality Control Analysis with Seurat
# =====================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

# Parse snakemake inputs
samples_file <- snakemake@params[["samples_file"]]
min_genes <- snakemake@params[["min_genes"]]
max_genes <- snakemake@params[["max_genes"]]
max_mito <- snakemake@params[["max_mito"]]
min_cells <- snakemake@params[["min_cells"]]

# Output files
qc_metrics_file <- snakemake@output[["metrics"]]
qc_plots_file <- snakemake@output[["plots"]]
seurat_obj_file <- snakemake@output[["seurat_obj"]]

# Load sample metadata
samples_df <- read.table(samples_file, header = TRUE, sep = "\t")

# Initialize list to store Seurat objects
seurat_list <- list()

# Load data for each sample
for (i in 1:nrow(samples_df)) {
  sample_id <- samples_df$sample_id[i]
  condition <- samples_df$condition[i]
  batch <- samples_df$batch[i]
  
  cat(sprintf("Loading sample: %s\n", sample_id))
  
  # Load 10x data
  h5_file <- sprintf("results/cellranger/%s/outs/filtered_feature_bc_matrix.h5", sample_id)
  
  if (file.exists(h5_file)) {
    data <- Read10X_h5(h5_file)
  } else {
    # Alternative: load from matrix directory
    data_dir <- sprintf("results/cellranger/%s/outs/filtered_feature_bc_matrix", sample_id)
    data <- Read10X(data_dir)
  }
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = data,
    project = sample_id,
    min.cells = min_cells,
    min.features = min_genes
  )
  
  # Add metadata
  seurat_obj$sample <- sample_id
  seurat_obj$condition <- condition
  seurat_obj$batch <- batch
  
  # Calculate mitochondrial percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Calculate ribosomal percentage
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  
  seurat_list[[sample_id]] <- seurat_obj
}

# Merge all samples
if (length(seurat_list) > 1) {
  seurat_combined <- merge(
    seurat_list[[1]],
    y = seurat_list[2:length(seurat_list)],
    add.cell.ids = names(seurat_list),
    project = "scRNAseq"
  )
} else {
  seurat_combined <- seurat_list[[1]]
}

# Extract QC metrics before filtering
qc_metrics <- data.frame(
  cell_id = colnames(seurat_combined),
  sample = seurat_combined$sample,
  condition = seurat_combined$condition,
  batch = seurat_combined$batch,
  nCount_RNA = seurat_combined$nCount_RNA,
  nFeature_RNA = seurat_combined$nFeature_RNA,
  percent_mt = seurat_combined$percent.mt,
  percent_ribo = seurat_combined$percent.ribo
)

# Save pre-filtering metrics
write.csv(qc_metrics, qc_metrics_file, row.names = FALSE)

# Generate QC plots
pdf(qc_plots_file, width = 12, height = 10)

# Violin plots
p1 <- VlnPlot(seurat_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3, pt.size = 0.1)
print(p1)

# Scatter plots
p2 <- FeatureScatter(seurat_combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
p3 <- FeatureScatter(seurat_combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(p2 + p3)

# Histograms
p4 <- ggplot(qc_metrics, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  geom_vline(xintercept = c(min_genes, max_genes), linetype = "dashed", color = "red") +
  theme_classic() +
  labs(title = "Distribution of Genes per Cell", x = "Number of Genes", y = "Count")

p5 <- ggplot(qc_metrics, aes(x = percent_mt)) +
  geom_histogram(bins = 50, fill = "coral", color = "black") +
  geom_vline(xintercept = max_mito, linetype = "dashed", color = "red") +
  theme_classic() +
  labs(title = "Distribution of Mitochondrial Content", x = "Percent Mitochondrial", y = "Count")

print(p4 / p5)

dev.off()

# Apply QC filters
seurat_filtered <- subset(
  seurat_combined,
  subset = nFeature_RNA > min_genes &
           nFeature_RNA < max_genes &
           percent.mt < max_mito
)

# Print filtering summary
cat("\n=== QC Filtering Summary ===\n")
cat(sprintf("Cells before filtering: %d\n", ncol(seurat_combined)))
cat(sprintf("Cells after filtering: %d\n", ncol(seurat_filtered)))
cat(sprintf("Cells removed: %d (%.1f%%)\n",
            ncol(seurat_combined) - ncol(seurat_filtered),
            100 * (ncol(seurat_combined) - ncol(seurat_filtered)) / ncol(seurat_combined)))

# Save filtered Seurat object
saveRDS(seurat_filtered, file = seurat_obj_file)

cat("\nQC analysis complete!\n")
