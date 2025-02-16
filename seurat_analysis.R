# seurat_analysis.R
# Single-Cell RNA-seq Analysis using Seurat

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Load dataset
DATASET_PATH <- "data/sample_10x/"
seurat_obj <- Read10X(DATASET_PATH)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = "SingleCell", min.cells = 3, min.features = 200)

# Quality Control Filtering
seurat_obj["percent.mt"] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify Highly Variable Features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scaling the Data
seurat_obj <- ScaleData(seurat_obj)

# PCA for Dimensionality Reduction
seurat_obj <- RunPCA(seurat_obj, npcs = 30)

# Find Neighbors and Clusters
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# UMAP for Visualization
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Plot UMAP
pdf("figures/umap_clusters.pdf")
DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("UMAP Clustering")
dev.off()

# Find Marker Genes
seurat_obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(seurat_obj.markers, "results/marker_genes.csv", row.names = FALSE)

# Plot Top Marker Genes
top_markers <- seurat_obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("figures/top_markers.pdf")
DoHeatmap(seurat_obj, features = top_markers$gene) + ggtitle("Top Marker Genes")
dev.off()

# Save Seurat object
saveRDS(seurat_obj, file = "results/seurat_obj.rds")

print("Seurat pipeline completed successfully.")
