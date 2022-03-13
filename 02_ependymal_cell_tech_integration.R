# Analysis of the Smart-seq2 spinal cord ependymal cell dataset from FOXj1-EGFP mice
# and integration of the Smart-seq2 and 10X ependymal cell datasets to study spinal cord ependymal cell heterogeneity

library(Seurat)
library(tidyverse)
library(viridis)

# Load the data
smartseq2 <- Read10X(data.dir = "~/path/to/STARsolo/Solo.out/Gene/filtered/")
smartseq2[10:15, 10:15]

colnames(smartseq2)
# Clean up column (cell) names to implement in the combined count matrix
colnames(smartseq2) <- gsub("-1_S(\\d+)", "", colnames(smartseq2))
colnames(smartseq2)
colnames(smartseq2) <- gsub("-", ".", colnames(smartseq2))
colnames(smartseq2)

smartseq2[1:5,1:3]

cellnames <- colnames(smartseq2)

# Create a Seurat object
smartseq2 <- CreateSeuratObject(counts = smartseq2, min.cells = 3, min.features = 500, project = "Smart-Seq2")
smartseq2

# View metadata (stored in object@meta.data)
smartseq2[[]]
names(smartseq2[[]])

table(smartseq2[["orig.ident"]]) # OK!

# Add the mouse, age, and tech column to metadata slot
# Extract plate name as mouse ID
#cellnames %>% str_split(., '\\.', simplify = T) %>% .[,2] %>% print()
mouse <- cellnames %>% str_extract(., 'B\\d{6}')
age <- rep("young", 915) # 915 cells
tech <- rep("Smart-Seq2", 915)
smartseq2 <- AddMetaData(smartseq2, metadata = age, col.name = "age")
smartseq2 <- AddMetaData(smartseq2, metadata = tech, col.name = "tech")
smartseq2 <- AddMetaData(smartseq2, metadata = mouse, col.name = "mouse")

names(smartseq2[[]])
table(smartseq2[["mouse"]])


## Data QC
VlnPlot(smartseq2, c("nFeature_RNA", "nCount_RNA"), pt.size = 0)
median(smartseq2$nFeature_RNA)
median(smartseq2$nCount_RNA)

FeatureScatter(smartseq2, feature1 = "nFeature_RNA", feature2 = 'nCount_RNA')

smartseq2[["pct_counts_mito"]] <- PercentageFeatureSet(smartseq2, pattern = "^mt-")
smartseq2[["pct_counts_ribo"]] <- PercentageFeatureSet(smartseq2, pattern = "^Rp[sl][[:digit:]]")

VlnPlot(smartseq2, features = c("pct_counts_mito", "pct_counts_ribo"), pt.size = 1, ncol = 2)
RidgePlot(smartseq2, features = c("nFeature_RNA", "nCount_RNA", "pct_counts_mito", "pct_counts_ribo"), ncol = 2)

FeatureScatter(smartseq2, feature1 = "nFeature_RNA", feature2 = "pct_counts_mito", group.by = "orig.ident")



# Normalisation
smartseq2 <- NormalizeData(smartseq2, normalization.method = "LogNormalize", scale.factor = 1e6)


# Finding highly variable features
smartseq2 <- FindVariableFeatures(object = smartseq2, selection.method = "vst", nfeatures = 2500,
                            loess.span = 0.3, clip.max = "auto", num.bin = 20,
                            binning.method = "equal_width", verbose = TRUE,
                            mean.cutoff = c(0.1, Inf), dispersion.cutoff = c(1, Inf))

# Print the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = smartseq2), 10)

# Get highly variable gene info
head(HVFInfo(object = smartseq2))

hvg <- HVFInfo(object = smartseq2)

# Sort variable genes by standarised variance
hvg_by_standarisedvariance <- hvg[order(hvg$variance.standardized, decreasing = TRUE),]

# Scaling
all.genes <- rownames(smartseq2)
smartseq2 <- ScaleData(smartseq2, features = all.genes, assay = "RNA",
                       vars.to.regress = NULL, model.use = "linear", use.umi = FALSE,
                       do.scale = TRUE, do.center = TRUE, scale.max = 10,
                       block.size = 1000, min.cells.to.block = 3000, verbose = TRUE) # results are stored in ecs[["RNA"]]@scale.data)


# PCA
# Linear dimensionality reduction
smartseq2 <- RunPCA(smartseq2, features = VariableFeatures(object = smartseq2),
                    npcs = 50, rev.pca = FALSE, weight.by.var = TRUE, verbose = TRUE,
                    ndims.print = 1:5, nfeatures.print = 10, reduction.name = "pca",
                    reduction.key = "PC_", seed.use = 42)

print(smartseq2[["pca"]], dims = 1:5, nfeatures = 5)


# Visualising PCs
VizDimLoadings(smartseq2, dims = 1:2, reduction = "pca")

DimPlot(smartseq2, reduction = "pca")

# Heatmap of PCs
DimHeatmap(smartseq2, dims = 1:6, balanced = TRUE)

# Determining the dimensiojnality of the dataset
ElbowPlot(smartseq2)

DimHeatmap(smartseq2, dims = 9:17, balanced = TRUE)
DimHeatmap(smartseq2, dims = 1:9, balanced = TRUE)


# Clustering
smartseq2 <- FindNeighbors(smartseq2, dims = 1:13, k.param = 10)

smartseq2 <- FindClusters(object = smartseq2, graph.name = NULL,
                    modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                    node.sizes = NULL, resolution = seq(0.2, 0.8, 0.1), algorithm = 1, n.start = 1000,
                    n.iter = 10, random.seed = 0, temp.file.location = NULL,
                    edge.file.name = NULL, verbose = TRUE)


# How many cells in each cluster?
table(x = Idents(smartseq2))

# What is the fraction of cells in each cluster?
prop.table(x = table(x = Idents(smartseq2)))

# Clustering tree
library(clustree)

names(smartseq2[[]])
clustree(smartseq2, prefix = "RNA_snn_res.")

# Set cell identities
Idents(object = smartseq2) <- smartseq2[["RNA_snn_res.0.7"]]
table(x = Idents(smartseq2))

# Build phylogenetic tree
smartseq2 <- BuildClusterTree(smartseq2)

Tool(object = smartseq2, slot = 'BuildClusterTree')

PlotClusterTree(smartseq2)


## Data visualisation

# UMAP
smartseq2 <- RunUMAP(object = smartseq2, dims = 1:13, assay = NULL, n.neighbors = 10, #30
                     n.components = 2L, metric = "correlation", n.epochs = NULL,
                     learning.rate = 1, min.dist = 0.1, spread = 1,
                     set.op.mix.ratio = 1, local.connectivity = 1L,
                     repulsion.strength = 1, negative.sample.rate = 5, a = NULL,
                     b = NULL, seed.use = 123, metric.kwds = NULL,
                     angular.rp.forest = FALSE, reduction.key = "UMAP_", verbose = TRUE)

DimPlot(object = smartseq2, reduction = "umap", pt.size = 1)
DimPlot(object = smartseq2, reduction = "umap", pt.size = 1, split.by = "mouse", ncol = 2)

FeaturePlot(smartseq2, features = c("Sntn"), reduction = "umap", cols = c("lightgrey", "red"), order = TRUE, pt.size = 1)


# Find markers
markers_smartseq2 <- FindAllMarkers(smartseq2, only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.25)
markers_smartseq2 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)


### Removing non-ependymal cells (likely doublets)
table(Idents(smartseq2))
smartseq2 <- subset(x= smartseq2, idents = c(0,1,2,3,4,5,6))

# Find highly variable features
smartseq2 <- FindVariableFeatures(smartseq2, selection.method = "vst", nfeatures = 2500)

# Scaling
all.genes <- rownames(smartseq2)
smartseq2 <- ScaleData(smartseq2, features = all.genes, assay = "RNA",
                       vars.to.regress = NULL, model.use = "linear", use.umi = FALSE,
                       do.scale = TRUE, do.center = TRUE, scale.max = 10,
                       block.size = 1000, min.cells.to.block = 3000, verbose = TRUE) # results are stored in ecs[["RNA"]]@scale.data)

# PCA
# Linear dimensionality reduction
smartseq2 <- RunPCA(smartseq2, features = VariableFeatures(object = smartseq2),
                    npcs = 50, rev.pca = FALSE, weight.by.var = TRUE, verbose = TRUE,
                    ndims.print = 1:5, nfeatures.print = 10, reduction.name = "pca",
                    reduction.key = "PC_", seed.use = 42)

print(smartseq2[["pca"]], dims = 1:5, nfeatures = 5)


# Visualising PCs
VizDimLoadings(smartseq2, dims = 1:2, reduction = "pca")

DimPlot(smartseq2, reduction = "pca")

# Heatmap of PCs
DimHeatmap(smartseq2, dims = 1:6, balanced = TRUE)

# Determining the dimensiojnality of the dataset
ElbowPlot(smartseq2)

DimHeatmap(smartseq2, dims = 9:17, balanced = TRUE)
DimHeatmap(smartseq2, dims = 1:9, balanced = TRUE)


# Clustering
smartseq2 <- FindNeighbors(smartseq2, dims = 1:13, k.param = 10)

smartseq2 <- FindClusters(object = smartseq2, graph.name = NULL,
                          modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                          node.sizes = NULL, resolution = seq(0.2, 0.8, 0.1), algorithm = 1, n.start = 1000,
                          n.iter = 10, random.seed = 0, temp.file.location = NULL,
                          edge.file.name = NULL, verbose = TRUE)


# What is the fraction of cells in each cluster?
prop.table(x = table(x = Idents(smartseq2)))

# Clustering tree
clustree(smartseq2, prefix = "RNA_snn_res.")

# Set cell identities
Idents(object = smartseq2) <- smartseq2[["RNA_snn_res.0.5"]]
table(x = Idents(smartseq2))


## Data visualisation

# UMAP
smartseq2 <- RunUMAP(object = smartseq2, dims = 1:13, assay = NULL, n.neighbors = 10, #30
                     n.components = 2L, metric = "correlation", n.epochs = NULL,
                     learning.rate = 1, min.dist = 0.1, spread = 1,
                     set.op.mix.ratio = 1, local.connectivity = 1L,
                     repulsion.strength = 1, negative.sample.rate = 5, a = NULL,
                     b = NULL, seed.use = 123, metric.kwds = NULL,
                     angular.rp.forest = FALSE, reduction.key = "UMAP_", verbose = TRUE)

DimPlot(object = smartseq2, reduction = "umap", pt.size = 1)
DimPlot(object = smartseq2, reduction = "umap", pt.size = 1, split.by = "mouse", ncol = 2)

FeaturePlot(smartseq2, features = c("Sntn"), reduction = "umap", cols = c("lightgrey", "red"), order = TRUE, pt.size = 1)


# Find markers
markers_smartseq2 <- FindAllMarkers(smartseq2, only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.25)

# Fetch QC data from the object's metadata
names(smartseq2[[]])
smartseq2$seurat_clusters <- smartseq2$RNA_snn_res.0.5

smartseq2_qc_data <- FetchData(object = smartseq2, vars = c("seurat_clusters", "nCount_RNA", "nFeature_RNA", "pct_counts_mito", "pct_counts_ribo"))
colnames(smartseq2)

# Boxplot of QC metrics by cluster ID
smartseq2_reads <- ggplot(data = smartseq2_qc_data, aes(x = seurat_clusters, y = nCount_RNA)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_linedraw()

smartseq2_features <- ggplot(data = smartseq2_qc_data, aes(x = seurat_clusters, y = nFeature_RNA)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_linedraw()

smartseq2_mito <- ggplot(data = smartseq2_qc_data, aes(x = seurat_clusters, y = pct_counts_mito)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_linedraw()

smartseq2_reads + smartseq2_features + smartseq2_mito


### Removing cluster 6 and cluster 7 (res 0.5)
table(Idents(smartseq2))
smartseq2 <- subset(x= smartseq2, idents = c(0,1,2,3,4,5))

# Find highly variable features
smartseq2 <- FindVariableFeatures(smartseq2, selection.method = "vst", nfeatures = 2500)

# Scaling
all.genes <- rownames(smartseq2)
smartseq2 <- ScaleData(smartseq2, features = all.genes, assay = "RNA",
                       vars.to.regress = NULL, model.use = "linear", use.umi = FALSE,
                       do.scale = TRUE, do.center = TRUE, scale.max = 10,
                       block.size = 1000, min.cells.to.block = 3000, verbose = TRUE) # results are stored in ecs[["RNA"]]@scale.data)


# PCA
# Linear dimensionality reduction
smartseq2 <- RunPCA(smartseq2, features = VariableFeatures(object = smartseq2),
                    npcs = 50, rev.pca = FALSE, weight.by.var = TRUE, verbose = TRUE,
                    ndims.print = 1:5, nfeatures.print = 10, reduction.name = "pca",
                    reduction.key = "PC_", seed.use = 42)

print(smartseq2[["pca"]], dims = 1:5, nfeatures = 5)


# Visualising PCs
VizDimLoadings(smartseq2, dims = 1:2, reduction = "pca")

DimPlot(smartseq2, reduction = "pca")

# Heatmap of PCs
DimHeatmap(smartseq2, dims = 1:6, balanced = TRUE)

# Determining the dimensionality of the dataset
ElbowPlot(smartseq2)

DimHeatmap(smartseq2, dims = 9:17, balanced = TRUE)
DimHeatmap(smartseq2, dims = 1:9, balanced = TRUE)


# Clustering
smartseq2 <- FindNeighbors(smartseq2, dims = 1:11, k.param = 10)

smartseq2 <- FindClusters(object = smartseq2, graph.name = NULL,
                          modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                          node.sizes = NULL, resolution = seq(0.2, 0.8, 0.1), algorithm = 1, n.start = 1000,
                          n.iter = 10, random.seed = 0, temp.file.location = NULL,
                          edge.file.name = NULL, verbose = TRUE)


# How many cells in each cluster?
table(x = Idents(smartseq2))

# What is the fraction of cells in each cluster?
prop.table(x = table(x = Idents(smartseq2)))

# Clustering tree
clustree(smartseq2, prefix = "RNA_snn_res.")

# Set cell identities
Idents(object = smartseq2) <- smartseq2[["RNA_snn_res.0.6"]]
table(x = Idents(smartseq2))

# Build phylogenetic tree
smartseq2 <- BuildClusterTree(smartseq2)

Tool(object = smartseq2, slot = 'BuildClusterTree')

PlotClusterTree(smartseq2)

## Data visualisation

# UMAP
smartseq2 <- RunUMAP(object = smartseq2, dims = 1:11, assay = NULL, n.neighbors = 10, #30
                     n.components = 2L, metric = "correlation", n.epochs = NULL,
                     learning.rate = 1, min.dist = 0.1, spread = 1,
                     set.op.mix.ratio = 1, local.connectivity = 1L,
                     repulsion.strength = 1, negative.sample.rate = 5, a = NULL,
                     b = NULL, seed.use = 1, metric.kwds = NULL,
                     angular.rp.forest = FALSE, reduction.key = "UMAP_", verbose = TRUE)

DimPlot(object = smartseq2, reduction = "umap", pt.size = 1)
DimPlot(object = smartseq2, reduction = "umap", pt.size = 1, split.by = "mouse", ncol = 2)

FeaturePlot(smartseq2, features = c("Jun"), reduction = "umap", cols = c("lightgrey", "red"), order = TRUE, pt.size = 1)


# Find markers
markers_smartseq2 <- FindAllMarkers(smartseq2, only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.25)

res06_c1_vs_c2 <- FindMarkers(smartseq2, ident.1 = 1, ident.2 = 2, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res06_c2_vs_c1 <- FindMarkers(smartseq2, ident.1 = 2, ident.2 = 1, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res06_c3_vs_c4 <- FindMarkers(smartseq2, ident.1 = 3, ident.2 = 4, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res06_c4_vs_c3 <- FindMarkers(smartseq2, ident.1 = 4, ident.2 = 3, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)


# Merge clusters 1 and 2 (lateral cells, immature) and 3 and 4 (lateral cells, activated)
Idents(smartseq2) <- smartseq2$RNA_snn_res.0.6
smartseq2 <- RenameIdents(smartseq2, `0` = "0", `1` = "1", `2` = "1",
                          `3` = "2", `4` = "2", `5` = "3",
                          `6` = "4")

DimPlot(smartseq2, reduction = "umap", label = FALSE, pt.size = 1)

smartseq2$merged_clusters <- Idents(smartseq2)

table(smartseq2$merged_clusters)

# Assign cell type
Idents(smartseq2) <- smartseq2$merged_clusters
smartseq2 <- RenameIdents(smartseq2, `0` = "Lateral (mature)", `1` = "Lateral (immature)",
                          `2` = "Lateral (activated)",
                          `3` = "Ventral",
                          `4` = "Proliferating")

DimPlot(smartseq2, reduction = "umap", label = FALSE, pt.size = 1)

smartseq2$cell_subtype <- Idents(smartseq2)

Idents(smartseq2) <- smartseq2$cell_subtype

# Assign colours to cell types
DimPlot(object = smartseq2, reduction = "umap", label = FALSE, pt.size = 1,
        cols = c( "Lateral (immature)" = "#fbd347", "Lateral (activated)" = "#f1443d",
                  "Lateral (mature)" = "#fdae46",
                  "Ventral" = "#43a5f2",
                  "Proliferating" = "#32cd32")) + NoLegend()


table(Idents(smartseq2))

subtype_markers <- FindAllMarkers(object = smartseq2, genes.use = NULL, only.pos = TRUE,
                                     logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)

table(Idents(smartseq2))

smartseq2$original_cluster <- smartseq2$cell_subtype

table(smartseq2$original_cluster) # This stores the cluster assigned by Seurat clustering before integrating datasets




saveRDS(smartseq2, file = "~/path/to/smartseq2_ependymal_cells.rds")


# QC plots
# Fetch QC data from the object's metadata
names(smartseq2[[]])
Idents(smartseq2) <- smartseq2$cell_subtype

smartseq2_qc_data <- FetchData(object = smartseq2, vars = c("cell_subtype", "mouse", "nCount_RNA", "nFeature_RNA"))
colnames(smartseq2_qc_data)

smartseq2_qc_data <- as.data.frame(smartseq2_qc_data)
head(smartseq2_qc_data)

# Boxplot of reads and features by ependymal cell subtype
smartseq2_reads <- ggplot(smartseq2_qc_data,
                          aes(factor(cell_subtype, levels = c("Lateral (immature)", "Proliferating ", "Lateral (mature)",
                                                              "Lateral (activated)", "Ventral")),
                                     y = nCount_RNA)) +
  geom_boxplot(outlier.size = 0.4, fill = "grey90") +
  #theme_linedraw() +
  theme_classic() +
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  theme(axis.text.x = element_blank())

smartseq2_features <- ggplot(smartseq2_qc_data,
                             aes(factor(cell_subtype, levels = c("Lateral (immature)", "Proliferating ", "Lateral (mature)",
                                                                 "Lateral (activated)", "Ventral")),
                                 y = nFeature_RNA)) +
  geom_boxplot(outlier.size = 0.4, fill = "grey90") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x='')

smartseq2_reads + smartseq2_features

# Boxplot of QC metrics by sample/mouse
reads_by_mouse <- ggplot(smartseq2_qc_data, aes(x = mouse, y = nCount_RNA)) +
  geom_boxplot(outlier.size = 0.4, fill = "grey90") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(x='')

features_by_mouse <- ggplot(smartseq2_qc_data, aes(x = mouse, y = nFeature_RNA)) +
  geom_boxplot(outlier.size = 0.4, fill = "grey90") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x='')

reads_by_mouse + features_by_mouse


# Plot together
plot_grid(reads_by_mouse, smartseq2_reads,
          features_by_mouse, smartseq2_features, ncol = 2, align = 'hv', axis = 'lr')













### Integrating Smart-seq2 and 10x ependymal cell datasets
ecs_smartseq2 <- readRDS(file = "~/Documents/Single-cell_RNA-seq/Smart-Seq2/STARsolo/smartseq2_ependymal_cells.rds")
ecs_10x <- readRDS(file = "~/Documents/Single-cell_RNA-seq/10X/data analysis/objects/young_ecs_subtypes_without_astrolike_cells.rds")

names(ecs_10x[[]])
table(ecs_10x$orig.ident) # add "mouse" and "tech" columns

mouse <- ecs_10x$orig.ident
ecs_10x <- AddMetaData(ecs_10x, metadata = mouse, col.name = "mouse")

ecs_10x # 1901 cells
tech <- rep("10X", 1901)
ecs_10x <- AddMetaData(ecs_10x, metadata = tech, col.name = "tech")

names(ecs_10x[[]])
names(ecs_smartseq2[[]])

# Idents are cell subtype
table(ecs_10x$cell_subtype)
table(Idents(ecs_10x))
table(ecs_smartseq2$cell_subtype)
table(Idents(ecs_smartseq2))

ecs_10x <- AddMetaData(ecs_10x, metadata = mouse, col.name = "mouse")

table(ecs_10x[["orig.ident"]])
table(ecs_10x[["mouse"]])

# The metadata column cell_subtype contains the original cluster info


ecs_smartseq2
DimPlot(ecs_smartseq2, pt.size = 1) + NoLegend()
DimPlot(ecs_10x, pt.size = 1)

# Print metadata
names(ecs_smartseq2[[]])
names(ecs_10x[[]])


# Merge the Smart-seq2 and 10x ependymal cell datasets
all_ecs_integrated <- merge(x = ecs_10x, y = ecs_smartseq2, merge.data = TRUE)

all_ecs_integrated # OK!

names(all_ecs_integrated[[]])
table(all_ecs_integrated[["tech"]])
table(all_ecs_integrated[["age"]])


# Split objects
all_ecs_integrated_list <- SplitObject(all_ecs_integrated, split.by = "tech")

all_ecs_integrated_list

# Normalisation
for (i in 1:length(all_ecs_integrated_list)) {
  all_ecs_integrated_list[[i]] <- NormalizeData(all_ecs_integrated_list[[i]], verbose = FALSE, normalization.method = "LogNormalize",
                                                scale.factor = 10000)
  all_ecs_integrated_list[[i]] <- FindVariableFeatures(all_ecs_integrated_list[[i]], selection.method = "vst", 
                                                       nfeatures = 2500, verbose = FALSE)
}

# Integrating
# Find anchors
all_ecs_integrated_anchors <- FindIntegrationAnchors(object.list = all_ecs_integrated_list, dims = 1:30)

all_ecs_integrated <- IntegrateData(anchorset = all_ecs_integrated_anchors, dims = 1:30)

all_ecs_integrated

DefaultAssay(all_ecs_integrated)

# Running the standard workflow for visualization and clustering
DefaultAssay(all_ecs_integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
all_ecs_integrated <- ScaleData(all_ecs_integrated, verbose = FALSE)
all_ecs_integrated <- RunPCA(all_ecs_integrated, npcs = 30, verbose = FALSE)

# Plot PCs
DimPlot(object = all_ecs_integrated, dims = c(1,2), reduction = "pca", group.by = "tech")
print(x = all_ecs_integrated[["pca"]], dims = 1, nfeatures = 10)


# What PCs explain most of the variance?
ElbowPlot(object = all_ecs_integrated)

DimHeatmap(object = all_ecs_integrated, dims = 9:17, balanced = TRUE, fast = FALSE, cells = 500)
DimHeatmap(object = all_ecs_integrated, dims = 1:6, balanced = TRUE, fast = FALSE, cells = 500)

# UMAP
all_ecs_integrated <- RunUMAP(object = all_ecs_integrated, dims = 1:15, seed.use = 12, learning.rate = 99,
                              n.neighbors = 30, spread = 1, min.dist = 0.1)

DimPlot(object = all_ecs_integrated, reduction = "umap", label = F, repel = TRUE, pt.size = 0.8) + NoLegend()

DimPlot(object = all_ecs_integrated, reduction = "umap", split.by = "tech", label = F, repel = TRUE, pt.size = 1) + NoLegend()
DimPlot(object = all_ecs_integrated, reduction = "umap", group.by =  "tech", label = FALSE, repel = TRUE, pt.size = 1)



# Clustering
DefaultAssay(all_ecs_integrated) <- "integrated"

names(all_ecs_integrated[[]])


# Remove earlier clustering info
columns.to.remove <- all_ecs_integrated@meta.data %>% 
  names() %>% 
  .[stringr::str_detect(., 'RNA_snn')]

for(i in columns.to.remove) {
  all_ecs_integrated[[i]] <- NULL
}

all_ecs_integrated$original_cluster <- NULL # the original cell subtype assignment is in cell_subtype

names(all_ecs_integrated[[]])


# Finding neighbours
all_ecs_integrated <- FindNeighbors(object = all_ecs_integrated, dims = 1:15, k.param = 10)

all_ecs_integrated <- FindClusters(object = all_ecs_integrated, graph.name = NULL,
                               modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                               node.sizes = NULL, resolution = seq(0.2, 0.8, 0.1), algorithm = 1, n.start = 1000,
                               n.iter = 10, random.seed = 0, temp.file.location = NULL,
                               edge.file.name = NULL, verbose = TRUE)

# Clustering trees
library(clustree)

names(all_ecs_integrated[[]])

clustree(all_ecs_integrated, prefix = "integrated_snn_res.")


# Set cell identities
Idents(object = all_ecs_integrated) <- all_ecs_integrated[["integrated_snn_res.0.7"]]
table(x = Idents(all_ecs_integrated))

DimPlot(object = all_ecs_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

# Plot UMAP
DimPlot(object = all_ecs_integrated, reduction = "umap", split.by = "tech", label = TRUE, repel = TRUE, pt.size = 1, cols = my_cols)
DimPlot(object = all_ecs_integrated, reduction = "umap", group.by = "tech", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(object = all_ecs_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(object = all_ecs_integrated, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 1, split.by = "mouse", ncol = 4)


# Normalise RNA counts for visualisation purposes
DefaultAssay(all_ecs_integrated) <- "RNA" # use RNA assay for visualisation and differential expression analysis

ecs_integrated <- NormalizeData(all_ecs_integrated, verbose = FALSE, scale.factor = 10000, normalization.method = "LogNormalize")
FeaturePlot(all_ecs_integrated, c("Sntn"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 1)


# Differential expression analysis
DefaultAssay(all_ecs_integrated) <- "RNA" # use RNA assay for visualisation and DE analysis

# Find markers
markers_ecs_res07 <- FindAllMarkers(all_ecs_integrated, only.pos = TRUE, test.use = "wilcox",
                                    logfc.threshold = 0.25, min.pct = 0.2, pseudocount.use = 1)

res07_c9_vs_c4 <- FindMarkers(all_ecs_integrated, ident.1 = 9, ident.2 = 4, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)


DimPlot(object = all_ecs_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.8)


# Removing microglial cells (cluster 10)
table(Idents(all_ecs_integrated))
all_ecs_integrated <- subset(x= all_ecs_integrated, idents = c(0,1,2,3,4,5,6,7,8,9,11))

DefaultAssay(all_ecs_integrated) <- "RNA"

# Split objects
all_ecs_integrated_list <- SplitObject(all_ecs_integrated, split.by = "tech")

all_ecs_integrated_list

# Normalisation
for (i in 1:length(all_ecs_integrated_list)) {
  all_ecs_integrated_list[[i]] <- NormalizeData(all_ecs_integrated_list[[i]], verbose = FALSE, normalization.method = "LogNormalize",
                                                scale.factor = 10000)
  all_ecs_integrated_list[[i]] <- FindVariableFeatures(all_ecs_integrated_list[[i]], selection.method = "vst", 
                                                       nfeatures = 2500, verbose = FALSE)
}

# Integrating
# Find anchors
all_ecs_integrated_anchors <- FindIntegrationAnchors(object.list = all_ecs_integrated_list, dims = 1:30)

all_ecs_integrated <- IntegrateData(anchorset = all_ecs_integrated_anchors, dims = 1:30)

all_ecs_integrated

DefaultAssay(all_ecs_integrated)

# Running the standard workflow for visualization and clustering
DefaultAssay(all_ecs_integrated) <- "integrated"

all_ecs_integrated <- ScaleData(all_ecs_integrated, verbose = FALSE)
all_ecs_integrated <- RunPCA(all_ecs_integrated, npcs = 30, verbose = FALSE)

# Plot PCs
DimPlot(object = all_ecs_integrated, dims = c(1,2), reduction = "pca", group.by = "tech")
print(x = all_ecs_integrated[["pca"]], dims = 1, nfeatures = 10)

# What PCs explain most of the variance?
ElbowPlot(object = all_ecs_integrated)

DimHeatmap(object = all_ecs_integrated, dims = 10:18, balanced = TRUE, fast = FALSE, cells = 500)
DimHeatmap(object = all_ecs_integrated, dims = 1:6, balanced = TRUE, fast = FALSE, cells = 500)

# UMAP
all_ecs_integrated <- RunUMAP(object = all_ecs_integrated, dims = 1:14, seed.use = 456, learning.rate = 90,
                              n.neighbors = 10, spread = 1, min.dist = 0.3)

DimPlot(object = all_ecs_integrated, reduction = "umap", label = F, repel = TRUE, pt.size = 0.8) + NoLegend()
DimPlot(object = all_ecs_integrated, reduction = "umap", split.by = "tech", label = F, repel = TRUE, pt.size = 1) + NoLegend()


# Clustering
DefaultAssay(all_ecs_integrated) <- "integrated"

names(all_ecs_integrated[[]])


# Remove earlier clustering info
columns.to.remove <- all_ecs_integrated@meta.data %>% 
  names() %>% 
  .[stringr::str_detect(., 'RNA_snn')]

for(i in columns.to.remove) {
  all_ecs_integrated[[i]] <- NULL
}

names(all_ecs_integrated[[]])


# Finding neighbours
all_ecs_integrated <- FindNeighbors(object = all_ecs_integrated, dims = 1:14, k.param = 10)

all_ecs_integrated <- FindClusters(object = all_ecs_integrated, graph.name = NULL,
                                   modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                                   node.sizes = NULL, resolution = seq(0.1, 0.8, 0.1), algorithm = 1, n.start = 1000,
                                   n.iter = 10, random.seed = 0, temp.file.location = NULL,
                                   edge.file.name = NULL, verbose = TRUE)

# Clustering trees
clustree(all_ecs_integrated, prefix = "integrated_snn_res.")


# Set cell identities
Idents(object = all_ecs_integrated) <- all_ecs_integrated[["integrated_snn_res.0.5"]]
table(x = Idents(all_ecs_integrated))

DimPlot(object = all_ecs_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

# Build phylogenetic tree with Seurat
DefaultAssay(all_ecs_integrated) <- "integrated"
all_ecs_integrated <- BuildClusterTree(all_ecs_integrated)

Tool(object = all_ecs_integrated, slot = 'BuildClusterTree')

PlotClusterTree(all_ecs_integrated)

# Remove node labels
# Pull the tree
data_tree <- Tool(object = all_ecs_integrated, slot = "BuildClusterTree")
# Plot the tree
ape::plot.phylo(x = data_tree, direction = "downwards")


# Plot UMAP
DimPlot(object = all_ecs_integrated, reduction = "umap", split.by = "tech", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(object = all_ecs_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(object = all_ecs_integrated, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 1, split.by = "mouse", ncol = 4)


# Normalise RNA counts for visualisation purposes
DefaultAssay(all_ecs_integrated) <- "RNA" # use RNA assay for visualisation and DE analysis

all_ecs_integrated <- NormalizeData(all_ecs_integrated, verbose = FALSE, scale.factor = 10000, normalization.method = "LogNormalize")
FeaturePlot(all_ecs_integrated, c("Pax6"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 0.8, split.by = "ident")

# Differential expression analysis
# Find markers
markers_ecs_res05 <- FindAllMarkers(all_ecs_integrated, only.pos = TRUE, test.use = "wilcox",
                                    logfc.threshold = 0.25, min.pct = 0.15, pseudocount.use = 1)

DimPlot(object = all_ecs_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

# Pair-wise differential expression analysis
res05_c4_vs_c0 <- FindMarkers(all_ecs_integrated, ident.1 = 4, ident.2 = 0, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res05_c0_vs_c4 <- FindMarkers(all_ecs_integrated, ident.1 = 0, ident.2 = 4, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1) # merge clusters 0 and 4
res05_c1_vs_c3 <- FindMarkers(all_ecs_integrated, ident.1 = 1, ident.2 = 3, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res05_c3_vs_c1 <- FindMarkers(all_ecs_integrated, ident.1 = 3, ident.2 = 1, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1) # merge clusters 1 and 3
LATimm_vs_LATmat <- FindMarkers(all_ecs_integrated, ident.1 = "Lateral (immature)", ident.2 = "Lateral (mature)", only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)

LATmat_vs_LATimm <- FindMarkers(all_ecs_integrated, ident.1 = "Lateral (mature)", ident.2 = "Lateral (immature)", only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)


# Merging and making sense of clusters
Idents(object = all_ecs_integrated) <- all_ecs_integrated[["integrated_snn_res.0.5"]]

# Assign cell subtype
all_ecs_integrated <- RenameIdents(all_ecs_integrated, `0` = "Lateral (immature)", `1` = "Lateral (activated)",
                                   `2` = "Lateral (mature)",
                                   `3` = "Lateral (activated)", `4` = "Lateral (immature)",
                                   `5` = "Ventral",
                                   `6` = "Dorsal",
                                   `7` = "Lateral (proliferating)")

DimPlot(all_ecs_integrated, reduction = "umap", label = FALSE, pt.size = 0.8) + NoLegend()
DimPlot(all_ecs_integrated, reduction = "umap", label = FALSE, pt.size = 1, group.by = "tech")
DimPlot(all_ecs_integrated, reduction = "umap", label = FALSE, pt.size = 1, split.by = "mouse", ncol = 4) # Figure S2C

# Figure 2A
DimPlot(object = all_ecs_integrated, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 1,
        cols = c( "Lateral (immature)" = "#fbd347", "Lateral (activated)" = "#f1443d",
                  "Lateral (mature)" = "#fdae46", "Lateral (proliferating)" = "#24ce30",
                  "Ventral" = "#43a5f2", "Dorsal" = "#a647b8")) + NoAxes()

table(Idents(all_ecs_integrated))

all_ecs_integrated$cell_subtype_integrated <- Idents(all_ecs_integrated)

all_ecs_integrated$cell_subtype_integrated <- factor(x = all_ecs_integrated$cell_subtype_integrated,
                                                     levels = c("Lateral (immature)", "Lateral (proliferating)",
                                                                "Lateral (mature)",
                                                                "Lateral (activated)",
                                                                "Dorsal", "Ventral"))


Idents(all_ecs_integrated) <- all_ecs_integrated$cell_subtype_integrated
ecs_palette <- c("#fbd347", "#24ce30", "#fdae46", "#f1443d", "#a647b8", "#43a5f2")

DimPlot(object = all_ecs_integrated, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 0.8,
        cols = ecs_palette) # Figure 2A

table(all_ecs_integrated$cell_subtype_integrated)



# Build phylogenetic tree with Seurat
DefaultAssay(all_ecs_integrated) <- "integrated"
all_ecs_integrated <- BuildClusterTree(all_ecs_integrated)

Tool(object = all_ecs_integrated, slot = 'BuildClusterTree')

PlotClusterTree(all_ecs_integrated)

# Pull the tree to remove node labels
data_tree <- Tool(object = all_ecs_integrated, slot = "BuildClusterTree")

# Plot the tree
ape::plot.phylo(x = data_tree, direction = "downwards")


Idents(all_ecs_integrated) <- all_ecs_integrated$cell_subtype_integrated

# Saving the object
saveRDS(all_ecs_integrated, file = "~/path/to/all_ecs_integrated_subtypes.rds")

DimPlot(object = all_ecs_integrated, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 0.8,
        cols = ecs_palette)


# Differential expression analysis
DefaultAssay(all_ecs_integrated) <- "RNA" # use RNA assay for visualisation and DE analysis

# Find markers
markers_ecs_subtypes <- FindAllMarkers(all_ecs_integrated, only.pos = TRUE, test.use = "wilcox",
                                    logfc.threshold = 0.25, min.pct = 0.15, pseudocount.use = 1)

table(Idents(all_ecs_integrated))



# Inferring ependymal functions: GO terms
# Figure 2E
library(cowplot)

# Dorsal cells
dorsal_goterms_selected <- read.csv(file = "~/path/to/selected_GO_terms_DOR-FDR.csv",
                                       header = TRUE)

names(dorsal_goterms_selected)

dorsal_goterms_selected <- arrange(dorsal_goterms_selected, adjusted_p_value) # sort by decreasing average log fold change

# Bar plot
dor_go <- ggplot(data = dorsal_goterms_selected, aes(x = reorder(stringr::str_to_sentence(term_name), negative_log10_of_adjusted_p_value),
                                              y = negative_log10_of_adjusted_p_value)) +
  geom_bar(stat = "identity", fill = "#56B4E9") +
  coord_flip() +
  theme_cowplot(12) +
  xlab(NULL)


# Ventral cells
ventral_goterms_selected <- read.csv(file = "~/path/to/selected_GO_terms_VEN-FDR.csv",
                                    header = TRUE)

ventral_goterms_selected <- arrange(ventral_goterms_selected, adjusted_p_value) # sort by decreasing average log fold change


# Bar plot
ven_go <- ggplot(data = ventral_goterms_selected, aes(x = reorder(stringr::str_to_sentence(term_name), negative_log10_of_adjusted_p_value),
                                           y = negative_log10_of_adjusted_p_value)) +
  geom_bar(stat = "identity", fill = "#56B4E9") +
  coord_flip() +
  theme_cowplot(12) +
  xlab(NULL)

# Lateral cells (mature)
lateral_mature_goterms_selected <- read.csv(file = "~/path/to/selected_GO_terms_LATmat_FDR.csv",
                                     header = TRUE)

lateral_mature_goterms_selected <- arrange(lateral_mature_goterms_selected, adjusted_p_value) # sort by decreasing average log fold change

names(lateral_mature_goterms_selected)

# Bar plot
latact_go <- ggplot(data = lateral_mature_goterms_selected, aes(x = reorder(stringr::str_to_sentence(term_name), negative_log10_of_adjusted_p_value),
                                            y = negative_log10_of_adjusted_p_value)) +
  geom_bar(stat = "identity", fill = "#56B4E9") +
  coord_flip() +
  theme_cowplot(12) +
  xlab(NULL)


# Lateral cells (activated)
lateral_activated_goterms_selected <- read.csv(file = "~/path/to/selected_GO_terms_LATact_FDR.csv",
                                            header = TRUE)

lateral_activated_goterms_selected <- arrange(lateral_activated_goterms_selected, adjusted_p_value) # sort by decreasing average log fold change


# Bar plot
latmat_go <- ggplot(data = lateral_activated_goterms_selected, aes(x = reorder(stringr::str_to_sentence(term_name), negative_log10_of_adjusted_p_value),
                                                   y = negative_log10_of_adjusted_p_value)) +
  geom_bar(stat = "identity", fill = "#56B4E9") +
  coord_flip() +
  theme_cowplot(12) +
  xlab(NULL)


# Plot together
plot_grid(dor_go, ven_go,
          latact_go, latmat_go,
          ncol = 2, align = 'v', axis = 'lr') # Figure 2E




# Supplementary info

# Data QC plots for supplementary figures

# Fetch QC data from the object's metadata
names(all_ecs_integrated[[]])

all_ecs_integrated_qc_data <- FetchData(object = all_ecs_integrated, vars = c("cell_subtype_integrated", "mouse", "nCount_RNA", "nFeature_RNA", "tech"))
colnames(all_ecs_integrated_qc_data)


all_ecs_integrated_qc_data <- as.data.frame(all_ecs_integrated_qc_data)
head(all_ecs_integrated_qc_data)

# Boxplot of QC metrics by cluster ID
ecs_integrated_reads <- ggplot(data = all_ecs_integrated_qc_data, aes(x = cell_subtype_integrated, y = nCount_RNA, fill = tech)) +
  geom_boxplot(outlier.size = 0.4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid = element_blank()) +
  facet_wrap(~tech, scales = "free") +
  labs(x='')

ecs_integrated_features <- ggplot(data = all_ecs_integrated_qc_data, aes(x = cell_subtype_integrated, y = nFeature_RNA, fill = tech)) +
  geom_boxplot(outlier.size = 0.4) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x='')

ecs_integrated_reads + ecs_integrated_features # Figure S2A




# Boxplot of QC metrics by sample/mouse
ecs_integrated_reads_by_mouse <- ggplot(data = all_ecs_integrated_qc_data, aes(x = mouse, y = nCount_RNA, fill = tech)) +
  geom_boxplot(outlier.size = 0.4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid = element_blank()) +
  facet_wrap(~tech, scales = "free") +
  labs(x='')

ecs_integrated_features_by_mouse <- ggplot(data = all_ecs_integrated_qc_data, aes(x = mouse, y = nFeature_RNA, fill = tech)) +
  geom_boxplot(outlier.size = 0.4) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x='')

ecs_integrated_reads_by_mouse + ecs_integrated_features_by_mouse



# Plot together
library(cowplot)
plot_grid(ecs_integrated_reads_by_mouse,
          ecs_integrated_features_by_mouse,
          ecs_integrated_reads,
          ecs_integrated_features, ncol = 2, align = 'hv', axis = 'lr')



# Median number of UMI/read counts
median(all_ecs_integrated$nCount_RNA)

# Median number of UMI/read counts
median(all_ecs_integrated$nFeature_RNA)


# UMAP plot split by mouse showing contribution of all samples to all clusters
DimPlot(object = all_ecs_integrated, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 0.8,
        cols = ecs_palette, split.by = "mouse", ncol = 4) # Figure 1A




# Stacked bar plots showing per mouse contribution to clusters
# How much does each mouse contribute to each cluster?
table(Idents(all_ecs_integrated), all_ecs_integrated@meta.data$mouse)
mouse_to_clusters <- table(Idents(all_ecs_integrated), all_ecs_integrated@meta.data$mouse)

head(mouse_to_clusters)

# Stacked barplot showing cells/tissue for each cluster
mouse_to_clusters <- as.data.frame(mouse_to_clusters)
head(mouse_to_clusters)

mouse_to_clusters <- mouse_to_clusters %>% 
  rename(cluster = Var1, mouse = Var2, count = Freq)

head(mouse_to_clusters)

# Percent stacked barplot
ggplot(mouse_to_clusters, aes(fill = mouse, y = count, x = cluster)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) # Figure S2C


# Stacked bar plot showing percentage of cells in each cluster
ggplot(cluster_percent, aes(fill=cluster, y=count, x=1)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = ecs_palette) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(breaks = seq(0, 2792, 500)) +
  #use positions to plot labels
  geom_text(aes(label = paste0(count)), 
            position = position_stack(vjust = 0.5), size = 4) +
labs(x='') # Figure 2A

round(prop.table(table((Idents(all_ecs_integrated))))*100, digits=1)




# Plotting gene expression patterns on UMAP plots with ggplot2
DefaultAssay(all_ecs_integrated) <- "RNA"

single_cell_object <- all_ecs_integrated

extract_plot_data <- function(single_cell_object, gene = NULL) {
  
  temp <- FeaturePlot(single_cell_object, features = c(gene), reduction = "umap", order = TRUE) 
  
  plot_df <- tibble(
    x = temp$data$UMAP_1,
    y = temp$data$UMAP_2,
    value = c(temp$data[,gene])
  )
  
  return(plot_df)
}


plot_df <- tibble(
  gene_name = c("Sntn")
) %>% mutate(map(gene_name, ~extract_plot_data(all_ecs_integrated, .x))) %>%
  unnest() %>% print()

plot_title <- plot_df$gene_name %>% unique()

ggplot(plot_df %>% filter(value == 0), aes(x, y)) +
  geom_point(color = 'lightgrey', size = 0.6) +
  geom_point(data = plot_df %>% filter(value != 0), aes(x, y, color = value), size = 0.6) +
  labs(x = 'UMAP 1', y = 'UMAP 2', title = NULL, subtitle = NULL) +
  viridis::scale_color_viridis(option = "viridis", name = 'Scaled\nexpression') +
  facet_wrap(~gene_name) +
  # labs(title = gene_name) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 12))




