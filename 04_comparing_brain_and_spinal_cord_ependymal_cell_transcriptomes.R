

## Comparing the transcriptomes of ependymal cell from across the central nervous system

# Integrationg the 10x dataset of spinal cord ependymal cells from young mice with ependymal cell transcriptomes from Zeisel et al 2018 (DOI: 10.1016/j.cell.2018.06.021)

# Load packages
library(Seurat)
library(tidyverse)
options(future.globals.maxSize = 4000 * 1024^5)


## Loading all the count matrices from Zeisel by tissue

# Amygdala
data_dir <- "~/path/to/zeisel_integration/outs_by_tissue/amygdala/"
sample_IDs <- list.files("~/path/to/zeisel_integration/outs_by_tissue/amygdala/")

zeisel_amygdala_matrices <- sapply(sample_IDs, function(i){
  amygdala <- Read10X(file.path(data_dir,i,"outs/filtered_feature_bc_matrix/"))
  colnames(amygdala) <- paste(i, colnames(amygdala), sep = "_")
  amygdala
})

head(colnames(zeisel_amygdala_matrices$`10X06_2`))

amygdala_data <- do.call("cbind", zeisel_amygdala_matrices)

amygdala_data[1:5, 1:5]

# Create Seurat object
amygdala <- CreateSeuratObject(counts = amygdala_data, min.cells = 3, min.features = 500)
amygdala # 18859 cells

# Prepare metadata for the amygdala
dataset <- rep("Zeisel", 18859)
tissue <- rep("amygdala", 18859)

amygdala <- AddMetaData(amygdala, dataset, col.name = "dataset")
amygdala <- AddMetaData(amygdala, tissue, col.name = "tissue")

names(amygdala[[]])
table(amygdala$tissue)

# Subset cells expressing Foxj1
amygdala_ecs <- subset(amygdala, subset = Foxj1 > 0)
amygdala_ecs


# Hypothalamus
data_dir <- "~/path/to/zeisel_integration/outs_by_tissue/hypothalamus/"
sample_IDs <- list.files("~/path/to/zeisel_integration/outs_by_tissue/hypothalamus/")

zeisel_hypothalamus_matrices <- sapply(sample_IDs, function(i){
  hypothalamus <- Read10X(file.path(data_dir,i,"outs/filtered_feature_bc_matrix/"))
  colnames(hypothalamus) <- paste(i, colnames(hypothalamus), sep = "_")
  hypothalamus
})

hypothalamus_data <- do.call("cbind", zeisel_hypothalamus_matrices)

hypothalamus_data[1:5, 1:5]
dim(hypothalamus_data)

# Create Seurat object
hypothalamus <- CreateSeuratObject(counts = hypothalamus_data, min.cells = 3, min.features = 500)
hypothalamus # 27624 cells

# Prepare metadata for the hypothalamus
dataset <- rep("Zeisel", 27624)
tissue <- rep("hypothalamus", 27624)

hypothalamus <- AddMetaData(hypothalamus, dataset, col.name = "dataset")
hypothalamus <- AddMetaData(hypothalamus, tissue, col.name = "tissue")

names(hypothalamus[[]])
table(hypothalamus$tissue)

# Subset cells expressing Foxj1
hypothalamus_ecs <- subset(hypothalamus, subset = Foxj1 > 0)
hypothalamus_ecs


# Hippocampus
data_dir <- "~/path/to/outs_by_tissue/hippocampus/"
sample_IDs <- list.files("~/path/to/zeisel_integration/outs_by_tissue/hippocampus/")

zeisel_hippocampus_matrices <- sapply(sample_IDs, function(i){
  hippocampus <- Read10X(file.path(data_dir,i,"outs/filtered_feature_bc_matrix/"))
  colnames(hippocampus) <- paste(i, colnames(hippocampus), sep = "_")
  hippocampus
})

hippocampus_data <- do.call("cbind", zeisel_hippocampus_matrices)

hippocampus_data[1:5, 1:5]
dim(hippocampus_data)

# Create Seurat object
hippocampus <- CreateSeuratObject(counts = hippocampus_data, min.cells = 3, min.features = 500)
hippocampus # 17444 cells

# Prepare metadata for the hippocampus
dataset <- rep("Zeisel", 17444)
tissue <- rep("hippocampus", 17444)

hippocampus <- AddMetaData(hippocampus, dataset, col.name = "dataset")
hippocampus <- AddMetaData(hippocampus, tissue, col.name = "tissue")

table(hippocampus$tissue)

# Subset cells expressing Foxj1
hippocampus_ecs <- subset(hippocampus, subset = Foxj1 > 0)
hippocampus_ecs


# Medulla
data_dir <- "~/path/to/zeisel_integration/outs_by_tissue/medulla/"
sample_IDs <- list.files("~/path/to/zeisel_integration/outs_by_tissue/medulla/")

zeisel_medulla_matrices <- sapply(sample_IDs, function(i){
  medulla <- Read10X(file.path(data_dir,i,"outs/filtered_feature_bc_matrix/"))
  colnames(medulla) <- paste(i, colnames(medulla), sep = "_")
  medulla
  
})

medulla_data <- do.call("cbind", zeisel_medulla_matrices)

medulla_data[1:5, 1:5]
dim(medulla_data)

# Create Seurat object
medulla <- CreateSeuratObject(counts = medulla_data, min.cells = 3, min.features = 500)
medulla # 17444 cells

# Prepare metadata for the medulla
dataset <- rep("Zeisel", 68631)
tissue <- rep("medulla", 68631)

medulla <- AddMetaData(medulla, dataset, col.name = "dataset")
medulla <- AddMetaData(medulla, tissue, col.name = "tissue")

table(medulla$tissue)

# Subset cells expressing Foxj1
medulla_ecs <- subset(medulla, subset = Foxj1 > 0)
medulla_ecs

# Midbrain
data_dir <- "~/path/to/zeisel_integration/outs_by_tissue/midbrain/"
sample_IDs <- list.files("~/path/to/zeisel_integration/outs_by_tissue/midbrain/")

zeisel_midbrain_matrices <- sapply(sample_IDs, function(i){
  midbrain <- Read10X(file.path(data_dir,i,"outs/filtered_feature_bc_matrix/"))
  colnames(midbrain) <- paste(i, colnames(midbrain), sep = "_")
  midbrain
  
})

midbrain_data <- do.call("cbind", zeisel_midbrain_matrices)

dim(midbrain_data)

# Create Seurat object
midbrain <- CreateSeuratObject(counts = midbrain_data, min.cells = 3, min.features = 500)
midbrain # 62044 cells

# Prepare metadata for the midbrain
dataset <- rep("Zeisel", 62044)
tissue <- rep("midbrain", 62044)

midbrain <- AddMetaData(midbrain, dataset, col.name = "dataset")
midbrain <- AddMetaData(midbrain, tissue, col.name = "tissue")

table(midbrain$tissue)

# Subset cells expressing Foxj1
midbrain_ecs <- subset(midbrain, subset = Foxj1 > 0)
midbrain_ecs


# Pons
data_dir <- "~/path/to/zeisel_integration/outs_by_tissue/pons/"
sample_IDs <- list.files("~/path/to/zeisel_integration/outs_by_tissue/pons/")

zeisel_pons_matrices <- sapply(sample_IDs, function(i){
  pons <- Read10X(file.path(data_dir,i,"outs/filtered_feature_bc_matrix/"))
  colnames(pons) <- paste(i, colnames(pons), sep = "_")
  pons
})

pons_data <- do.call("cbind", zeisel_pons_matrices)

dim(pons_data)

# Create Seurat object
pons <- CreateSeuratObject(counts = pons_data, min.cells = 3, min.features = 500)
pons # 57641 cells

# Prepare metadata for the amygdala
dataset <- rep("Zeisel", 57641)
tissue <- rep("pons", 57641)

pons <- AddMetaData(pons, dataset, col.name = "dataset")
pons <- AddMetaData(pons, tissue, col.name = "tissue")

table(pons$tissue)

# Subset cells expressing Foxj1
pons_ecs <- subset(pons, subset = Foxj1 > 0)
pons_ecs


# Spinal cord
data_dir <- "~/path/to/zeisel_integration/outs_by_tissue/spinal_cord/"
sample_IDs <- list.files("~/path/to/zeisel_integration/outs_by_tissue/spinal_cord/")

zeisel_spinal_cord_matrices <- sapply(sample_IDs, function(i){
  spinal_cord <- Read10X(file.path(data_dir,i,"outs/filtered_feature_bc_matrix/"))
  colnames(spinal_cord) <- paste(i, colnames(spinal_cord), sep = "_")
  spinal_cord
  
})

spinal_cord_data <- do.call("cbind", zeisel_spinal_cord_matrices)

dim(spinal_cord_data)

# Create Seurat object
spinal_cord <- CreateSeuratObject(counts = spinal_cord_data, min.cells = 3, min.features = 500)
spinal_cord # 25918 cells

# Prepare metadata for the spinal cord
dataset <- rep("Zeisel", 25918)
tissue <- rep("spinal_cord", 25918)

spinal_cord <- AddMetaData(spinal_cord, dataset, col.name = "dataset")
spinal_cord <- AddMetaData(spinal_cord, tissue, col.name = "tissue")

table(spinal_cord$tissue)


# Subset cells expressing Foxj1
spinal_cord_ecs <- subset(spinal_cord, subset = Foxj1 > 0)
spinal_cord_ecs


# Striatum
data_dir <- "~/path/to/zeisel_integration/outs_by_tissue/striatum/"
sample_IDs <- list.files("~/path/to/zeisel_integration/outs_by_tissue/striatum/")

zeisel_striatum_matrices <- sapply(sample_IDs, function(i){
  striatum <- Read10X(file.path(data_dir,i,"outs/filtered_feature_bc_matrix/"))
  colnames(striatum) <- paste(i, colnames(striatum), sep = "_")
  striatum
  
})

striatum_data <- do.call("cbind", zeisel_striatum_matrices)

dim(striatum_data)

# Create Seurat object
striatum <- CreateSeuratObject(counts = striatum_data, min.cells = 3, min.features = 500)
striatum # 36128 cells

# Prepare metadata for the amygdala
dataset <- rep("Zeisel", 36128)
tissue <- rep("striatum", 36128)

striatum <- AddMetaData(striatum, dataset, col.name = "dataset")
striatum <- AddMetaData(striatum, tissue, col.name = "tissue")

table(striatum$tissue)

# Subset cells expressing Foxj1
striatum_ecs <- subset(striatum, subset = Foxj1 > 0)
striatum_ecs


# Thalamus
data_dir <- "~/path/to/zeisel_integration/outs_by_tissue/thalamus/"
sample_IDs <- list.files("~/path/to/zeisel_integration/outs_by_tissue/thalamus/")

zeisel_thalamus_matrices <- sapply(sample_IDs, function(i){
  thalamus <- Read10X(file.path(data_dir,i,"outs/filtered_feature_bc_matrix/"))
  colnames(thalamus) <- paste(i, colnames(thalamus), sep = "_")
  thalamus
  
})

thalamus_data <- do.call("cbind", zeisel_thalamus_matrices)

dim(thalamus_data)

# Create Seurat object
thalamus <- CreateSeuratObject(counts = thalamus_data, min.cells = 3, min.features = 500)
thalamus # 52699 cells

# Prepare metadata for the amygdala
dataset <- rep("Zeisel", 52699)
tissue <- rep("thalamus", 52699)

thalamus <- AddMetaData(thalamus, dataset, col.name = "dataset")
thalamus <- AddMetaData(thalamus, tissue, col.name = "tissue")

table(thalamus$tissue)

# Subset cells expressing Foxj1
thalamus_ecs <- subset(thalamus, subset = Foxj1 > 0)
thalamus_ecs





## Merge all ependymal cell transcriptomes from Zeisel
zeisel_ecs <- merge(x = amygdala_ecs, y = c(hypothalamus_ecs, hippocampus_ecs, 
                                        medulla_ecs, midbrain_ecs, pons_ecs, spinal_cord_ecs,
                                        striatum_ecs, thalamus_ecs), merge.data = TRUE)

zeisel_ecs # OK!

names(zeisel_ecs[[]])
table(zeisel_ecs[["tissue"]])

zeisel_ecs[["pct_counts_mito"]] <- PercentageFeatureSet(zeisel_ecs, pattern = "^mt-")
zeisel_ecs[["pct_counts_ribo"]] <- PercentageFeatureSet(zeisel_ecs, pattern = "^Rp[sl][[:digit:]]")


## Processing the Zeisel dataset

# Load all cells expressing Foxj1 from Zeisel
zeisel_ecs <- readRDS(file = "~/path/to/zeisel_ecs.rds")

# Normalise
zeisel_ecs <- NormalizeData(zeisel_ecs, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
zeisel_ecs <- FindVariableFeatures(zeisel_ecs, selection.method = "vst", nfeatures = 2200)

# Scale data
all.genes <- rownames(zeisel_ecs)
zeisel_ecs <- ScaleData(zeisel_ecs, features = all.genes)

# PCA
zeisel_ecs <- RunPCA(zeisel_ecs, ndims = 1:30,
                         ndims.print = 1:5, nfeatures.print = 10)

# Plot PCs
DimPlot(object = zeisel_ecs, dims = c(1,2), reduction = "pca")
print(x = zeisel_ecs[["pca"]], dims = 1:4, nfeatures = 10)

VizDimLoadings(zeisel_ecs, dims = 1:2, reduction = "pca")


# What PCs explain most of the variance?
ElbowPlot(object = zeisel_ecs)

DimHeatmap(object = zeisel_ecs, dims = 13:21, balanced = TRUE, fast = FALSE, cells = 500) # 20 PCs

# UMAP
zeisel_ecs <- RunUMAP(object = zeisel_ecs, dims = 1:20, seed.use = 5,
                          n.neighbors = 30, learning.rate = 1)

DimPlot(object = zeisel_ecs, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.8) + NoLegend()
DimPlot(object = zeisel_ecs, reduction = "umap", group.by =  "tissue", label = TRUE, repel = TRUE, pt.size = 0.8)
DimPlot(object = zeisel_ecs, reduction = "umap", split.by =  "tissue", label = TRUE, repel = TRUE, pt.size = 0.8, ncol = 3)
DimPlot(object = zeisel_ecs, reduction = "umap", group.by =  "dataset", label = FALSE, repel = TRUE, pt.size = 0.8)


# Clustering
zeisel_ecs <- FindNeighbors(object = zeisel_ecs, dims = 1:20, k.param = 30)

zeisel_ecs <- FindClusters(object = zeisel_ecs, graph.name = NULL,
                               modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                               node.sizes = NULL, resolution = seq(0.1, 0.8, 0.1), algorithm = 1, n.start = 1000,
                               n.iter = 10, random.seed = 0, temp.file.location = NULL,
                               edge.file.name = NULL, verbose = TRUE)

# Clustering tree
library(clustree)

names(zeisel_ecs[[]])

clustree(zeisel_ecs, prefix = "RNA_snn_res.")


# Set cell identities
Idents(object = zeisel_ecs) <- zeisel_ecs[["RNA_snn_res.0.4"]]
table(x = Idents(zeisel_ecs))

DimPlot(object = zeisel_ecs, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.8)

# Plot UMAP
DimPlot(object = zeisel_ecs, reduction = "umap", group.by = "tissue", label = FALSE, repel = TRUE, pt.size = 0.8)
DimPlot(object = zeisel_ecs, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(object = zeisel_ecs, reduction = "umap", split.by = "tissue", pt.size = 0.8, ncol = 5)


# How many cells in each cluster?
table(Idents(zeisel_ecs))

# From which tissue?
table(Idents(zeisel_ecs), zeisel_ecs@meta.data$tissue)


FeaturePlot(zeisel_ecs, c("Foxj1", "Ecrg4", "Mia", "Ccdc153", "Tmem212", "Rarres2", "Ccnd2"),
            order = TRUE, cols = c("lightgrey", "red"), pt.size = 0.8)


# Differential expression analysis
# Find markers
markers_zeisel_ecs <- FindAllMarkers(zeisel_ecs, only.pos = TRUE, test.use = "wilcox",
                                         logfc.threshold = 0.25, min.pct = 0.15, pseudocount.use = 1)

c01246_vs_c3 <- FindMarkers(zeisel_ecs, ident.1 = c(0,1,2,6,4),
                        ident.2 = 3, only.pos = TRUE,
                        logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)

c3_vs_c10 <- FindMarkers(zeisel_ecs, ident.1 = 3,
                        ident.2 = 10, only.pos = TRUE,
                        logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)

## Filtering out non-ependymal cells

# Subsettingependymal cells (clusters 3 and 10)
# Cells in clusters 3 and 10 express the ependymal markers Ccdc153, Tmem212, Rarres2, Ecrg4, Mia
zeisel_ecs_filtered <- subset(zeisel_ecs, idents = c(3,10))

zeisel_ecs_filtered

# Normalise
zeisel_ecs_filtered <- NormalizeData(zeisel_ecs_filtered, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
zeisel_ecs_filtered <- FindVariableFeatures(zeisel_ecs_filtered, selection.method = "vst", nfeatures = 2200)

# Scale data
all.genes <- rownames(zeisel_ecs_filtered)
zeisel_ecs_filtered <- ScaleData(zeisel_ecs_filtered, features = all.genes)


# PCA
zeisel_ecs_filtered <- RunPCA(zeisel_ecs_filtered, ndims = 1:30,
                     ndims.print = 1:5, nfeatures.print = 10)

# Plot PCs
DimPlot(object = zeisel_ecs_filtered, dims = c(1,3), reduction = "pca")
print(x = zeisel_ecs_filtered[["pca"]], dims = 1:6, nfeatures = 10)

VizDimLoadings(zeisel_ecs_filtered, dims = 1:2, reduction = "pca")


# What PCs explain most of the variance?
ElbowPlot(object = zeisel_ecs_filtered)

DimHeatmap(object = zeisel_ecs_filtered, dims = 13:21, balanced = TRUE, fast = FALSE, cells = 500) # 20 PCs

# UMAP
zeisel_ecs_filtered <- RunUMAP(object = zeisel_ecs_filtered, dims = 1:20, seed.use = 1,
                               n.neighbors = 10, learning.rate = 100)

DimPlot(object = zeisel_ecs_filtered, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 1)
DimPlot(object = zeisel_ecs_filtered, reduction = "umap", group.by =  "tissue", label = FALSE, repel = TRUE, pt.size = 1)

FeaturePlot(zeisel_ecs_filtered, c("Krt15", "Arx", "Sntn", "Ucma", "Bhmt", "Vtn", "Mia"),
            order = TRUE, cols = c("lightgrey", "red"), pt.size = 1.2, ncol = 3)

# Clustering
zeisel_ecs_filtered <- FindNeighbors(object = zeisel_ecs_filtered, dims = 1:20, k.param = 10)

zeisel_ecs_filtered <- FindClusters(object = zeisel_ecs_filtered, graph.name = NULL,
                           modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                           node.sizes = NULL, resolution = seq(0.3, 1, 0.1), algorithm = 1, n.start = 1000,
                           n.iter = 10, random.seed = 0, temp.file.location = NULL,
                           edge.file.name = NULL, verbose = TRUE)

# Clustering tree
clustree(zeisel_ecs_filtered, prefix = "RNA_snn_res.")

# Set cell identities
Idents(object = zeisel_ecs_filtered) <- zeisel_ecs_filtered[["RNA_snn_res.0.7"]]
table(x = Idents(zeisel_ecs_filtered))

DimPlot(object = zeisel_ecs_filtered, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)


# Plot UMAP
DimPlot(object = zeisel_ecs_filtered, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1) + NoLegend()
DimPlot(object = zeisel_ecs_filtered, reduction = "umap", split.by = "tissue", pt.size = 1, ncol = 3, label = TRUE)
DimPlot(object = zeisel_ecs_filtered, reduction = "umap", group.by = "tissue", pt.size = 1, label = FALSE)

FeaturePlot(zeisel_ecs_filtered, c("Sntn"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(zeisel_ecs_filtered, c("Six3"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 1)


# Differential expression analysis
# Find markers
markers_zeisel_ecs_filtered <- FindAllMarkers(zeisel_ecs_filtered, only.pos = TRUE, test.use = "wilcox",
                                     logfc.threshold = 0.25, min.pct = 0.15, pseudocount.use = 1)

## Filtering out non-ependymal cells
# Removing cells expressing  oligodendrocyte marker genes (cluster 6) and re-running the analysis
zeisel_ecs_filtered2 <- subset(zeisel_ecs_filtered, idents = c(0,1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17))

zeisel_ecs_filtered2

# Normalise
zeisel_ecs_filtered2 <- NormalizeData(zeisel_ecs_filtered2, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
zeisel_ecs_filtered2 <- FindVariableFeatures(zeisel_ecs_filtered2, selection.method = "vst", nfeatures = 2200)

# Scale data
all.genes <- rownames(zeisel_ecs_filtered2)
zeisel_ecs_filtered2 <- ScaleData(zeisel_ecs_filtered2, features = all.genes)

# PCA
zeisel_ecs_filtered2 <- RunPCA(zeisel_ecs_filtered2, ndims = 1:30,
                              ndims.print = 1:5, nfeatures.print = 10)

# Plot PCs
DimPlot(object = zeisel_ecs_filtered2, dims = c(1,2), reduction = "pca", pt.size = 1)
print(x = zeisel_ecs_filtered2[["pca"]], dims = 1:4, nfeatures = 10)

VizDimLoadings(zeisel_ecs_filtered2, dims = 1:2, reduction = "pca")


# What PCs explain most of the variance?
ElbowPlot(object = zeisel_ecs_filtered2)

DimHeatmap(object = zeisel_ecs_filtered2, dims = 13:24, balanced = TRUE, fast = FALSE, cells = 500)

# UMAP
zeisel_ecs_filtered2 <- RunUMAP(object = zeisel_ecs_filtered2, dims = 1:15, seed.use = 1,
                                n.neighbors = 10, learning.rate = 100, min.dist = 0.1)


DimPlot(object = zeisel_ecs_filtered2, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(object = zeisel_ecs_filtered2, reduction = "umap", group.by =  "tissue", label = FALSE, repel = TRUE, pt.size = 1)

FeaturePlot(zeisel_ecs_filtered2, c("Arx", "Sntn", "Ucma", "Bhmt"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 1, ncol = 2)

# Clustering
zeisel_ecs_filtered2 <- FindNeighbors(object = zeisel_ecs_filtered2, dims = 1:15, k.param = 10)

zeisel_ecs_filtered2 <- FindClusters(object = zeisel_ecs_filtered2, graph.name = NULL,
                                    modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                                    node.sizes = NULL, resolution = seq(0.2, 1, 0.1), algorithm = 1, n.start = 1000,
                                    n.iter = 10, random.seed = 0, temp.file.location = NULL,
                                    edge.file.name = NULL, verbose = TRUE)

zeisel_ecs_filtered2 <- FindClusters(object = zeisel_ecs_filtered2, graph.name = NULL,
                                    modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                                    node.sizes = NULL, resolution = 1.5, algorithm = 1, n.start = 1000,
                                    n.iter = 10, random.seed = 0, temp.file.location = NULL,
                                    edge.file.name = NULL, verbose = TRUE)

# Clustering tree
clustree(zeisel_ecs_filtered2, prefix = "RNA_snn_res.")


# Set cell identities
Idents(object = zeisel_ecs_filtered2) <- zeisel_ecs_filtered2[["RNA_snn_res.0.8"]]
table(x = Idents(zeisel_ecs_filtered2))

DimPlot(object = zeisel_ecs_filtered2, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

# Plot UMAP
DimPlot(object = zeisel_ecs_filtered2, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(object = zeisel_ecs_filtered2, reduction = "umap", split.by = "tissue", pt.size = 1, ncol = 3, label = TRUE)
DimPlot(object = zeisel_ecs_filtered2, reduction = "umap", group.by  = "orig.ident", pt.size = 1, label = FALSE)


mypal <- DiscretePalette(n = 35, palette = "polychrome")
DimPlot(object = zeisel_ecs, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1, cols = mypal)
DimPlot(object = zeisel_ecs, reduction = "umap", split.by = "tissue", pt.size = 1, ncol = 3, label = TRUE, repel = TRUE, cols = mypal)

# How many cells in each cluster?
table(Idents(zeisel_ecs))

# From which tissue?
table(Idents(zeisel_ecs), zeisel_ecs@meta.data$tissue)

FeaturePlot(zeisel_ecs, c("Dpysl2", "Dsg2", "Abca1", "Mtpn"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(zeisel_ecs, c("Cdh26", "Cdh12"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(zeisel_ecs, c("Gfap"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(zeisel_ecs, c("Six3"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 1)


# Differential expression analysis
# Find markers
markers_zeisel_ecs <- FindAllMarkers(zeisel_ecs, only.pos = TRUE, test.use = "wilcox",
                                              logfc.threshold = 0.25, min.pct = 0.15, pseudocount.use = 1)

## Filtering out non-ependymal cells
# Removing cells expressing  oligodendrocyte (cluster 17) and neuronal (cluster 16) genes
zeisel_ecs_filtered3 <- subset(zeisel_ecs, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,18))

zeisel_ecs_filtered3

# Normalise
zeisel_ecs_filtered3 <- NormalizeData(zeisel_ecs_filtered3, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
zeisel_ecs_filtered3 <- FindVariableFeatures(zeisel_ecs_filtered3, selection.method = "vst", nfeatures = 2200)

# Scale data
all.genes <- rownames(zeisel_ecs_filtered3)
zeisel_ecs_filtered3 <- ScaleData(zeisel_ecs_filtered3, features = all.genes)

# PCA
zeisel_ecs_filtered3 <- RunPCA(zeisel_ecs_filtered3, ndims = 1:30,
                               ndims.print = 1:5, nfeatures.print = 10)

# Plot PCs
DimPlot(object = zeisel_ecs_filtered3, dims = c(1,2), reduction = "pca", pt.size = 1)
print(x = zeisel_ecs_filtered3[["pca"]], dims = 1:4, nfeatures = 10)

VizDimLoadings(zeisel_ecs_filtered3, dims = 1:2, reduction = "pca")


# What PCs explain most of the variance?
ElbowPlot(object = zeisel_ecs_filtered3)

DimHeatmap(object = zeisel_ecs_filtered3, dims = 13:24, balanced = TRUE, fast = FALSE, cells = 500)

# UMAP
zeisel_ecs_filtered3 <- RunUMAP(object = zeisel_ecs_filtered3, dims = 1:15, seed.use = 6,
                                n.neighbors = 10, learning.rate = 100, min.dist = 0.1)

DimPlot(object = zeisel_ecs_filtered3, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(object = zeisel_ecs_filtered3, reduction = "umap", group.by =  "tissue", label = FALSE, repel = TRUE, pt.size = 1)


# Clustering
zeisel_ecs_filtered3 <- FindNeighbors(object = zeisel_ecs_filtered3, dims = 1:15, k.param = 10)

zeisel_ecs_filtered3 <- FindClusters(object = zeisel_ecs_filtered3, graph.name = NULL,
                                     modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                                     node.sizes = NULL, resolution = seq(0.2, 1, 0.1), algorithm = 1, n.start = 1000,
                                     n.iter = 10, random.seed = 0, temp.file.location = NULL,
                                     edge.file.name = NULL, verbose = TRUE)


# Clustering tree
clustree(zeisel_ecs_filtered3, prefix = "RNA_snn_res.")
clustree(zeisel_ecs, prefix = "RNA_snn_res.")

# Set cell identities
Idents(object = zeisel_ecs_filtered3) <- zeisel_ecs_filtered3[["RNA_snn_res.0.8"]]
table(x = Idents(zeisel_ecs_filtered3))

DimPlot(object = zeisel_ecs_filtered3, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

# Plot UMAP
DimPlot(object = zeisel_ecs_filtered3, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(object = zeisel_ecs_filtered3, reduction = "umap", split.by = "tissue", pt.size = 1, ncol = 3, label = TRUE)
DimPlot(object = zeisel_ecs_filtered3, reduction = "umap", group.by  = "orig.ident", pt.size = 1, label = FALSE)

mypal <- DiscretePalette(n = 35, palette = "polychrome")
DimPlot(object = zeisel_ecs_filtered3, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1, cols = mypal)
DimPlot(object = zeisel_ecs_filtered3, reduction = "umap", split.by = "tissue", pt.size = 1, ncol = 3, label = TRUE, repel = TRUE, cols = mypal)

# How many cells in each cluster?
table(Idents(zeisel_ecs_filtered3))

# From which tissue?
table(Idents(zeisel_ecs_filtered3), zeisel_ecs_filtered3@meta.data$tissue)

# Differential expression analysis
# Find markers
markers_zeisel_ecs <- FindAllMarkers(zeisel_ecs_filtered3, only.pos = TRUE, test.use = "wilcox",
                                     logfc.threshold = 0.25, min.pct = 0.15, pseudocount.use = 1)

# Save cluster ID to celltype
table(Idents(zeisel_ecs_filtered3))
zeisel_ecs_filtered3$celltype <- Idents(zeisel_ecs_filtered3)


# Exploring cell clusters in Zeisel
c0_vs_c1 <- FindMarkers(zeisel_ecs, ident.1 = 0,
                        ident.2 = 1, only.pos = TRUE,
                        logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1) # no significantly different
c1_vs_c0 <- FindMarkers(zeisel_ecs, ident.1 = 1,
                        ident.2 = 0, only.pos = TRUE,
                        logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1) # a few genes differentially expressed in cluster 1 but mostly related to % cells expressing


## Filtering out non-ependymal cells
# Removing choroid plexus cells and cells expressing oligodendrocyte genes and high fraction of mitochondrial genes
zeisel_ecs_filtered4 <- subset(zeisel_ecs, idents = c(0,1,2,3,4,5,7,8,9,10,11,12,14,15))

zeisel_ecs_filtered4

# Normalise
zeisel_ecs_filtered4 <- NormalizeData(zeisel_ecs_filtered4, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
zeisel_ecs_filtered4 <- FindVariableFeatures(zeisel_ecs_filtered4, selection.method = "vst", nfeatures = 2200)

# Scale data
all.genes <- rownames(zeisel_ecs_filtered4)
zeisel_ecs_filtered4 <- ScaleData(zeisel_ecs_filtered4, features = all.genes)


# PCA
zeisel_ecs_filtered4 <- RunPCA(zeisel_ecs_filtered4, ndims = 1:30,
                               ndims.print = 1:5, nfeatures.print = 10)


# What PCs explain most of the variance?
ElbowPlot(object = zeisel_ecs_filtered4)

DimHeatmap(object = zeisel_ecs_filtered4, dims = 13:24, balanced = TRUE, fast = FALSE, cells = 500)

# UMAP
zeisel_ecs_filtered4 <- RunUMAP(object = zeisel_ecs_filtered4, dims = 1:15, seed.use = 6,
                                n.neighbors = 10, learning.rate = 100, min.dist = 0.1)

DimPlot(object = zeisel_ecs_filtered4, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(object = zeisel_ecs_filtered4, reduction = "umap", group.by =  "tissue", label = FALSE, repel = TRUE, pt.size = 1)


# Clustering
zeisel_ecs_filtered4 <- FindNeighbors(object = zeisel_ecs_filtered4, dims = 1:15, k.param = 10)

zeisel_ecs_filtered4 <- FindClusters(object = zeisel_ecs_filtered4, graph.name = NULL,
                                     modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                                     node.sizes = NULL, resolution = seq(0.2, 1, 0.1), algorithm = 1, n.start = 1000,
                                     n.iter = 10, random.seed = 0, temp.file.location = NULL,
                                     edge.file.name = NULL, verbose = TRUE)


# Clustering tree
clustree(zeisel_ecs_filtered4, prefix = "RNA_snn_res.")

# Set cell identities
Idents(object = zeisel_ecs_filtered4) <- zeisel_ecs_filtered4[["RNA_snn_res.0.9"]]
table(x = Idents(zeisel_ecs_filtered4))

DimPlot(object = zeisel_ecs_filtered4, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1, cols = mypal)


# Plot UMAP
DimPlot(object = zeisel_ecs_filtered4, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(object = zeisel_ecs_filtered4, reduction = "umap", split.by = "tissue", pt.size = 1, ncol = 5, label = TRUE, cols = mypal)
DimPlot(object = zeisel_ecs_filtered4, reduction = "umap", group.by  = "orig.ident", pt.size = 1, label = FALSE)


# How many cells in each cluster?
table(Idents(zeisel_ecs_filtered4))

# From which tissue?
table(Idents(zeisel_ecs_filtered4), zeisel_ecs_filtered4@meta.data$tissue)

FeaturePlot(zeisel_ecs_filtered4, c("Cdh26", "Cdh12"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(zeisel_ecs_filtered4, c("Gfap"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 1)



# Differential expression analysis
# Find markers
markers_zeisel_ecs <- FindAllMarkers(zeisel_ecs_filtered4, only.pos = TRUE, test.use = "wilcox",
                                     logfc.threshold = 0.25, min.pct = 0.15, pseudocount.use = 1)

# Save cluster ID to celltype
table(Idents(zeisel_ecs_filtered4))
zeisel_ecs_filtered4$celltype <- Idents(zeisel_ecs_filtered4)




## Integrating Zeisel and Rodrigo Albors datasets

# Loading the datasets
rodrigo <- readRDS(file = "~/path/to/young_ecs_subtypes_without_astrolike_cells.rds") # from 01_central_canal_cells_young_10x.R
zeisel <- readRDS(file = "~/path/to/zeisel_ecs_filtered_04_clustering.rds")

# keeping common genes
is.common <- rownames(rodrigo) %in% rownames(zeisel)
rodrigo <- rodrigo[is.common,]

zeisel <- zeisel[match(rownames(rodrigo), rownames(zeisel)),]

identical(rownames(zeisel), rownames(rodrigo))


names(rodrigo[[]])
rodrigo
DimPlot(rodrigo, label = TRUE, repel = TRUE) + NoLegend()

DimPlot(zeisel, label = TRUE, repel = TRUE) + NoLegend()
DimPlot(object = zeisel, reduction = "umap", label = F, repel = TRUE, pt.size = 0.8, group.by = "tissue", cols = pal) # Figure 4

# Add dataset and tissue info to metadata
dataset <- rep("Rodrigo", 1901)
tissue <- rep("spinal_cord", 1901)

rodrigo <- AddMetaData(rodrigo, dataset, col.name = "dataset")
rodrigo <- AddMetaData(rodrigo, tissue, col.name = "tissue")

table(rodrigo$tissue)

rodrigo

# Merging data
ecs <- merge(x = zeisel, y = rodrigo, merge.data = TRUE)

ecs # OK!

names(ecs[[]])

# Removing earlier clustering info
columns.to.remove <- ecs@meta.data %>% 
  names() %>% 
  .[stringr::str_detect(., 'integrated')]

for(i in columns.to.remove) {
  ecs[[i]] <- NULL
}

names(ecs[[]])

columns.to.remove <- ecs@meta.data %>% 
  names() %>% 
  .[stringr::str_detect(., 'RNA_snn')]

for(i in columns.to.remove) {
  ecs[[i]] <- NULL
}

names(ecs[[]])

# Split objects
ecs_list <- SplitObject(ecs, split.by = "dataset")

ecs_list

# Normalisation
for (i in 1:length(ecs_list)) {
  ecs_list[[i]] <- NormalizeData(ecs_list[[i]], verbose = FALSE, normalization.method = "LogNormalize",
                                 scale.factor = 10000)
  ecs_list[[i]] <- FindVariableFeatures(ecs_list[[i]], selection.method = "vst", 
                                        nfeatures = 2500, verbose = FALSE)
}


# Selecting integration features
ecs_features <- SelectIntegrationFeatures(object.list = ecs_list, nfeatures = 2000)

# Identify anchors and integrate the data sets
ecs_anchors <- FindIntegrationAnchors(object.list = ecs_list, normalization.method = "LogNormalize",
                                      anchor.features =  ecs_features)
ecs_integrated <- IntegrateData(anchorset = ecs_anchors, normalization.method = "LogNormalize")

DefaultAssay(ecs_integrated)

# Running the standard Seurat workflow for visualization
DefaultAssay(ecs_integrated) <- "integrated"
ecs_integrated <- ScaleData(ecs_integrated, verbose = FALSE)
ecs_integrated <- RunPCA(ecs_integrated, ndims = 1:30,
                         ndims.print = 1:5, nfeatures.print = 10)

# Plot PCs
DimPlot(object = ecs_integrated, dims = c(1,2), reduction = "pca")
print(x = ecs_integrated[["pca"]], dims = 1, nfeatures = 10)

# What PCs explain most of the variance?
ElbowPlot(object = ecs_integrated)

DimHeatmap(object = ecs_integrated, dims = 13:21, balanced = TRUE, fast = FALSE, cells = 500)
DimHeatmap(object = ecs_integrated, dims = 19:27, balanced = TRUE, fast = FALSE, cells = 500)

# UMAP
ecs_integrated <- RunUMAP(object = ecs_integrated, dims = 1:19, seed.use = 0,
                          n.neighbors = 10, learning.rate = 1, min.dist = 0.3) 

DimPlot(object = ecs_integrated, reduction = "umap", label = T, repel = TRUE, pt.size = 0.8) + NoLegend()

DimPlot(object = ecs_integrated, reduction = "umap", split.by = "dataset", label = FALSE, repel = TRUE, pt.size = 0.8)
DimPlot(object = ecs_integrated, reduction = "umap", split.by =  "tissue", label = TRUE, repel = TRUE, 
        pt.size = 0.8, ncol = 3, cols = mypal, group.by = "dataset") + NoLegend()
DimPlot(object = ecs_integrated, reduction = "umap", split.by =  "tissue", label = TRUE, repel = TRUE, 
        pt.size = 0.8, ncol = 3, cols = mypal) + NoLegend()
DimPlot(object = ecs_integrated, reduction = "umap", group.by =  "dataset", label = FALSE, repel = TRUE, pt.size = 0.8)
DimPlot(object = ecs_integrated, reduction = "umap", group.by =  "tissue", label = FALSE, repel = TRUE, pt.size = 0.8)

table(ecs_integrated[["dataset"]])


# Clustering
DefaultAssay(ecs_integrated) <- "integrated"

ecs_integrated <- FindNeighbors(object = ecs_integrated, dims = 1:19, k.param = 10)

ecs_integrated <- FindClusters(object = ecs_integrated, graph.name = NULL,
                               modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                               node.sizes = NULL, resolution = seq(0.2, 0.8, 0.1), algorithm = 1, n.start = 1000,
                               n.iter = 10, random.seed = 0, temp.file.location = NULL,
                               edge.file.name = NULL, verbose = TRUE)

ecs_integrated <- FindClusters(object = ecs_integrated, graph.name = NULL,
                               modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                               node.sizes = NULL, resolution = 1.5, algorithm = 1, n.start = 1000,
                               n.iter = 10, random.seed = 0, temp.file.location = NULL,
                               edge.file.name = NULL, verbose = TRUE)

ecs_integrated <- FindClusters(object = ecs_integrated, graph.name = NULL,
                               modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                               node.sizes = NULL, resolution = 0.1, algorithm = 1, n.start = 1000,
                               n.iter = 10, random.seed = 0, temp.file.location = NULL,
                               edge.file.name = NULL, verbose = TRUE)

ecs_integrated <- FindClusters(object = ecs_integrated, graph.name = NULL,
                               modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                               node.sizes = NULL, resolution = 1, algorithm = 1, n.start = 1000,
                               n.iter = 10, random.seed = 0, temp.file.location = NULL,
                               edge.file.name = NULL, verbose = TRUE)

ecs_integrated <- FindClusters(object = ecs_integrated, graph.name = NULL,
                               modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                               node.sizes = NULL, resolution = 2, algorithm = 1, n.start = 1000,
                               n.iter = 10, random.seed = 0, temp.file.location = NULL,
                               edge.file.name = NULL, verbose = TRUE)


# Clustering tree
names(ecs_integrated[[]])
clustree(ecs_integrated, prefix = "integrated_snn_res.")
clustree(ecs_integrated, prefix = "integrated_snn_res.", node_colour_aggr = "mean", node_colour = "Sntn") + scale_color_viridis()

# Set cell identities
mypal <- DiscretePalette(n = 35, palette = "polychrome")
Idents(object = ecs_integrated) <- ecs_integrated[["integrated_snn_res.1.5"]]

table(x = Idents(ecs_integrated))

# Plot UMAP
DimPlot(object = ecs_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1, cols = mypal)
DimPlot(object = ecs_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(object = ecs_integrated, reduction = "umap", split.by = "tissue", pt.size = 0.8, ncol = 3, cols = mypal)
DimPlot(object = ecs_integrated, reduction = "umap", split.by = "tissue", group.by = "tissue", pt.size = 0.8, ncol = 5)
DimPlot(object = ecs_integrated, reduction = "umap", split.by = "tissue", pt.size = 0.8, ncol = 5)
DimPlot(object = ecs_integrated, reduction = "umap", split.by = "dataset", pt.size = 1, group.by = "tissue")
pal <- brewer.pal(n=9, name = 'Set3')
DimPlot(object = ecs_integrated, reduction = "umap", group.by = "tissue", pt.size = 0.6, cols = pal)
DimPlot(object = ecs_integrated, reduction = "umap", group.by = "dataset", pt.size = 0.6)

# Normalise RNA counts for visualisation purposes
DefaultAssay(ecs_integrated) <- "RNA" # use RNA assay for visualisation and DE analysis

ecs_integrated <- NormalizeData(ecs_integrated, verbose = FALSE)
FeaturePlot(ecs_integrated, c("Sntn", "Arx", "Bhmt", "Cdh26"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 0.8)
FeaturePlot(ecs_integrated, c("Sntn"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 0.8, split.by = "tissue", ncol = 3)
FeaturePlot(ecs_integrated, c("Mki67"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 0.8, split.by = "tissue", ncol = 3)

# Differential expression analysis
DefaultAssay(ecs_integrated) <- "RNA" # use RNA assay for visualisation and DE analysis

# Find markers
markers_ecs_integrated <- FindAllMarkers(ecs_integrated, only.pos = TRUE, test.use = "wilcox",
                                         logfc.threshold = 0.25, min.pct = 0.25, pseudocount.use = 1)

FeaturePlot(ecs_integrated, c("Tagln"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 0.8)

tissue_levels <- c("hippocampus", "striatum", "amygdala", "thalamus", "hypothalamus", "midbrain", "pons", "medulla", "spinal_cord")

ecs_integrated$tissue <- factor(ecs_integrated$tissue, levels = tissue_levels)

DimPlot(object = ecs_integrated, reduction = "umap", split.by = "tissue", group.by = "dataset", pt.size = 0.6, ncol = 4)
DimPlot(object = ecs_integrated, reduction = "umap", group.by = "tissue", pt.size = 0.6) +
  scale_color_brewer(type = 'qual', palette = 'Set3', direction = -1) # Figure S5B


# From which tissue?
table(Idents(ecs_integrated), ecs_integrated@meta.data$tissue)

c0_vs_1 <- FindMarkers(ecs_integrated, ident.1 = 0,
                       ident.2 = 1, only.pos = TRUE,
                       logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
c1_vs_c0 <- FindMarkers(ecs_integrated, ident.1 = 1,
                        ident.2 = 0, only.pos = TRUE,
                        logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)

c10_vs_c8 <- FindMarkers(ecs_integrated, ident.1 = 10,
                         ident.2 = 8, only.pos = TRUE,
                         logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
c8_vs_c10 <- FindMarkers(ecs_integrated, ident.1 = 8,
                         ident.2 = 10, only.pos = TRUE,
                         logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
c11_vs_c3 <- FindMarkers(ecs_integrated, ident.1 = 11,
                         ident.2 = 3, only.pos = TRUE,
                         logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
c3_vs_c11 <- FindMarkers(ecs_integrated, ident.1 = 3,
                         ident.2 = 11, only.pos = TRUE,
                         logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
c17_vs_c0 <- FindMarkers(ecs_integrated, ident.1 = 17,
                         ident.2 = 0, only.pos = TRUE,
                         logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
c0_vs_c17 <- FindMarkers(ecs_integrated, ident.1 = 0,
                         ident.2 = 17, only.pos = TRUE,
                         logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
c3_vs_c7 <- FindMarkers(ecs_integrated, ident.1 = 3,
                        ident.2 = 7, only.pos = TRUE,
                        logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
c0_vs_c17 <- FindMarkers(ecs_integrated, ident.1 = 0,
                         ident.2 = 17, only.pos = TRUE,
                         logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)

FeaturePlot(ecs_integrated, c("Wkn", "Hydin"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 0.8)
FeaturePlot(ecs_integrated, c("Cdh26", "Zic1", "Cldn11"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 0.8)
FeaturePlot(ecs_integrated, c("Tnnt3", "Cfap126", "Mif1"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 0.8)

# Merge clusters that do not seem biologically relevant
Idents(ecs_integrated) <- ecs_integrated$integrated_snn_res.1

DimPlot(object = ecs_integrated, reduction = "umap", label = F, repel = TRUE, pt.size = 0.6)

# Merging from clustering resolution 1.0
ecs_integrated <- RenameIdents(ecs_integrated, `0` = "0", `1` = "0", `2` = "1",
                               `3` = "2", `4` = "3", `5` = "4",
                               `6` = "2", `7` = "2", `8` = "1", `9` = "5",
                               `10` = "1", `11` = "6", `12` = "7", `13` = "8", `14` = "9",
                               `15` = "10", `16` = "11", `17` = "12")


DimPlot(object = ecs_integrated, reduction = "umap", label = F, repel = T, pt.size = 0.6)
DimPlot(object = ecs_integrated, reduction = "umap", label = F, repel = TRUE, pt.size = 0.6, split.by = "tissue", ncol = 3)

ecs_integrated$merged_clusters <- Idents(ecs_integrated)

table(ecs_integrated$merged_clusters)

names(ecs_integrated[[]])

Idents(ecs_integrated) <- ecs_integrated$tissue
tissuepal <- brewer.pal(n=9, name = 'Set1')

DimPlot(object = ecs_integrated, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 0.8, cols = pal)

table(Idents(ecs_integrated))

# Assign anatomical region
ecs_integrated <- RenameIdents(ecs_integrated, `amygdala` = "forebrain", `hypothalamus` = "forebrain", `hippocampus` = "forebrain",
                               `medulla` = "hindbrain", `midbrain` = "midbrain", `pons` = "hindbrain",
                               `spinal_cord` = "spinal_cord", `striatum` = "forebrain", `thalamus` = "forebrain")

ecs_integrated$tissue_region <- Idents(ecs_integrated)
Idents(ecs_integrated) <- ecs_integrated$tissue_region
Idents(ecs_integrated) <- ecs_integrated$tissue

DimPlot(object = ecs_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.6)

# Differential expression analysis
Idents(ecs_integrated) <- ecs_integrated$merged_clusters
DimPlot(object = ecs_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(object = ecs_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1, split.by = "dataset")

DefaultAssay(ecs_integrated) <- "RNA"

# Find markers
DefaultAssay(ecs_integrated) <- "RNA"
markers_ecs_integrated <- FindAllMarkers(ecs_integrated, only.pos = TRUE, test.use = "wilcox",
                                         logfc.threshold = 0.25, min.pct = 0.25, pseudocount.use = 1)

table(Idents(ecs_integrated))

## Filtering out non-ependymal cells (microglia and mixed astrocytes and ependymal cells)
Idents(ecs_integrated) <- ecs_integrated$integrated_snn_res.1
DimPlot(object = ecs_integrated, reduction = "umap", label = T, repel = TRUE, pt.size = 0.6, cols = mypal)

tissue_levels <- c("hippocampus", "striatum", "amygdala", "thalamus", "hypothalamus", "midbrain", "pons", "medulla", "spinal_cord")
ecs_integrated$tissue <- factor(ecs_integrated$tissue, levels = tissue_levels)


table(Idents(ecs_integrated))
ecs_integrated <- subset(x= ecs_integrated, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,15,16))

DefaultAssay(ecs_integrated) <- "RNA"

# Split objects
ecs_integrated_list <- SplitObject(ecs_integrated, split.by = "dataset")

ecs_integrated_list

# Normalisation
for (i in 1:length(ecs_integrated_list)) {
  ecs_integrated_list[[i]] <- NormalizeData(ecs_integrated_list[[i]], verbose = FALSE, normalization.method = "LogNormalize",
                                            scale.factor = 10000)
  ecs_integrated_list[[i]] <- FindVariableFeatures(ecs_integrated_list[[i]], selection.method = "vst",
                                                   nfeatures = 2500, verbose = FALSE)
}

# Re-running the integration steps

# Find anchors
ecs_integrated_anchors <- FindIntegrationAnchors(object.list = ecs_integrated_list, dims = 1:30)

ecs_integrated <- IntegrateData(anchorset = ecs_integrated_anchors, dims = 1:30)

ecs_integrated

DefaultAssay(ecs_integrated)

# Running the standard workflow for visualization and clustering
DefaultAssay(ecs_integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
ecs_integrated <- ScaleData(ecs_integrated, verbose = FALSE)
ecs_integrated <- RunPCA(ecs_integrated, npcs = 30, verbose = FALSE)

# Plot PCs
DimPlot(object = ecs_integrated, dims = c(1,2), reduction = "pca", group.by = "tech")
print(x = ecs_integrated[["pca"]], dims = 1, nfeatures = 10)


# What PCs explain most of the variance?
ElbowPlot(object = ecs_integrated)

DimHeatmap(object = ecs_integrated, dims = 10:18, balanced = TRUE, fast = FALSE, cells = 500)
DimHeatmap(object = ecs_integrated, dims = 1:9, balanced = TRUE, fast = FALSE, cells = 500)

# UMAP
ecs_integrated <- RunUMAP(object = ecs_integrated, dims = 1:22, seed.use = 13, learning.rate = 1,
                          n.neighbors = 10, spread = 1, min.dist = 0.3)

DimPlot(object = ecs_integrated, reduction = "umap", label = F, repel = TRUE, pt.size = 0.8, cols = mypal) + NoLegend()


# Removing earlier clustering info
columns.to.remove <- ecs_integrated@meta.data %>%
  names() %>%
  .[stringr::str_detect(., 'integrated')]

for(i in columns.to.remove) {
  ecs_integrated[[i]] <- NULL
}

names(ecs_integrated[[]])


# Clustering
DefaultAssay(ecs_integrated) <- "integrated"

names(ecs_integrated[[]])

ecs_integrated <- FindNeighbors(object = ecs_integrated, dims = 1:22, k.param = 10)

ecs_integrated <- FindClusters(object = ecs_integrated, graph.name = NULL,
                               modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                               node.sizes = NULL, resolution = seq(0.2, 0.8, 0.1), algorithm = 1, n.start = 1000,
                               n.iter = 10, random.seed = 0, temp.file.location = NULL,
                               edge.file.name = NULL, verbose = TRUE)

ecs_integrated <- FindClusters(object = ecs_integrated, graph.name = NULL,
                               modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                               node.sizes = NULL, resolution = 1, algorithm = 1, n.start = 1000,
                               n.iter = 10, random.seed = 0, temp.file.location = NULL,
                               edge.file.name = NULL, verbose = TRUE)

ecs_integrated <- FindClusters(object = ecs_integrated, graph.name = NULL,
                               modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                               node.sizes = NULL, resolution = 1.5, algorithm = 1, n.start = 1000,
                               n.iter = 10, random.seed = 0, temp.file.location = NULL,
                               edge.file.name = NULL, verbose = TRUE)


# Clustering tree
clustree(ecs_integrated, prefix = "integrated_snn_res.")
clustree(ecs_integrated, prefix = "integrated_snn_res.", node_colour_aggr = "mean", node_colour = "Sntn") + scale_color_viridis()

# Set cell identities
Idents(object = ecs_integrated) <- ecs_integrated[["integrated_snn_res.0.8"]]

table(x = Idents(ecs_integrated))

DimPlot(object = ecs_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.8)

# Plot UMAP
DimPlot(object = ecs_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

pal <- brewer.pal(n=9, name = 'Set3')
DimPlot(object = ecs_integrated, reduction = "umap", group.by = "tissue", pt.size = 0.6, cols = pal)# Figure S5B
DimPlot(object = ecs_integrated, reduction = "umap", group.by = "dataset", pt.size = 0.6)


# Normalise RNA counts for visualisation purposes
DefaultAssay(ecs_integrated) <- "RNA"

ecs_integrated <- NormalizeData(ecs_integrated, verbose = FALSE)
FeaturePlot(ecs_integrated, c("Sntn", "Arx", "Bhmt", "Cdh26"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 0.8)

# Differential expression analysis
# Find markers
markers_ecs_integrated <- FindAllMarkers(ecs_integrated, only.pos = TRUE, test.use = "wilcox",
                                         logfc.threshold = 0.25, min.pct = 0.25, pseudocount.use = 1)

FeaturePlot(ecs_integrated, c("Cftr"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 0.8)
DimPlot(object = ecs_integrated, reduction = "umap", group.by = "tissue", pt.size = 0.6, cols = pal)


DefaultAssay(ecs_integrated) <- "RNA"
DimPlot(object = ecs_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
FeaturePlot(ecs_integrated, c("Nkx6-1"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 0.8) # Figure S5C

c7_vs_c8 <- FindMarkers(ecs_integrated, ident.1 = 7,
                        ident.2 = 8, only.pos = TRUE,
                        logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
c8_vs_c7 <- FindMarkers(ecs_integrated, ident.1 = 8,
                        ident.2 = 7, only.pos = TRUE,
                        logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1) # merge
c9_vs_c2 <- FindMarkers(ecs_integrated, ident.1 = 9,
                        ident.2 = 2, only.pos = TRUE,
                        logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
c1_vs_c0 <- FindMarkers(ecs_integrated, ident.1 = 1,
                        ident.2 = 0, only.pos = TRUE,
                        logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
c0_vs_c1 <- FindMarkers(ecs_integrated, ident.1 = 0,
                        ident.2 = 1, only.pos = TRUE,
                        logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)

c12_vs_c7 <- FindMarkers(ecs_integrated, ident.1 = 12,
                         ident.2 = 7, only.pos = TRUE,
                         logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
c7_vs_c12 <- FindMarkers(ecs_integrated, ident.1 = 7,
                         ident.2 = 12, only.pos = TRUE,
                         logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
c1_vs_c0 <- FindMarkers(ecs_integrated, ident.1 = 1,
                        ident.2 = 0, only.pos = TRUE,
                        logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)


# Merge clusters that do not seem biologically relevant
Idents(ecs_integrated) <- ecs_integrated$integrated_snn_res.0.8
DimPlot(object = ecs_integrated, reduction = "umap", label = T, repel = TRUE, pt.size = 0.8, cols = mypal)

ecs_integrated <- RenameIdents(ecs_integrated, `0` = "0", `1` = "1", `2` = "2",
                               `3` = "3", `4` = "0", `5` = "4",
                               `6` = "3", `7` = "5", `8` = "2", `9` = "2",
                               `10` = "6", `11` = "7", `12` = "8", `13` = "9", `14` = "10")

DimPlot(object = ecs_integrated, reduction = "umap", label = F, repel = TRUE, pt.size = 0.8)
DimPlot(object = ecs_integrated, reduction = "umap", label = F, repel = TRUE, pt.size = 0.6, split.by = "tissue", ncol = 3)

ecs_integrated$merged_clusters <- Idents(ecs_integrated)

table(ecs_integrated$merged_clusters)

names(ecs_integrated[[]])

Idents(ecs_integrated) <- ecs_integrated$tissue
tissuepal <- rev(brewer.pal(n=9, name = 'Set3'))

DimPlot(object = ecs_integrated, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 0.8, cols = tissuepal)
DimPlot(object = ecs_integrated, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 0.8, group.by = "dataset")

table(Idents(ecs_integrated))

# Differential expression analysis
Idents(ecs_integrated) <- ecs_integrated$merged_clusters
DimPlot(object = ecs_integrated, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 0.8, cols = mypal)
DimPlot(object = ecs_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.6, cols = mypal, split.by = "tissue", ncol = 3)


# Find markers
DefaultAssay(ecs_integrated) <- "RNA"
markers_ecs_integrated <- FindAllMarkers(ecs_integrated, only.pos = TRUE, test.use = "wilcox",
                                         logfc.threshold = 0.25, min.pct = 0.25, pseudocount.use = 1)

c5_vs_c1 <- FindMarkers(ecs_integrated, ident.1 = 5,
                        ident.2 = 1, only.pos = TRUE,
                        logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
c1_vs_c5 <- FindMarkers(ecs_integrated, ident.1 = 1,
                        ident.2 = 5, only.pos = TRUE,
                        logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)

c5_vs_c2 <- FindMarkers(ecs_integrated, ident.1 = 5,
                        ident.2 = 2, only.pos = TRUE,
                        logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)

# Assign cell type
Idents(ecs_integrated) <- ecs_integrated$merged_clusters
ecs_integrated <- RenameIdents(ecs_integrated,
                               `0` = "Brain ependymal cells (likely E1 cells)",
                               `1` = "Pons, medulla, and spinal cord cells",
                               `2` = "Lateral spinal cord cells",
                               `3` = "Lateral spinal cord cells (activated)",
                               `4` = "Ventral spinal cord cells",
                               `5` = "Lateral spinal cord cells (mature)",
                               `6` = "Dorsal thalamus, midbrain, and spinal cord cells",
                               `7` = "Dorsal pons and medulla cells",
                               `8` = "Unknown brain cells",
                               `9` = "Ventral pons and medulla cells",
                               `10` = "Proliferating cells")

ecs_integrated$cell_type_integrated <- Idents(ecs_integrated)

table(Idents(ecs_integrated))


# Find markers
markers_ecs_integrated <- FindAllMarkers(ecs_integrated, only.pos = TRUE, test.use = "wilcox",
                                         logfc.threshold = 0.25, min.pct = 0.25, pseudocount.use = 1)


# Plotting with ggplot2
extract_plot_data <- function(single_cell_object, gene = NULL) {
  
  temp <- FeaturePlot(single_cell_object, features = c(gene), reduction = "umap") 
  
  plot_df <- tibble(
    x = temp$data$UMAP_1,
    y = temp$data$UMAP_2,
    value = c(temp$data[,gene])
  )
  
  return(plot_df)
}

plot_df <- tibble(
  gene_name = c("Cpe") 
) %>% mutate(map(gene_name, ~extract_plot_data(ecs_integrated, .x))) %>%
  unnest() %>% print()

plot_title <- plot_df$gene_name %>% unique()

ggplot(plot_df %>% filter(value == 0), aes(x, y)) +
  geom_point(color = '#d3d3d3', size = 0.8) +
  geom_point(data = plot_df %>% filter(value != 0), aes(x, y, color = value), size = 0.8) +
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



# Plots showing from which tissues cells in each cluster come from

# How many cells in each cluster?
table(Idents(ecs_integrated))

# From which tissue?
table(Idents(ecs_integrated), ecs_integrated@meta.data$tissue)
tissue_counts <- table(Idents(ecs_integrated), ecs_integrated@meta.data$tissue)

head(tissue_counts)

# From which dataset?
table(Idents(ecs_integrated), ecs_integrated@meta.data$dataset)
dataset_counts <- table(Idents(ecs_integrated), ecs_integrated@meta.data$dataset)

# Stacked barplot showing cells/tissue for each cluster
tissue_counts <- as.data.frame(tissue_counts)
head(tissue_counts)

tissue_counts <- tissue_counts %>% 
  rename(cluster = Var1, tissue = Var2, count = Freq)

head(tissue_counts)


# Percent stacked barplot
ggplot(tissue_counts, aes(fill = tissue, y = count, x = cluster)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_brewer(palette = "Set3", direction = -1) +
  theme_bw() + RotatedAxis()


Idents(ecs_integrated) <- ecs_integrated$merged_clusters
saveRDS(ecs_integrated, file = "~/path/to/zeisel_rodrigo_ecs_integrated_featuresmatched.rds")

# Figure 4A
Idents(ecs_integrated) <- ecs_integrated$merged_clusters

DimPlot(object = ecs_integrated, reduction = "umap", label = F, repel = TRUE, pt.size = 0.6,
        cols = c(`0` = "#1db751", `1` = "#f081a6",  #"#fb6568", 
                 `2` = "#cb181d", `3` = "#fb3235", `4` = "#42a4f2",
                 `5` = "#fb65c5", `6` = "#fbaf2a",
                 `7` = "#fedf32", `8` = "lightgrey",
                 `9` = "#2a62ff", `10` = "#a94bb9"))


