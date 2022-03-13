
## Analysis of single-cell transcriptomes from cells of the spinal cord central canal region from young mice

## Loading libraries
library(Seurat)
library(tidyverse)
library(cowplot)

## Loading the count matrices (CellRanger's output) and creating a Seurat object

# Load the data from sample 10X01_1
young1.data <- Read10X(data.dir = "~/path/to/cellranger_count/outs/10X01_1/filtered_feature_bc_matrix/")
young1.data[10:15, 10:15]

# Create a Seurat object
young1 <- CreateSeuratObject(counts = young1.data, project = "10X01_1") # project is sample

# Load doublet info and add to metadata slot
doublet_info <- read_csv("~/path/to/cellranger_count/outs/10X01_1/scrublet_output_table.csv") %>%
  as.data.frame()

cell_names_1 <- WhichCells(young1)
rownames(doublet_info) <- cell_names_1

head(doublet_info)

young1 <- AddMetaData(young1, metadata = doublet_info, col.name = c("doublet_score", "doublet"))

# View metadata (stored in object@meta.data)
young1[[]]
names(young1[[]])

table(young1[["orig.ident"]]) # OK!
table(young1[["doublet"]]) # OK!

# Load the data from sample 10X03_1
young2.data <- Read10X(data.dir = "~/path/to/cellranger_count/outs/10X03_1/filtered_feature_bc_matrix/")

young2 <- CreateSeuratObject(counts = young2.data, project = "10X03_1")

# Load doublet info and add to metadata slot
doublet_info <- read_csv("~/path/to/cellranger_count/outs/10X03_1/scrublet_output_table.csv") %>%
  as.data.frame()

cell_names_2 <- WhichCells(young2)
rownames(doublet_info) <- cell_names_2

head(doublet_info)

young2 <- AddMetaData(young2, metadata = doublet_info, col.name = c("doublet_score", "doublet"))

table(young2[["doublet"]]) # OK!

# Load the data from sample 10X05_2
young3.data <- Read10X(data.dir = "~/path/to/cellranger_count/outs/10X05_2/filtered_feature_bc_matrix/")

young3 <- CreateSeuratObject(counts = young3.data, project = "10X05_2")

# Load doublet info and add to metadata slot
doublet_info <- read_csv("~/path/to/cellranger_count/outs/10X05_2/scrublet_output_table.csv") %>%
  as.data.frame()

cell_names_3 <- WhichCells(young3)
rownames(doublet_info) <- cell_names_3
head(doublet_info)

young3 <- AddMetaData(young3, metadata = doublet_info, col.name = c("doublet_score", "doublet"))
table(young3[["doublet"]]) # OK!


# Load the data from sample 10X06_2
young4.data <- Read10X(data.dir = "~/path/to/cellranger_count/outs/10X06_2/filtered_feature_bc_matrix/")

young4 <- CreateSeuratObject(counts = young4.data, project = "10X06_2")

# Load doublet info and add to metadata slot
doublet_info <- read_csv("~/path/to/cellranger_count/outs/10X06_2/scrublet_output_table.csv") %>%
  as.data.frame()

cell_names_4 <- WhichCells(young4)
rownames(doublet_info) <- cell_names_4
head(doublet_info)

young4 <- AddMetaData(young4, metadata = doublet_info, col.name = c("doublet_score", "doublet"))
table(young4[["doublet"]]) # OK!


# Merge all objects from young mice
young_merged <- merge(x = young1, y = c(young2, young3, young4))

young_merged

# Add the age column to metadata slot
age <- rep("young", 12279) # 12279 cells
young_merged <- AddMetaData(young_merged, metadata = age, col.name = "age")

names(young_merged[[]])
table(young_merged$doublet)


# Pull the raw expression matrix and metadata to create a new Seurat object with genes expressed in fewer than 3  cells filtered out
raw_data <- as.matrix(GetAssayData(young_merged, slot = "counts"))
metadata <- young_merged@meta.data
head(metadata)

young <- CreateSeuratObject(counts = raw_data, meta.data = metadata, min.cells = 3, min.features = 500)
young


## Data QC

# Quick visualisation
VlnPlot(young, c("nFeature_RNA", "nCount_RNA"), pt.size = 0)

FeatureScatter(young, feature1 = "nFeature_RNA", feature2 = 'nCount_RNA')

young[["pct_counts_mito"]] <- PercentageFeatureSet(young, pattern = "^mt-")
young[["pct_counts_ribo"]] <- PercentageFeatureSet(young, pattern = "^Rp[sl][[:digit:]]")

VlnPlot(young, features = c("pct_counts_mito", "pct_counts_ribo"), pt.size = 0, ncol = 2)
RidgePlot(young, features = c("nFeature_RNA", "nCount_RNA", "pct_counts_mito", "pct_counts_ribo"), ncol = 2)

FeatureScatter(young, feature1 = "nFeature_RNA", feature2 = "pct_counts_mito", group.by = "orig.ident")
FeatureScatter(young, feature1 = "nCount_RNA", feature2 = "pct_counts_mito", group.by = "orig.ident")


# Which cells have high percentage of mitochondrial genes?
median(young$pct_counts_mito)
mad(young$pct_counts_mito)
(median(young$pct_counts_mito) + 3 * mad(young$pct_counts_mito)) # 20.30244

high_mito_cells <- WhichCells(object = young, expression = pct_counts_mito > 20.30244)

length(x = high_mito_cells)  # 1562 cells

# Remove cells with high percent of mito genes
young <- subset(young, subset = pct_counts_mito < 20.30244)
young # 9781 cells


# Fetch QC data from the object's metadata
names(young[[]])

qc_data <- FetchData(object = young, vars = c("orig.ident", "nCount_RNA", "nFeature_RNA", "pct_counts_mito", "pct_counts_ribo"))
colnames(qc_data)

# Boxplot of QC metrics by mouse
reads <- ggplot(data = qc_data, aes(x = orig.ident, y = nCount_RNA)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_linedraw() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

features <- ggplot(data = qc_data, aes(x = orig.ident, y = nFeature_RNA)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_linedraw() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

mito <- ggplot(data = qc_data, aes(x = orig.ident, y = pct_counts_mito)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_linedraw() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

ribo <- ggplot(data = qc_data, aes(x = orig.ident, y = pct_counts_ribo)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_linedraw() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

reads + features+ mito + ribo



## Processing

# Normalisation
young <- NormalizeData(young, normalization.method = "LogNormalize", scale.factor = 10000)

# Finding variable features
young <- FindVariableFeatures(young, selection.method = "vst", nfeatures = 2500)

# Scaling the data
all.genes <- rownames(young)
young <- ScaleData(young, features = all.genes)

# PCA
young <- RunPCA(young, features = VariableFeatures(object = young))

ElbowPlot(young)

DimHeatmap(young, dims = 12:20, cells = 500, balanced = TRUE)
DimHeatmap(young, dims = 21:26, cells = 500, balanced = TRUE) # 21 PCs

# UMAP
young <- RunUMAP(object = young, dims = 1:21, seed.use = 33,
                    learning.rate = 1, n.neighbors = 30, min.dist = 0.3)
DimPlot(young, reduction = "umap", pt.size = 0.8, label = TRUE)


# Clustering
young <- FindNeighbors(object = young, dims = 1:21, k.param = 10)

young <- FindClusters(object = young, graph.name = NULL,
                         modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                         node.sizes = NULL, resolution = seq(0.1, 0.8, 0.1), algorithm = 1, n.start = 1000,
                         n.iter = 10, random.seed = 0, temp.file.location = NULL,
                         edge.file.name = NULL, verbose = TRUE)

young <- FindClusters(object = young, graph.name = NULL,
                         modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                         node.sizes = NULL, resolution = 1, algorithm = 1, n.start = 1000,
                         n.iter = 10, random.seed = 0, temp.file.location = NULL,
                         edge.file.name = NULL, verbose = TRUE)

young <- FindClusters(object = young, graph.name = NULL,
                         modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                         node.sizes = NULL, resolution = 1.2, algorithm = 1, n.start = 1000,
                         n.iter = 10, random.seed = 0, temp.file.location = NULL,
                         edge.file.name = NULL, verbose = TRUE)

# Clustering trees
library(clustree)

names(young[[]])

clustree(young, prefix = "RNA_snn_res.")


# Set cell identities
Idents(object = young) <- young[["RNA_snn_res.1"]]
young$seurat_clusters <- young[["RNA_snn_res.1"]]
my_cols <- DiscretePalette(n = 37, palette = "polychrome") # options: "alphabet", "alphabet2", "glasbey", "polychrome", and "stepped"
DimPlot(object = young, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.8, cols = my_cols)
DimPlot(object = young, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.8)

DimPlot(young, reduction = "umap", label = FALSE, pt.size = 0.6, split.by = "orig.ident", cols = my_cols, ncol = 2)
DimPlot(young, reduction = "umap", label = FALSE, pt.size = 0.6, split.by = "orig.ident", ncol = 2)

# Find markers
Idents(object = young) <- young[["RNA_snn_res.1"]]
res1 <- FindAllMarkers(object = young, genes.use = NULL, only.pos = TRUE,
                       logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.20, pseudocount.use = 1)
res1_top10 <- res1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)f

# Cells in cluster 25 are suspected CSF-contacting neurons and ependymal cell doublets > remove
res1_c32_vs_c25 <- FindMarkers(young, ident.1 = 32, ident.2 = 25, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
res1_c25_vs_c32 <- FindMarkers(young, ident.1 = 25, ident.2 = 32, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)

# Cells in cluster 26 are suspected astrocyte and microglia doublet > remove
res1_c26_vs_c6 <- FindMarkers(young, ident.1 = 26, ident.2 = 6, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)

# Cells in cluster 18 are a mix of neuronal and microglial genes (with none expressed broadly across the cluster) > remove

# Cells in cluster 11 are of low quality (and most of them derive from only one mouse) > remove
res1_c11_vs_c7 <- FindMarkers(young_filtered, ident.1 = 11, ident.2 = 7, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
res1_c7_vs_c11 <- FindMarkers(young_filtered, ident.1 = 7, ident.2 = 11, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)

VlnPlot(young_filtered, c("nCount_RNA", "nFeature_RNA", "pct_counts_mito"), pt.size = 0, cols = my_cols, ncol = 1, idents = c("7","11"))


FeaturePlot(young_filtered, c("Rgs5"), order = TRUE, cols = c("lightgrey", "red"), pt.size = 0.6, reduction = "umap")
VlnPlot(young_filtered, c("nCount_RNA", "nFeature_RNA", "pct_counts_mito"), pt.size = 0, cols = my_cols, ncol = 1)


## Filtering out doublets and low-quality cells
young_filtered <- subset(x= young, idents = c(0,1,2,3,4,5,6,7,8,9,10,12,13,14,15,16, 17,
                                                       19,20,21,22,23,24,27,28,29,30,31,32))

# Find highly variable features
young_filtered <- FindVariableFeatures(young_filtered, selection.method = "vst", nfeatures = 2500)

# Scaling the data
all.genes <- rownames(young_filtered)
young_filtered <- ScaleData(young_filtered, features = all.genes)

# PCA
young_filtered <- RunPCA(young_filtered, features = VariableFeatures(object = young_filtered))

# Find which PCs contribute to most of the variability
ElbowPlot(young_filtered)

DimHeatmap(young_filtered, dims = 13:21, cells = 2000, balanced = TRUE)

# UMAP
young_filtered <- RunUMAP(object = young_filtered, dims = 1:21, seed.use = 4,
                          learning.rate = 1, n.neighbors = 30, min.dist = 0.3)
DimPlot(young_filtered, reduction = "umap", pt.size = 0.8, cols = my_cols) + NoLegend()




# Clustering
young_filtered <- FindNeighbors(object = young_filtered, dims= 1:21, k.param = 10)

young_filtered <- FindClusters(object = young_filtered, graph.name = NULL,
                               modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                               node.sizes = NULL, resolution = seq(0.1, 1.2, 0.1), algorithm = 1, n.start = 1000,
                               n.iter = 10, random.seed = 0, temp.file.location = NULL,
                               edge.file.name = NULL, verbose = TRUE)


# Clustering tree
clustree(young_filtered, prefix = "RNA_snn_res.")

# Set cell identities
Idents(object = young_filtered) <- young_filtered[["RNA_snn_res.1"]]
my_cols <- DiscretePalette(n = 31, palette = "polychrome")

DimPlot(object = young_filtered, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.6, cols = my_cols)

FeaturePlot(young_filtered, c("Sntn"), cols = c("lightgrey", "red"), order = TRUE,
            pt.size = 0.8, ncol = 1)


# How does clustering vary by mouse?
table(Idents(young_filtered), young_filtered@meta.data$orig.ident)
table(Idents(young_filtered))

prop.table(x = table(Idents(young_filtered), young_filtered@meta.data$orig.ident), margin = 2)


# Clusters QC
VlnPlot(young_filtered, features = c("nFeature_RNA", "nCount_RNA"), pt.size = 0, ncol = 1, cols = my_cols)

# Find markers
Idents(object = young_filtered) <- young_filtered[["RNA_snn_res.1"]]
res1 <- FindAllMarkers(object = young_filtered, genes.use = NULL, only.pos = TRUE,
                       logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.20, pseudocount.use = 1)
res1_top10 <- res1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

## Removing clusters with mixed signatures
young_filtered <- subset(x= young_filtered, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
                                                       16,17,18,19,20,21,23,24,25,27,28,29))

# Find highly variable features
young_filtered <- FindVariableFeatures(young_filtered, selection.method = "vst", nfeatures = 2500)

# Scaling the data
all.genes <- rownames(young_filtered)
young_filtered <- ScaleData(young_filtered, features = all.genes)

# PCA
young_filtered <- RunPCA(young_filtered, features = VariableFeatures(object = young_filtered))

# Find which PCs contribute to most of the variability
ElbowPlot(young_filtered)

DimHeatmap(young_filtered, dims = 13:21, cells = 2000, balanced = TRUE)

# UMAP
young_filtered <- RunUMAP(object = young_filtered, dims = 1:21, seed.use = 0,
                          learning.rate = 1, n.neighbors = 30, min.dist = 0.3)
DimPlot(young_filtered, reduction = "umap", pt.size = 0.8, cols = my_cols) + NoLegend()


# Clustering
young_filtered <- FindNeighbors(object = young_filtered, dims= 1:21, k.param = 10)

young_filtered <- FindClusters(object = young_filtered, graph.name = NULL,
                               modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                               node.sizes = NULL, resolution = seq(0.1, 1.2, 0.1), algorithm = 1, n.start = 1000,
                               n.iter = 10, random.seed = 0, temp.file.location = NULL,
                               edge.file.name = NULL, verbose = TRUE)

young_filtered <- FindClusters(object = young_filtered, graph.name = NULL,
                               modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                               node.sizes = NULL, resolution = 1.4, algorithm = 1, n.start = 1000,
                               n.iter = 10, random.seed = 0, temp.file.location = NULL,
                               edge.file.name = NULL, verbose = TRUE)


# Clustering tree
clustree(young_filtered, prefix = "RNA_snn_res.")

# Set cell identities
Idents(object = young_filtered) <- young_filtered[["RNA_snn_res.1"]]
my_cols <- DiscretePalette(n = 31, palette = "polychrome")

DimPlot(object = young_filtered, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.1, cols = my_cols) # Figure S1B
DimPlot(object = young_filtered, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 0.1, cols = my_cols, 
        split.by = "orig.ident", ncol = 2) + NoLegend() # Figure S1C

FeaturePlot(young_filtered, c("Vtn"), cols = c("lightgrey", "red"), order = TRUE,
            pt.size = 0.8, ncol = 1)

DimPlot(young_filtered, reduction = "umap", label = FALSE, pt.size = 1, split.by = "orig.ident", cols = my_cols, ncol = 2)

# How does clustering vary by mouse?
table(Idents(young_filtered), young_filtered@meta.data$orig.ident)
table(Idents(young_filtered))

prop.table(x = table(Idents(young_filtered), young_filtered@meta.data$orig.ident), margin = 2)

# Clusters QC
VlnPlot(young_filtered, features = c("nFeature_RNA", "nCount_RNA"), pt.size = 0, ncol = 1, cols = my_cols)

# Find markers
Idents(object = young_filtered) <- young_filtered[["RNA_snn_res.1"]]
res1 <- FindAllMarkers(object = young_filtered, genes.use = NULL, only.pos = TRUE,
                       logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.20, pseudocount.use = 1)
res1_top10 <- res1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

# Assign cell types
Idents(object = young_filtered) <- young_filtered[["RNA_snn_res.1"]]
young_filtered <- RenameIdents(young_filtered, `0` = "Microglia", `1` = "Microglia", `2` = "Microglia",
                               `3` = "Ependymal cells", `4` = "Ependymal cells", `5` = "Astrocytes", `6` = "Microglia",
                               `7` = "Oligodendrocytes", `8` = "Ependymal cells", `9` = "Microglia", `10` = "Microglia",
                               `11` = "Ependymal cells", `12` = "Ependymal cells", `13` = "Pericytes", `14` = "Vascular endothelial cells",
                               `15` = "Microglia", `16` = "Macrophages", `17` = "Vascular leptomeningeal cells", `18` = "Vascular leptomeningeal cells",
                               `19` = "T cells", `20` = "Microglia", `21` = "Microglia",
                               `22` = "Ependymal cells", `23` = "Vascular leptomeningeal cells", `24` = "Schwann cells",
                               `25` = "Vascular endothelial cells", `26` = "CSF-contacting neurons")

young_filtered$cell_type <- Idents(young_filtered)
table(young_filtered$cell_type)

young_filtered$cell_type <- factor(x = young_filtered$cell_type, levels = c("Astrocytes", "Ependymal cells", "Oligodendrocytes", "Schwann cells",
                                                                            "Microglia", "Macrophages", "Pericytes", "Vascular leptomeningeal cells", "Vascular endothelial cells",
                                                                            "T cells", "CSF-contacting neurons"))
young_palette <- c("#cb181d", "#fb3235", "#fb6568", "#fb65c5", "#1480d0", "#20b6c7", "#fe6e00", "#fbaf2a", "#fedf32", "#8b5e3c", "#a94bb9")


Idents(young_filtered) <- young_filtered$cell_type

DimPlot(young_filtered, pt.size = 1) + NoLegend()

# Figure 1C
DimPlot(object = young_filtered, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 0.1,
        cols = c( "Oligodendrocytes" = "#fb6568", "Ependymal cells" = "#fb3235", "Astrocytes" = "#cb181d",
                  "T cells" = "#8b5e3c",
                  "Pericytes" = "#fe6e00", "Vascular endothelial cells" = "#fedf32", "Vascular leptomeningeal cells" = "#fbaf2a",
                  "Microglia" = "#1480d0", "Macrophages" = "#20b6c7",
                  "CSF-contacting neurons" = "#a94bb9",
                  "Schwann cells" = "#fb65c5"))


DimPlot(young_filtered, reduction = "umap", label = FALSE, pt.size = 0.1, cols = young_palette) # Figure 1


# Find markers
young_celltype_markers <- FindAllMarkers(object = young_filtered, genes.use = NULL, only.pos = TRUE,
                                         logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)

# Ependymal cell markers
ependymal_markers <- FindMarkers(young_filtered, ident.1 = "Ependymal cells", only.pos = TRUE,
                                 logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)

# Figure 1D
selected_cell_type_markers <- c("Atp1b2","Aldoc","Ecrg4","Nnat", "Mia", "Cspg5","Olig1","Mpz",
                                "Sostdc1","Ctss","Hexb","Mrc1" ,"Pf4", "Vtn","Rgs5","Apod","Dcn","Cldn5","Itm2a",
                                "Nkg7","Trbc1", "Pkd2l1", "Resp18")

# Dot plot of top 2 markers
DotPlot(young_filtered, features = rev(selected_cell_type_markers), dot.scale = 6, col.min = 0,
        idents = c("Astrocytes", "Ependymal cells", "Oligodendrocytes", "Schwann cells", "Microglia", "Macrophages", "Pericytes",
                   "Vascular leptomeningeal cells", "Vascular endothelial cells", "T cells", "CSF-contacting neurons")) + 
  RotatedAxis() +
  scale_colour_viridis_c(option = "viridis") + 
  coord_flip()

table(Idents(young_filtered))


# Saving the object
Idents(young_filtered) <- young_filtered$cell_type
saveRDS(young_filtered, file = "~/Documents/Single-cell_RNA-seq/10x/data analysis/objects/young_filtered_cell_types.rds")
young_filtered <- readRDS(file = "~/Documents/Single-cell_RNA-seq/10x/data analysis/objects/young_filtered_cell_types.rds")

DimPlot(object = young_filtered, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 0.1,
        cols = c( "Oligodendrocytes" = "#fb6568", "Ependymal cells" = "#fb3235", "Astrocytes" = "#cb181d",
                  "T cells" = "#8b5e3c",
                  "Pericytes" = "#fe6e00", "Vascular endothelial cells" = "#fedf32", "Vascular leptomeningeal cells" = "#fbaf2a",
                  "Microglia" = "#1480d0", "Macrophages" = "#20b6c7",
                  "CSF-contacting neurons" = "#a94bb9",
                  "Schwann cells" = "#fb65c5"))



## QC plots of the final, filtered Seurat object

# Median number of UMIs per cell
median(young_filtered$nCount_RNA)

# Median
median(young_filtered$nFeature_RNA)


# Supplementary figures

# Fetch QC data from the object's metadata
names(young_filtered[[]])

young_qc_data <- FetchData(object = young_filtered, vars = c("orig.ident", "cell_type", "nCount_RNA", "nFeature_RNA", "pct_counts_mito"))

young_qc_data <- as.data.frame(young_qc_data)
head(young_qc_data)

# Boxplot of QC metrics by cluster ID
young_reads <- ggplot(data = young_qc_data, aes(x = cell_type, y = nCount_RNA)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

young_features <- ggplot(data = young_qc_data, aes(x = cell_type, y = nFeature_RNA)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

young_mito <- ggplot(data = young_qc_data, aes(x = cell_type, y = pct_counts_mito)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_linedraw()

young_reads + young_features


# Figure S1A

# Boxplot of QC metrics by sample/mouse
umis_by_mouse <- ggplot(data = young_qc_data, aes(x = orig.ident, y = nCount_RNA)) +
  geom_boxplot(outlier.size = 0.3, fill = "grey90") +
  #theme_linedraw() +
  theme_classic() +
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(axis.text.x = element_blank()) +
  labs(x='')

features_by_mouse <- ggplot(data = young_qc_data, aes(x = orig.ident, y = nFeature_RNA)) +
  geom_boxplot(outlier.size = 0.3, fill = "grey90") +
  #theme_linedraw() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x='')

umis_by_mouse + features_by_mouse

# Figure S1A
top5_young_umis <- ggplot(subset(young_qc_data, cell_type %in% c("Microglia", "Ependymal cells", "Astrocytes", "Oligodendrocytes", "Vascular leptomeningeal cells")),
                          aes(factor(x = cell_type, levels = c("Microglia", "Ependymal cells", "Astrocytes", "Oligodendrocytes", "Vascular leptomeningeal cells")),
                              y = nCount_RNA)) +
  geom_boxplot(outlier.size = 0.3, fill = "grey90") + # use outlier.shape = NA to not show outliers
  #theme_linedraw() +
  theme_classic() +
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(axis.text.x = element_blank()) +
  labs(x='', y='')

top5_young_features <- ggplot(subset(young_qc_data, cell_type %in% c("Microglia", "Ependymal cells", "Astrocytes", "Oligodendrocytes", "Vascular leptomeningeal cells")),
                              aes(factor(x = cell_type, levels = c("Microglia", "Ependymal cells", "Astrocytes", "Oligodendrocytes", "Vascular leptomeningeal cells")),
                                  y = nFeature_RNA)) +
  geom_boxplot(outlier.size = 0.3, fill = "grey90") +
  #theme_linedraw() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x='', y='')

top5_young_umis + top5_young_features

plot_grid(top5_young_umis, top5_young_features, ncol = 1, align = 'v', axis = 'lr')


# Plot together
plot_grid(umis_by_mouse, top5_young_umis,
          features_by_mouse, top5_young_features, ncol = 2, align = 'hv', axis = 'lr')



# Figure S1C

# How much does each mouse contribute to each cluster?
Idents(young_filtered) <- young_filtered$RNA_snn_res.1
mouse_to_clusters <- table(Idents(young_filtered), young_filtered@meta.data$orig.ident)

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
  theme_bw()




# Figure 1F

# Selected GO:BP terms enriched in ependymal cells for plotting
ependymal_goterms_selected <- read.csv(file = "~/Documents/Single-cell_RNA-seq/10x/data analysis/gprofiler/ependymal_cells_selected_GO-BP_FDR.csv",
                                       header = TRUE)

names(ependymal_goterms_selected)

ependymal_goterms_selected <- arrange(ependymal_goterms_selected, adjusted_p_value) # sort by decreasing average log fold change

names(ependymal_goterms_selected)

# Bar plot
ggplot(data = ependymal_goterms_selected, aes(x = reorder(stringr::str_to_sentence(term_name), negative_log10_of_adjusted_p_value),
                                              y = negative_log10_of_adjusted_p_value)) +
  geom_bar(stat = "identity", fill = "#56B4E9") +
  coord_flip() +
  theme_cowplot(12) +
  xlab(NULL)





## Subsetting ependymal cells
young <- readRDS(file = "~/Documents/Single-cell_RNA-seq/10X/data analysis/objects/young_lognorm.rds") # res = 1.0

Idents(young) <- young$RNA_snn_res.1

my_cols <- DiscretePalette(n = 33, palette = "polychrome")
DimPlot(object = young, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.8, cols = my_cols)

# Genes specific to cluster 27? (Astrocyte-like ependymal cells)
res1_c27_vs_c6 <- FindMarkers(young, ident.1 = 27, ident.2 = 6, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res1_c6_vs_c27 <- FindMarkers(young, ident.1 = 6, ident.2 = 27, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)

FeaturePlot(young, c("Nkx6-1"), cols = c("lightgrey", "red"), order = TRUE,
            pt.size = 1, ncol = 1)

young_ecs <- subset(x= young, idents = c(2,4,8,10,12,27))

young_ecs


# Find highly variable features
young_ecs <- FindVariableFeatures(young_ecs, selection.method = "vst", nfeatures = 2500)

# Scale
all.genes <- rownames(young_ecs)
young_ecs <- ScaleData(young_ecs)

# PCA
young_ecs <- RunPCA(young_ecs, features = VariableFeatures(object = young_ecs))

# Print PCs
print(x = young_ecs[["pca"]], dims = 1, nfeatures = 10)

# Visualise PCs
VizDimLoadings(object = young_ecs, dims = 1:2, reduction = "pca")

# Plot PCs
DimPlot(object = young_ecs, dims = c(1, 2), reduction = "pca")
DimPlot(object = young_ecs, dims = c(2, 3), reduction = "pca")


# Find which PCs contribute to most of the variability
ElbowPlot(young_ecs)

DimHeatmap(young_ecs, dims = 12:17, cells = 2000, balanced = TRUE) # 15 PCs

young_ecs <- RunUMAP(object = young_ecs, dims = 1:15, seed.use = 33,
                     learning.rate = 1, n.neighbors = 10, min.dist = 0.3)
DimPlot(young_ecs, reduction = "umap", pt.size = 0.8)
DimPlot(young_ecs, reduction = "umap", split.by = "orig.ident", pt.size = 0.8, ncol = 2)


## Clustering
young_ecs <- FindNeighbors(object = young_ecs, dims= 1:15, k.param = 10)

young_ecs <- FindClusters(object = young_ecs, graph.name = NULL,
                          modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                          node.sizes = NULL, resolution = seq(0.1, 1.2, 0.1), algorithm = 1, n.start = 1000,
                          n.iter = 10, random.seed = 0, temp.file.location = NULL,
                          edge.file.name = NULL, verbose = TRUE)

# Clustering tree
clustree(young_ecs, prefix = "RNA_snn_res.")


# Set cell identities
Idents(object = young_ecs) <- young_ecs[["RNA_snn_res.0.9"]]
my_cols_ecs <- DiscretePalette(n = 16, palette = "polychrome")
DimPlot(object = young_ecs, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1, cols = my_cols_ecs)
DimPlot(object = young_ecs, reduction = "tsne", label = TRUE, repel = TRUE, pt.size = 1, cols = my_cols_ecs) + NoLegend()

# Find markers
res09_ecs <- FindAllMarkers(object = young_ecs, genes.use = NULL, only.pos = TRUE,
                            logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)

# Removing perivascular macrophages, microglia, vascular endothelial cells, and astrocyte-like cells (likely doublets)
young_ecs <- subset(x= young_ecs, idents = c(0,1,2,3,4,5,6,7,8,14))

# Find highly variable features
young_ecs <- FindVariableFeatures(young_ecs, selection.method = "vst", nfeatures = 2500)

# Scale
all.genes <- rownames(young_ecs)
young_ecs <- ScaleData(young_ecs, features = all.genes)

# PCA
young_ecs <- RunPCA(young_ecs, features = VariableFeatures(object = young_ecs))

# Find which PCs contribute to most of the variability
ElbowPlot(young_ecs)

DimHeatmap(young_ecs, dims = 7:14, balanced = TRUE) # 11 PCs

young_ecs <- RunUMAP(object = young_ecs, dims = 1:11, seed.use = 4,
                     learning.rate = 95, n.neighbors = 10, min.dist = 0.3)
DimPlot(young_ecs, reduction = "umap", pt.size = 1, cols = my_cols_ecs)
DimPlot(ecs_clean, reduction = "umap", split.by = "orig.ident", pt.size = 1, cols = my_cols_ecs)


## Clustering
young_ecs <- FindNeighbors(object = young_ecs, dims= 1:11, k.param = 10)

young_ecs <- FindClusters(object = young_ecs, graph.name = NULL,
                          modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                          node.sizes = NULL, resolution = seq(0.1, 1.2, 0.1), algorithm = 1, n.start = 1000,
                          n.iter = 10, random.seed = 0, temp.file.location = NULL,
                          edge.file.name = NULL, verbose = TRUE)

# Clustering tree
clustree(young_ecs, prefix = "RNA_snn_res.")

# Set cell identities
Idents(object = young_ecs) <- young_ecs[["RNA_snn_res.0.4"]]
my_cols_ecs <- DiscretePalette(n = 10, palette = "polychrome")
DimPlot(object = young_ecs, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1, cols = my_cols_ecs)
DimPlot(object = young_ecs, reduction = "tsne", label = FALSE, repel = TRUE, pt.size = 1, cols = my_cols_ecs)

FeaturePlot(young_ecs, c("Bhmt"), cols = c("lightgrey", "red"), order = TRUE, pt.size = 1)


# How does clustering vary by mouse?
table(Idents(young_ecs), young_ecs@meta.data$orig.ident)
table(Idents(young_ecs))

prop.table(x = table(Idents(young_ecs), young_ecs@meta.data$orig.ident), margin = 2)


# Build phylogenetic tree
young_ecs <- BuildClusterTree(young_ecs, reorder = TRUE)

Tool(object = young_ecs, slot = 'BuildClusterTree')

PlotClusterTree(young_ecs)



# Find markers
res04_young_ecs <- FindAllMarkers(object = young_ecs, genes.use = NULL, only.pos = TRUE,
                                  logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)

res04_c1_vs_c0 <- FindMarkers(young_ecs, ident.1 = 1, ident.2 = 0, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res04_c0_vs_c1 <- FindMarkers(young_ecs, ident.1 = 0, ident.2 = 1, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res04_c4_vs_c2 <- FindMarkers(young_ecs, ident.1 = 4, ident.2 = 2, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res04_c2_vs_c4 <- FindMarkers(young_ecs, ident.1 = 2, ident.2 = 4, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res04_c6_vs_c4 <- FindMarkers(young_ecs, ident.1 = 6, ident.2 = 4, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res04_c4_vs_c6 <- FindMarkers(young_ecs, ident.1 = 4, ident.2 = 6, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res04_c6_vs_c1 <- FindMarkers(young_ecs, ident.1 = 6, ident.2 = 1, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res04_c1_vs_c6 <- FindMarkers(young_ecs, ident.1 = 1, ident.2 = 6, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res04_c0_vs_c1 <- FindMarkers(young_ecs, ident.1 = 0, ident.2 = 1, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res04_c1_vs_c0 <- FindMarkers(young_ecs, ident.1 = 1, ident.2 = 0, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res04_c3_vs_c4 <- FindMarkers(young_ecs, ident.1 = 3, ident.2 = 4, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res04_c4_vs_c3 <- FindMarkers(young_ecs, ident.1 = 4, ident.2 = 3, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)


# Merge cluster 0 and cluster 1 (Lateral cells, immature)
young_ecs <- RenameIdents(young_ecs, `0` = "0", `1` = "0", `2` = "1",
                          `3` = "2", `4` = "1", `5` = "3",
                          `6` = "0", `7` = "4")

DimPlot(young_ecs, reduction = "umap", label = FALSE, pt.size = 1, cols = my_cols_ecs)

young_ecs$merged_clusters <- Idents(young_ecs)

table(young_ecs$merged_clusters)

# Assign cell type
Idents(young_ecs) <- young_ecs$merged_clusters
young_ecs <- RenameIdents(young_ecs, `0` = "Lateral (immature)", `1` = "Lateral (activated)",
                          `2` = "Lateral (mature)",
                          `3` = "Ventral", `4` = "Dorsal")

DimPlot(young_ecs, reduction = "umap", label = FALSE, pt.size = 1)

young_ecs$cell_subtype <- Idents(young_ecs)

# Assign colours to cell types
DimPlot(object = young_ecs, reduction = "umap", label = FALSE, pt.size = 1,
        cols = c( "Lateral (immature)" = "#fb6568", "Lateral (activated)" = "#cb181d",
                  "Dorsal" = "#1db751",
                  "Lateral (mature)" = "#f4d941",
                  "Ventral" = "#20b6c7"))
table(Idents(young_ecs))

young_ecs_subtypes <- FindAllMarkers(object = young_ecs, genes.use = NULL, only.pos = TRUE,
                                     logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.20, pseudocount.use = 1)

saveRDS(young_ecs, file = "~/path/to/young_ecs_subtypes_without_astrolike_cells.rds")



