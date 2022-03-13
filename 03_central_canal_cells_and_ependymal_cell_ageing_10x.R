## Single-cell RNA-sequencing analysis of cells from the spinal cord central canal region from young and old mice (10x)

# Load libraries
library(Seurat)
library(tidyverse)
library(viridis)
options(future.globals.maxSize = 4000 * 1024^5)

## Loading and filtering the 10x dataset with cells from the the central canal region of aged mice

# Load the data from sample 10X01_1
old1.data <- Read10X(data.dir = "~/path/to/cellranger_count/outs/10X01_2/filtered_feature_bc_matrix/")
old1.data[10:15, 10:15]

# Create a Seurat object
old1 <- CreateSeuratObject(counts = old1.data, project = "10X01_2") # project is experiment

# Load doublet info and add to metadata slot
doublet_info <- read_csv("~/path/to/cellranger_count/outs/10X01_2/scrublet_output_table.csv") %>%
  as.data.frame()

cell_names_1 <- WhichCells(old1)
rownames(doublet_info) <- cell_names_1

head(doublet_info)

old1 <- AddMetaData(old1, metadata = doublet_info, col.name = c("doublet_score", "doublet"))

# View metadata (stored in object@meta.data)
old1[[]]
names(old1[[]])

table(old1[["orig.ident"]]) # OK!
table(old1[["doublet"]]) # OK!

# Load the data from sample 10X03_1
old2.data <- Read10X(data.dir = "~/path/to/cellranger_count/outs/10X03_2/filtered_feature_bc_matrix/")

old2 <- CreateSeuratObject(counts = old2.data, project = "10X03_2")

# Load doublet info and add to metadata slot
doublet_info <- read_csv("~/path/to/cellranger_count/outs/10X03_2/scrublet_output_table.csv") %>%
  as.data.frame()

cell_names_2 <- WhichCells(old2)
rownames(doublet_info) <- cell_names_2

head(doublet_info)

old2 <- AddMetaData(old2, metadata = doublet_info, col.name = c("doublet_score", "doublet"))

table(old2[["doublet"]]) # OK!

# Load the data from sample 10X05_2
old3.data <- Read10X(data.dir = "~/path/to/cellranger_count/outs/10X05_1/filtered_feature_bc_matrix/")

old3 <- CreateSeuratObject(counts = old3.data, project = "10X05_1") # project is experiment

# Load doublet info and add to metadata slot
doublet_info <- read_csv("~/path/to/cellranger_count/outs/10X05_1/scrublet_output_table.csv") %>%
  as.data.frame()

cell_names_3 <- WhichCells(old3)
rownames(doublet_info) <- cell_names_3
head(doublet_info)

old3 <- AddMetaData(old3, metadata = doublet_info, col.name = c("doublet_score", "doublet"))
table(old3[["doublet"]]) # OK!


# Load the data from sample 10X06_2
old4.data <- Read10X(data.dir = "~/path/to/cellranger_count/outs/10X06_1/filtered_feature_bc_matrix/")

old4 <- CreateSeuratObject(counts = old4.data, project = "10X06_1") # project is experiment



# Merge all objects from old mice
old_merged <- merge(x = old1, y = c(old2, old3, old4))

old_merged

# Add the age column to metadata slot
age <- rep("old", 12175) # 12175 cells
old_merged <- AddMetaData(old_merged, metadata = age, col.name = "age")

names(old_merged[[]])
table(old_merged$orig.ident)


# Pull the raw expression matrix and metadata to create a new Seurat object with genes expressed in fewer than 3  cells filtered out
raw_data <- as.matrix(GetAssayData(old_merged, slot = "counts"))
metadata <- old_merged@meta.data
head(metadata)

old <- CreateSeuratObject(counts = raw_data, meta.data = metadata, min.cells = 3, min.features = 500)
old

## Data QC

# Quick visualisation
VlnPlot(old, c("nFeature_RNA", "nCount_RNA"), pt.size = 0)

FeatureScatter(old, feature1 = "nFeature_RNA", feature2 = 'nCount_RNA')

old[["pct_counts_mito"]] <- PercentageFeatureSet(old, pattern = "^mt-")
old[["pct_counts_ribo"]] <- PercentageFeatureSet(old, pattern = "^Rp[sl][[:digit:]]")

VlnPlot(old, features = c("pct_counts_mito", "pct_counts_ribo"), pt.size = 0, ncol = 2)

FeatureScatter(old, feature1 = "nFeature_RNA", feature2 = "pct_counts_mito", group.by = "orig.ident")
FeatureScatter(old, feature1 = "nCount_RNA", feature2 = "pct_counts_mito", group.by = "orig.ident")


# Which cells have high percentage of mitochondrial genes?
median(old$pct_counts_mito)
mad(old$pct_counts_mito)
(median(old$pct_counts_mito) + 3 * mad(old$pct_counts_mito)) # 14.23161 (keeping all cells, percentage is lower than in cells from young mice)


# Fetch QC data from the object's metadata
names(old[[]])

qc_data <- FetchData(object = old, vars = c("orig.ident", "nCount_RNA", "nFeature_RNA", "pct_counts_mito", "pct_counts_ribo"))
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


saveRDS(old, file = "~/path/to/old-filtered.rds")






# Loading the filtered datasets 10x datasets of cells from the central canal region 
# of young (from script 01_central_canal_cells_young_10x) and aged mice (above)
young <- readRDS(file = "~/path/to/young_filtered_cell_types.rds")
old <-  readRDS(file = "~/path/to/old-filtered.rds")


# Merge datasets
allcells <- merge(x = young, y = old, merge.data = TRUE)

allcells # OK!

names(allcells[[]])

table(allcells[["age"]])

allcells$age <- factor(x = allcells$age, levels = c("young", "old"))

allcells <- NormalizeData(allcells, normalization.method = "LogNormalize", scale.factor = 10000)

allcells <- FindVariableFeatures(allcells, selection.method = "vst", nfeatures = 2500)

all.genes <- rownames(allcells)
allcells <- ScaleData(allcells, features = all.genes)

allcells <- RunPCA(allcells, features = VariableFeatures(object = allcells))

ElbowPlot(allcells)

DimHeatmap(allcells, dims = 24:32, cells = 500, balanced = TRUE) # 31 PCs

allcells <- RunUMAP(object = allcells, dims = 1:31, seed.use = 17,
                    learning.rate = 1, n.neighbors = 20, min.dist = 0.3)

DimPlot(allcells, reduction = "umap", pt.size = 0.6, label = F, repel = TRUE) + NoLegend()
DimPlot(allcells, reduction = "umap", split.by = "age", pt.size = 0.6)
DimPlot(allcells, reduction = "umap", split.by = "orig.ident", pt.size = 0.6, ncol = 4)



## Clustering
allcells <- FindNeighbors(object = allcells, dims = 1:31, k.param = 20)

allcells <- FindClusters(object = allcells, graph.name = NULL,
                         modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                         node.sizes = NULL, resolution = seq(0.1, 1, 0.1), algorithm = 1, n.start = 1000,
                         n.iter = 10, random.seed = 0, temp.file.location = NULL,
                         edge.file.name = NULL, verbose = TRUE)

allcells <- FindClusters(object = allcells, graph.name = NULL,
                         modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                         node.sizes = NULL, resolution = 1.2, algorithm = 1, n.start = 1000,
                         n.iter = 10, random.seed = 0, temp.file.location = NULL,
                         edge.file.name = NULL, verbose = TRUE)

# Clustering trees
library(clustree)

names(allcells[[]])

clustree(allcells, prefix = "RNA_snn_res.")


# Set cell identities
Idents(object = allcells) <- allcells[["RNA_snn_res.1"]]
my_cols <- DiscretePalette(n = 30, palette = "polychrome") # options: "alphabet", "alphabet2", "glasbey", "polychrome", and "stepped"
DimPlot(object = allcells, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.6, cols = my_cols)
DimPlot(object = allcells, reduction = "umap", split.by = "age", label = FALSE, repel = TRUE, pt.size = 0.6, cols = my_cols)

DimPlot(allcells, reduction = "umap", label = T, pt.size = 0.6, split.by = "orig.ident", cols = my_cols, ncol = 4)
DimPlot(allcells, reduction = "umap", label = FALSE, pt.size = 0.4, split.by = "age", cols = my_cols)
DimPlot(allcells, reduction = "umap", label = FALSE, pt.size = 0.4, split.by = "ident", cols = my_cols, ncol = 5)


# How does clustering vary by mouse?
table(Idents(allcells), allcells@meta.data$orig.ident)
table(Idents(allcells))

prop.table(x = table(Idents(allcells), allcells@meta.data$orig.ident), margin = 2)

# Clusters QC
VlnPlot(allcells, features = c("nFeature_RNA", "nCount_RNA"), pt.size = 0, ncol = 1)

DimPlot(allcells, reduction = "umap", label = TRUE, pt.size = 0.6, cols = my_cols)



# Find markers
res1 <- FindAllMarkers(object = allcells, genes.use = NULL, only.pos = TRUE,
                       logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res1_top10 <- res1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

VlnPlot(allcells, c("nCount_RNA", "nFeature_RNA", "pct_counts_mito"), pt.size = 0, cols = my_cols, ncol = 1)

Idents(object = allcells) <- allcells[["RNA_snn_res.1"]]

res1_c11_vs_c8 <- FindMarkers(allcells, ident.1 = 11, ident.2 = 8, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
res1_c8_vs_c11 <- FindMarkers(allcells, ident.1 = 8, ident.2 = 11, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
res1_c13_vs_c8 <- FindMarkers(allcells, ident.1 = 13, ident.2 = 8, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
res1_c8_vs_c13 <- FindMarkers(allcells, ident.1 = 8, ident.2 = 13, only.pos = TRUE,
                               logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
res1_c24_vs_c18 <- FindMarkers(allcells, ident.1 = 24, ident.2 = 18, only.pos = TRUE,
                               logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1) # mixed gene signature (remove)
res1_c18_vs_c24 <- FindMarkers(allcells, ident.1 = 18, ident.2 = 24, only.pos = TRUE,
                               logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
res1_c16_vs_c20 <- FindMarkers(allcells, ident.1 = 16, ident.2 = 20, only.pos = TRUE,
                               logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
res1_c20_vs_c16 <- FindMarkers(allcells, ident.1 = 20, ident.2 = 16, only.pos = TRUE,
                               logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)



# Assign cell type
Idents(object = allcells) <- allcells[["RNA_snn_res.1"]]

allcells <- RenameIdents(allcells, `0` = "Microglia", `1` = "Microglia", `2` = "Microglia",
                      `3` = "Microglia", `4` = "Ependymal cells", `5` = "Microglia", `6` = "Ependymal cells",
                      `7` = "Oligodendrocytes", `8` = "Astrocytes", `9` = "Microglia", `10` = "Ependymal cells",
                      `11` = "Astrocytes", `12` = "T cells", `13` = "Astrocytes", `14` = "Ependymal cells",
                      `15` = "Macrophages", `16` = "Vascular leptomeningeal cells", `17` = "Pericytes", `18` = "Vascular endothelial cells",
                      `19` = "Microglia", `20` = "Vascular leptomeningeal cells", `21` = "Ependymal cells",
                      `22` = "Proliferating cells", `23` = "Oligodendrocytes", `24` = "Vascular endothelial cells",
                      `25` = "B cells", `26` = "CSF-contacting neurons", `27` = "Schwann cells")

allcells$cell_type <- Idents(allcells)
table(allcells$cell_type)

allcells$cell_type <- factor(x = allcells$cell_type, levels = c("Astrocytes", "Ependymal cells", "Oligodendrocytes", "Schwann cells",
                                                                "Proliferating cells",
                                                                "Microglia",
                                                                "Pericytes", "Vascular leptomeningeal cells",
                                                                "Vascular endothelial cells",
                                                                "Macrophages",
                                                                "T cells", "B cells",
                                                                "CSF-contacting neurons"))

Idents(allcells) <- allcells$cell_type

DimPlot(allcells, pt.size = 1, label = T)

# Figure 3B
DimPlot(object = allcells, reduction = "umap", label = F, repel = TRUE, pt.size = 0.4,
        cols = c( "Oligodendrocytes" = "#fb6568", "Ependymal cells" = "#fb3235", "Astrocytes" = "#cb181d",
                  "Proliferating cells" = "#f081a6",
                  "Macrophages" = "#1db751", "T cells" = "#a6dba0", "B cells" = "#91c558",
                  "Pericytes" = "#fbaf2a", "Vascular endothelial cells" = "#fff7be", "Vascular leptomeningeal cells" = "#fedf32",
                  "Schwann cells" = "#fb65c5",
                  "Microglia" = "#1480d0", "Perivascular macrophages" = "#42a4f2",
                  "CSF-contacting neurons" = "#a94bb9"), split.by = "age")

saveRDS(allcells, file = "~/path/to/allcells_cell_types.rds")


# Find markers
allcells_celltype_markers <- FindAllMarkers(object = allcells, genes.use = NULL, only.pos = TRUE,
                       logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)



# Supplementary info

# Data QC plots for supplementary figures

# Fetch QC data from the object's metadata
names(allcells[[]])

allcells_qc <- FetchData(object = allcells, vars = c("cell_type", "orig.ident", "nCount_RNA", "nFeature_RNA", "age"))
colnames(allcells_qc)


allcells_qc <- as.data.frame(allcells_qc)
head(allcells_qc)

# Boxplot of QC metrics by cluster ID
allcells_reads <- ggplot(data = allcells_qc, aes(x = cell_type, y = nCount_RNA, fill = age)) +
  geom_boxplot(outlier.size = 0.2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # theme(axis.text.x=element_blank(),
  #       panel.grid = element_blank()) +
  labs(x='') + NoLegend()

allcells_features <- ggplot(data = allcells_qc, aes(x = cell_type, y = nFeature_RNA, fill = age)) +
  geom_boxplot(outlier.size = 0.2) +
  theme_classic() +
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(axis.text.x=element_blank(),
        panel.grid = element_blank()) +
  labs(x='')

#allcells_reads + allcells_features # Figure S5A

# Plot together
library(cowplot)
plot_grid(allcells_reads,
          allcells_features,
          ncol = 1, align = 'hv', axis = 'lr')



# Boxplot of QC metrics by sample/mouse
allcells_reads_by_mouse <- ggplot(data = allcells_qc, aes(x = orig.ident, y = nCount_RNA, fill = age)) +
  geom_boxplot(outlier.size = 0.2) +
  theme_classic() +
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(axis.text.x=element_blank(),
        panel.grid = element_blank()) +
  labs(x='') + NoLegend()

allcells_features_by_mouse <- ggplot(data = allcells_qc, aes(x = orig.ident, y = nFeature_RNA, fill = age)) +
  geom_boxplot(outlier.size = 0.2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x='') + NoLegend()

#allcells_reads_by_mouse + allcells_features_by_mouse

plot_grid(allcells_reads_by_mouse,
          allcells_features,
          allcells_features_by_mouse,
          allcells_reads,
          ncol = 2, align = 'hv', axis = 'lr')


# Boxplot of QC metrics by sample/mouse
allcells_qc <- FetchData(object = allcells, vars = c("orig.ident", "nCount_RNA", "nFeature_RNA", "age"))
colnames(allcells_qc)

df <- allcells_qc %>% 
  group_by(orig.ident) %>% 
  mutate(total_UMIs = sum(nCount_RNA), 
         total_features = sum(nFeature_RNA))

ggplot(data = df, aes(x = age, y = total_features, fill = age)) +
  geom_boxplot(outlier.size = 0.2) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        panel.grid = element_blank()) +
  labs(x='')

ggplot(data = df, aes(x = age, y = total_features, fill = age)) +
  geom_bar(position = 'dodge', stat = 'summary', fun.y = 'mean', width = 0.8) +
  geom_errorbar(stat = 'summary', position = position_dodge(0.8), width = 0.25) +
  #scale_fill_viridis_d(direction = -1) +
  geom_point(aes(x = age), shape = 21,
             position =  position_jitterdodge(jitter.width = 0.5, 
                                              dodge.width=0.8),
             size=1, alpha=0.8) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45,hjust=1)) # Figure 3D
  

allcells_reads_by_mouse <- ggplot(data = allcells_qc, aes(x = orig.ident, y = nCount_RNA, fill = age)) +
  geom_boxplot(outlier.size = 0.2) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        panel.grid = element_blank()) +
  labs(x='')

allcells_features_by_mouse <- ggplot(data = allcells_qc, aes(x = orig.ident, y = nFeature_RNA, fill = age)) +
  geom_boxplot(outlier.size = 0.2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x='')





# Median number of UMI/read counts
median(allcells$nCount_RNA)

# Median number of UMI/read counts
median(allcells$nFeature_RNA)


# UMAP plot split by mouse showing contribution of all samples to all clusters
Idents(object = allcells) <- allcells[["RNA_snn_res.1"]]
my_cols <- DiscretePalette(n = 30, palette = "polychrome")

DimPlot(object = allcells, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 0.3,
        split.by = "orig.ident", ncol = 4, cols = my_cols) # Figure S4C






## Subsetting ependymal cells
allcells <- readRDS(file = "~/path/to/allcells_cell_types.rds")

my_cols <- DiscretePalette(n = 30, palette = "polychrome") # options: "alphabet", "alphabet2", "glasbey", "polychrome", and "stepped"
DimPlot(object = allcells, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 0.6, cols = my_cols) + NoLegend()

# Subset clusters assigned to ependymal cells
ecs <- subset(x= allcells, idents = "Ependymal cells")

# Find highly variable features
ecs <- FindVariableFeatures(ecs, selection.method = "vst", nfeatures = 2500)

# Scale
all.genes <- rownames(ecs)
ecs <- ScaleData(ecs, features = all.genes)

# PCA
ecs <- RunPCA(ecs, features = VariableFeatures(object = ecs))

# Find which PCs contribute to most of the variability
ElbowPlot(ecs)

DimHeatmap(ecs, dims = 12:17, cells = 2000, balanced = TRUE)
DimHeatmap(ecs, dims = 17:22, cells = 2000, balanced = TRUE) # 16 PCs

ecs <- RunUMAP(object = ecs, dims = 1:16, seed.use = 1,
               learning.rate = 95, n.neighbors = 10, min.dist = 0.3)

DimPlot(ecs, reduction = "umap", pt.size = 0.6) + NoLegend()
DimPlot(ecs, reduction = "umap", group.by = "age", pt.size = 0.6)
DimPlot(ecs, reduction = "umap", split.by = "age", pt.size = 0.6)


## Clustering
ecs <- FindNeighbors(object = ecs, dims= 1:16, k.param = 10)

ecs <- FindClusters(object = ecs, graph.name = NULL,
                    modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                    node.sizes = NULL, resolution = seq(0.1, 1, 0.1), algorithm = 1, n.start = 1000,
                    n.iter = 10, random.seed = 0, temp.file.location = NULL,
                    edge.file.name = NULL, verbose = TRUE)

# Clustering tree
library(clustree)

names(ecs[[]])

clustree(ecs, prefix = "RNA_snn_res.")

# Set cell identities
Idents(object = ecs) <- ecs[["RNA_snn_res.0.6"]]
my_cols_ecs <- DiscretePalette(n = 15, palette = "polychrome")
DimPlot(object = ecs, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.6, cols = my_cols_ecs)
DimPlot(object = ecs, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.6, cols = my_cols_ecs, split.by = "age") + NoLegend()


# Find markers
res06_ecs <- FindAllMarkers(object = ecs, genes.use = NULL, only.pos = TRUE,
                            logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)

# Removing microglia and clusters with mixed gene signature
ecs <- subset(x= ecs, idents = c(0,1,2,3,4,5,6,8,10)) # "~/.../10x/data analysis/de_genes/ecs_ageing_clustering_01_res06.csv"

# Find highly variable features
ecs <- FindVariableFeatures(ecs, selection.method = "vst", nfeatures = 2500)

# Scale
all.genes <- rownames(ecs)
ecs <- ScaleData(ecs, features = all.genes)

# PCA
ecs <- RunPCA(ecs, features = VariableFeatures(object = ecs))

# Find which PCs contribute to most of the variability
ElbowPlot(ecs)

DimHeatmap(ecs, dims = 7:15, cells = 2000, balanced = TRUE) # 13 PCs

ecs <- RunUMAP(object = ecs, dims = 1:13, seed.use = 40, # 33 40 42 64
               learning.rate = 99, n.neighbors = 10, min.dist = 0.3)
DimPlot(ecs, reduction = "umap", pt.size = 1) + NoLegend()
DimPlot(ecs, reduction = "umap", split.by = "age", pt.size = 1)

## Clustering
ecs <- FindNeighbors(object = ecs, dims= 1:13, k.param = 10)

ecs <- FindClusters(object = ecs, graph.name = NULL,
                    modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                    node.sizes = NULL, resolution = seq(0.1, 1, 0.1), algorithm = 1, n.start = 1000,
                    n.iter = 10, random.seed = 0, temp.file.location = NULL,
                    edge.file.name = NULL, verbose = TRUE)

# Clustering tree
clustree(ecs, prefix = "RNA_snn_res.")

# Set cell identities
Idents(object = ecs) <- ecs[["RNA_snn_res.0.6"]]
my_cols_ecs <- DiscretePalette(n = 12, palette = "polychrome")
DimPlot(object = ecs, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.8, cols = my_cols_ecs)
DimPlot(object = ecs, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.8, cols = my_cols_ecs, split.by = "age")

FeaturePlot(ecs, c("Sntn"), cols = c("lightgrey", "red"), order = TRUE, pt.size = 0.8, split.by = "age")


# Find markers
res06_ecs <- FindAllMarkers(object = ecs, genes.use = NULL, only.pos = TRUE,
                            logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)

res04_c3_vs_c0 <- FindMarkers(ecs, ident.1 = 3, ident.2 = 0, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)


# Removing ependymal cells mixed with astrocyte cells and 3 choroid plexus cells
ecs <- subset(x= ecs, idents = c(0,1,2,3,4,5,6,7,9)) # "~/.../10x/data analysis/de_genes/ecs_ageing_clustering_02_res06.csv"

# Find highly variable genes
ecs <- FindVariableFeatures(ecs, selection.method = "vst", nfeatures = 2500)

# Scale
all.genes <- rownames(ecs)
ecs <- ScaleData(ecs, features = all.genes)

# PCA
ecs <- RunPCA(ecs, features = VariableFeatures(object = ecs))

# Find which PCs contribute to most of the variability
ElbowPlot(ecs)

DimHeatmap(ecs, dims = 7:15, cells = 2000, balanced = TRUE) # 13 PCs

ecs <- RunUMAP(object = ecs, dims = 1:13, seed.use = 14,
               learning.rate = 95, n.neighbors = 10, min.dist = 0.3)
DimPlot(ecs, reduction = "umap", pt.size = 1) + NoLegend()
DimPlot(ecs, reduction = "umap", split.by = "age", pt.size = 1)

FeaturePlot(ecs, c("Sntn"), cols = c("lightgrey", "red"), order = TRUE, pt.size = 0.8, split.by = "age")


## Clustering
ecs <- FindNeighbors(object = ecs, dims= 1:13, k.param = 10)

ecs <- FindClusters(object = ecs, graph.name = NULL,
                    modularity.fxn = 1, initial.membership = NULL, weights = NULL,
                    node.sizes = NULL, resolution = seq(0.1, 1, 0.1), algorithm = 1, n.start = 1000,
                    n.iter = 10, random.seed = 0, temp.file.location = NULL,
                    edge.file.name = NULL, verbose = TRUE)

# Clustering tree
clustree(ecs, prefix = "RNA_snn_res.")

# Colour the nodes by the expression of a specific gene
clustree(ecs, node_colour = "Sntn", node_colour_aggr = "median")

# Set cell identities
Idents(object = ecs) <- ecs[["RNA_snn_res.0.5"]]
my_cols_ecs <- DiscretePalette(n = 12, palette = "polychrome")
DimPlot(object = ecs, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.8, cols = my_cols_ecs)
DimPlot(object = ecs, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.8, cols = my_cols_ecs, split.by = "age")

DimPlot(ecs, reduction = "umap", label = FALSE, pt.size = 0.6, split.by = "age", cols = my_cols_ecs)
DimPlot(ecs, reduction = "umap", label = TRUE, pt.size = 0.6, cols = my_cols_ecs)


# Find markers
res08_ecs <- FindAllMarkers(object = ecs, genes.use = NULL, only.pos = TRUE,
                            logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
res05_ecs <- FindAllMarkers(object = ecs, genes.use = NULL, only.pos = TRUE,
                            logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)

DimPlot(ecs, reduction = "umap", label = TRUE, pt.size = 0.6, cols = my_cols_ecs)


res08_c2_vs_c8 <- FindMarkers(ecs, ident.1 = 2, ident.2 = 8, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res08_c8_vs_c2 <- FindMarkers(ecs, ident.1 = 8, ident.2 = 2, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)

res05_c3_vs_c0 <- FindMarkers(ecs, ident.1 = 3, ident.2 = 0, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
res05_c0_vs_c3 <- FindMarkers(ecs, ident.1 = 0, ident.2 = 3, only.pos = TRUE,
                              logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1) # 3 are activated cells


FeaturePlot(ecs, c("Bhmt"), cols = c("lightgrey", "red"), order = TRUE, pt.size = 0.8)


# Merge clusters 0, 3, and 4 (Lateral, immature)
my_cols_ecs <- DiscretePalette(n = 11, palette = "polychrome")
DimPlot(ecs, reduction = "umap", label = FALSE, pt.size = 0.8, cols = my_cols_ecs)
DimPlot(ecs, reduction = "umap", label = FALSE, pt.size = 0.8, cols = my_cols_ecs, split.by = "age")

# Merge clusters 1, 2, 5 and 6 (Lateral, immature)
Idents(ecs) <- ecs$RNA_snn_res.0.5
ecs <- RenameIdents(ecs, `0` = "0", `1` = "1", `2` = "1",
                    `3` = "2", `4` = "3", `5` = "1",
                    `6` = "1", `7` = "0", `8` = "4")

DimPlot(ecs, reduction = "umap", label = FALSE, pt.size = 1, cols = my_cols_ecs)

ecs$merged_clusters <- Idents(ecs)

table(ecs$merged_clusters)

Idents(ecs) <- ecs$merged_clusters
table(Idents(ecs))


# Assign cell type
ecs <- RenameIdents(ecs, `0` = "Lateral (mature)", `1` = "Lateral (immature)", `2` = "Lateral (activated)",
                    `3` = "Ventral", `4` = "Dorsal")

ecs$cell_type <- Idents(ecs)

ecs$cell_type <- factor(x = ecs$cell_type, levels = c("Lateral (immature)", "Lateral (mature)", "Lateral (activated)",
                                                      "Dorsal",
                                                      "Ventral"))

Idents(ecs) <- ecs$cell_type

DimPlot(ecs, reduction = "umap", label = FALSE, pt.size = 1)
DimPlot(ecs, reduction = "umap", label = FALSE, pt.size = 1, cols = my_cols_ecs, split.by = "age")

names(ecs[[]])

# Figure 3C
DimPlot(object = ecs, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 1,
        cols = c( "Lateral (immature)" = "#fbd347", "Lateral (activated)" = "#f1443d",
                  "Lateral (mature)" = "#fdae46",
                  "Ventral" = "#43a5f2", "Dorsal" = "#a647b8")) + NoAxes()


FeaturePlot(ecs, features = "Sntn", cols = c("lightgrey", "red"), order = TRUE, pt.size = 0.8, split.by = "age")



### Age-related changes in ependymal cells
table(Idents(ecs))

Idents(ecs) <- ecs$cell_type

# Append age to cell identity
ecs$cell_type_age <- paste(Idents(ecs), ecs$age, sep = "_")

table(ecs$cell_type_age)
table(ecs$cell_type) # pass this to Ident to remove the age tag

Idents(ecs) <- ecs$cell_type_age

table(Idents(ecs))

Idents(ecs) <- ecs$cell_type

saveRDS(ecs, file = "~/path/to/ageing_ecs_cell_types.rds")
ecs <- readRDS(file = "~/path/to/ageing_ecs_cell_types.rds") # Aggregated dataset



### Preparing the processed data for ArrayExpress
ecs
ecs_counts <- GetAssayData(object = ecs, slot = "counts")
write.table(ecs_counts, file = "~/path/to/ependymal_cells_ageing_raw_counts.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE)
ecs_metadata <- ecs[[]]
write.table(ecs_metadata, file = "~/path/to/ependymal_cells_ageing_metadata.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE)

ecs_umap <- Embeddings(ecs, reduction = c("umap"))
head(ecs_umap)
write.table(ecs_umap, file = "~/path/to/ependymal_cells_ageing_umap.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE)


# Split by sample and retrieve the count matrix and object metadata
table(ecs$orig.ident)


## All cells in the ageing dataset
allcells <- readRDS(file = "~/path/to/allcells_cell_types.rds")

DimPlot(allcells) + NoLegend()

# Extract raw counts and metadata
allcells_counts <- GetAssayData(object = allcells, slot = "counts")
write.table(allcells_counts, file = "~/path/to/all_cells_ageing_raw_counts.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE)
allcells_metadata <- allcells[[]]
write.table(allcells_metadata, file = "~/path/to/all_cells_ageing_metadata.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE)
allcells_umap <- Embeddings(allcells, reduction = c("umap"))
head(allcells_umap)
write.table(allcells_umap, file = "~/path/to/all_cells_ageing_umap.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE)

####






### Age-related gene expression changes in the ependymal cell population

Idents(ecs) <- ecs$age
table(Idents(ecs))

# Differentially regulated genes with ageing
ageing_genes <- FindMarkers(ecs, ident.1 = "old", ident.2 = "young",
                            only.pos = FALSE, logfc.threshold = 0.2, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
write.csv(ageing_genes, file = "~/path/to/ependymal_cells_age-regulated_genes.csv")


add_gene_names <- function(x) {
  
  data <- x %>%
    rownames_to_column() %>%
    rename(gene = rowname) #%>% print()
}


old_up <- ageing_genes %>% add_gene_names() %>% top_n(10, avg_logFC) #%>% print()
old_down <- ageing_genes %>% add_gene_names() %>% top_n(-10, avg_logFC) #%>% print()


Idents(ecs) <- ecs$age
young_vs_old <- FindMarkers(ecs, ident.1 = "young", ident.2 = "old",
                            only.pos = TRUE, logfc.threshold = 0.2, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)
old_vs_young <- FindMarkers(ecs, ident.1 = "old", ident.2 = "young",
                            only.pos = TRUE, logfc.threshold = 0.2, test.use = "wilcox", min.pct = 0.25, pseudocount.use = 1)

Idents(ecs) <- ecs$cell_type_age

### How do specific ependymal cell subtypes change with age?
mature_young <- FindMarkers(ecs, ident.1 = "Lateral (mature)_young", ident.2 = "Lateral (mature)_old",
                          only.pos = TRUE, logfc.threshold = 0.2, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
mature_old <- FindMarkers(ecs, ident.1 = "Lateral (mature)_old", ident.2 = "Lateral (mature)_young",
                        only.pos = TRUE, logfc.threshold = 0.2, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)

ventral_young <- FindMarkers(ecs, ident.1 = "Ventral_young", ident.2 = "Ventral_old",
                             only.pos = TRUE, logfc.threshold = 0.2, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
ventral_old <- FindMarkers(ecs, ident.1 = "Ventral_old", ident.2 = "Ventral_young",
                           only.pos = TRUE, logfc.threshold = 0.2, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)

dorsal_young <- FindMarkers(ecs, ident.1 = "Dorsal_young", ident.2 = "Dorsal_old",
                            only.pos = TRUE, logfc.threshold = 0.2, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
dorsal_old <- FindMarkers(ecs, ident.1 = "Dorsal_old", ident.2 = "Dorsal_young",
                          only.pos = TRUE, logfc.threshold = 0.2, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)

immature_young <- FindMarkers(ecs, ident.1 = "Lateral (immature)_young", ident.2 = "Lateral (immature)_old",
                                 only.pos = TRUE, logfc.threshold = 0.2, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)
immature_old <- FindMarkers(ecs, ident.1 = "Lateral (immature)_old", ident.2 = "Lateral (immature)_young",
                               only.pos = TRUE, logfc.threshold = 0.2, test.use = "wilcox", min.pct = 0.15, pseudocount.use = 1)


# Common genes enriched in cells from young mice
mature_young <- tibble::rownames_to_column(mature_young, var = 'gene')
dorsal_young <- tibble::rownames_to_column(dorsal_young, var = 'gene')
ventral_young <- tibble::rownames_to_column(ventral_young, var = 'gene')
immature_young <- tibble::rownames_to_column(immature_young, var = 'gene')

common_young <- intersect(intersect(intersect(mature_young$gene, dorsal_young$gene),ventral_young$gene), immature_young$gene)
common_young <- intersect(mature_young$gene, dorsal_young$gene)

# Common genes enriched in cells from old mice
mature_old <- tibble::rownames_to_column(mature_old, var = 'gene')
dorsal_old <- tibble::rownames_to_column(dorsal_old, var = 'gene')
ventral_old <- tibble::rownames_to_column(ventral_old, var = 'gene')
immature_old <- tibble::rownames_to_column(immature_old, var = 'gene')


FeaturePlot(ecs, features = common_old, pt.size = 1, order = TRUE, cols = c('lightgrey', 'red'), split.by = 'age')
FeaturePlot(ecs, features = "Caly", pt.size = 1, order = TRUE, cols = c('lightgrey', 'red'), split.by = 'age')



## How does the celullar composition of the ependymal cell population change with age?

# Calculating ependymal cell subtype fractions
library(tidyr)

Idents(ecs) <- "cell_type" # pass cluster ID as cell identity
table(Idents(ecs))

cell_fractions <- table(Idents(ecs), ecs@meta.data$orig.ident)  %>% as.data.frame() %>% print()
cell_fractions <- plyr::rename(cell_fractions, c("Var1" = "cluster", "Var2" = "mouse", "Freq" = "value"))
cell_fractions <- as_tibble(cell_fractions)

head(cell_fractions)

str(cell_fractions)

# Add age info
mouse_ids <- c("10X01_1", "10X01_2", "10X03_1", "10X03_2", "10X05_1", "10X05_2", "10X06_1", "10X06_2")
mouse_age <- c("young", "old", "young", "old", "old", "young", "old", "young")

age <- NULL

age <- plyr::mapvalues(x = cell_fractions$mouse, from = mouse_ids, to = mouse_age)
cell_fractions <- add_column(cell_fractions, age, .after = "mouse")

# Calculate fractions
#cell_fractions %>%  group_by(mouse) %>% dplyr::summarise(vs = value/sum(value))

cell_fractions <- cell_fractions %>% group_by(mouse) %>% 
  mutate(total_cells = sum(value), cluster_pct = value/total_cells*100) %>% print()

cell_pcts <- cell_fractions %>% group_by(mouse) %>% 
  dplyr::mutate(total_cells = sum(value), 
                cluster_pct = value / total_cells*100)


# # Summarise the data
stats <- data_summary(cell_fractions, varname="cluster_pct",
                      groupnames=c("cluster", "age"))
head(stats)

# Bar plot with all the values superimposed
p <- ggplot(cell_pcts, aes(cluster, cluster_pct, fill = age))

p + geom_bar(position = 'dodge', stat = 'summary', fun.y = 'mean', width = 0.8) +
  geom_errorbar(stat = 'summary', position = position_dodge(0.8), width = 0.25) +
  geom_point(aes(x = cluster),
             position =  position_jitterdodge(jitter.width = 0.5, 
                                              dodge.width=0.8),
             size=1, alpha=0.8) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45,hjust=1)) # Figure 3D



# Statistical analysis
table(Idents(ecs))

mature_young_fractions <- dplyr::filter(cell_fractions, cluster == "Lateral (mature)", age == "young") %>%
  pull(cluster_pct) %>% print()

mature_old_fractions <- dplyr::filter(cell_fractions, cluster == "Lateral (mature)", age == "old") %>%
  pull(cluster_pct) %>% print()

wilcox.test(mature_young_fractions, mature_old_fractions,
            alternative = "two.sided",
            pool.sd = TRUE,
            paired = FALSE,
            var.equal = TRUE)


immature_young_fractions <- dplyr::filter(cell_fractions, cluster == "Lateral (immature)", age == "young") %>%
  pull(cluster_pct) %>% print()

immature_old_fractions <- dplyr::filter(cell_fractions, cluster == "Lateral (immature)", age == "old") %>%
  pull(cluster_pct) %>% print()

wilcox.test(immature_young_fractions, immature_old_fractions,
            alternative = "two.sided",
            pool.sd = TRUE,
            paired = FALSE,
            var.equal = TRUE)

activated_young_fractions <- dplyr::filter(cell_fractions, cluster == "Lateral (activated)", age == "young") %>%
  pull(cluster_pct) %>% print()

activated_old_fractions <- dplyr::filter(cell_fractions, cluster == "Lateral (activated)", age == "old") %>%
  pull(cluster_pct) %>% print()

wilcox.test(activated_young_fractions, activated_old_fractions,
            alternative = "two.sided",
            pool.sd = TRUE,
            paired = FALSE,
            var.equal = TRUE)

dorsal_young_fractions <- dplyr::filter(cell_fractions, cluster == "Dorsal", age == "young") %>%
  pull(cluster_pct) %>% print()

dorsal_old_fractions <- dplyr::filter(cell_fractions, cluster == "Dorsal", age == "old") %>%
  pull(cluster_pct) %>% print()

wilcox.test(dorsal_young_fractions, dorsal_old_fractions,
            alternative = "two.sided",
            pool.sd = TRUE,
            paired = FALSE,
            var.equal = TRUE)

ventral_young_fractions <- dplyr::filter(cell_fractions, cluster == "Ventral", age == "young") %>%
  pull(cluster_pct) %>% print()

ventral_old_fractions <- dplyr::filter(cell_fractions, cluster == "Ventral", age == "old") %>%
  pull(cluster_pct) %>% print()

wilcox.test(ventral_young_fractions, ventral_old_fractions,
            alternative = "two.sided",
            pool.sd = TRUE,
            paired = FALSE,
            var.equal = TRUE)



# Plotting gene expression patterns on UMAP plots with ggplot2
single_cell_object <- ecs

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
) %>% mutate(map(gene_name, ~extract_plot_data(ecs, .x))) %>%
  unnest() %>% print()

plot_title <- plot_df$gene_name %>% unique()

ggplot(plot_df %>% filter(value == 0), aes(x, y)) +
  geom_point(color = 'lightgrey', size = 1.1) +
  geom_point(data = plot_df %>% filter(value != 0), aes(x, y, color = value), size = 1) +
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



# Plotting young and old using the viridis coulour palette
FeaturePlot(ecs, features = c("Sntn"), split.by = "age", cols = c('grey90',viridis::viridis(1000)),
            pt.size = 0.6, reduction = "umap", order = TRUE, combine = TRUE)


