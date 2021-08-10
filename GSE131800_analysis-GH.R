library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggpubr)
library(hdf5r)
library(GEOquery)
# working with the GEO record
#seems to download to a tmp directory. and if you want the series matrix file you 
GSE131800 <- getGEO(GEO = 'GSE131800')
#look at the structure of the file. 
str(GSE131800)
#create a df of sample characteristics/ phenoData
GSE131800_series_matrix <- GSE131800$GSE131800_series_matrix.txt.gz@phenoData@data
write.csv(GSE131800_series_matrix, file = 'GSE131800_series_matrix.csv', row.names = FALSE)
# useGEOquery to retrieve supplemental files for GSE131800
getGEOSuppFiles(GEO = 'GSE131800', makeDirectory = FALSE, baseDir = '.')
# untar the supplemental file
untar(tarfile = 'GSE131800_RAW.tar', exdir = '.')

library(stringr)

dataFiles <- list.files(".", pattern = 'filtered', full.names = TRUE)
prefixes <- str_extract(dataFiles, pattern = 'WT[:alnum:]*')
#associate sampleID with corresponding file
names(dataFiles) <- prefixes
# Read in first data file and create Seurat object WTd0
h5Input.data_WTd0 <- Read10X_h5(filename = dataFiles[1])
lung_WTd0 <- CreateSeuratObject(counts = h5Input.data_WTd0$`Gene Expression`, min.cells = 3, min.genes = 200, project = "WTd0")
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes; no mt genes found
lung_WTd0[["percent.mt"]] <- PercentageFeatureSet(object = lung_WTd0, pattern = "mt-")
VlnPlot(object = lung_WTd0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(lung_WTd0, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(lung_WTd0, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(lung_WTd0, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data
lung_WTd0 <- subset(lung_WTd0,  subset = nFeature_RNA > 1000 & nFeature_RNA < 3000)
lung_WTd0 <- NormalizeData(lung_WTd0, verbose = FALSE)
lung_WTd0 <- FindVariableFeatures(lung_WTd0, selection.method = "vst", nfeatures = 2000)

# Read in first data file and create Seurat object WTBleoD21
h5Input.data_WTBleoD21 <- Read10X_h5(filename = dataFiles[3])
lung_WTBleoD21 <- CreateSeuratObject(counts = h5Input.data_WTBleoD21$`Gene Expression`, min.cells = 3, min.genes = 200, project = "WTBleoD21")
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes; no mt genes found
lung_WTBleoD21[["percent.mt"]] <- PercentageFeatureSet(object = lung_WTBleoD21, pattern = "mt-")
VlnPlot(object = lung_WTBleoD21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(lung_WTBleoD21, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(lung_WTBleoD21, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(lung_WTBleoD21, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data
lung_WTBleoD21 <- subset(lung_WTBleoD21,  subset = nFeature_RNA > 1000 & nFeature_RNA < 3000)
lung_WTBleoD21 <- NormalizeData(lung_WTBleoD21, verbose = FALSE)
lung_WTBleoD21 <- FindVariableFeatures(lung_WTBleoD21, selection.method = "vst", nfeatures = 2000)

# save seurat objects
saveRDS(lung_WTd0, file = "lung_WTd0.RDS")
saveRDS(lung_WTBleoD21, file = 'lung_WTBleoD21.RDS')
# load objects
lung_WTd0 <- readRDS("lung_WTd0.RDS")
lung_WTBleoD21 <- readRDS("lung_WTBleoD21.RDS")
# Assign Groups
lung_WTd0$group = "lung_WTd0"
lung_WTBleoD21$group = "lung_WTBleoD21"
lung_WTBleoD21$group
#Find Anchors for data
AEC.anchors <- FindIntegrationAnchors(object.list = list(lung_WTd0, lung_WTBleoD21), dims = 1:30)

#Combine data
AEC.combined <- IntegrateData(anchorset = AEC.anchors, dims = 1:30)


#Perform an integrated analysis
DefaultAssay(object = AEC.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
AEC.combined <- ScaleData(object = AEC.combined, features=rownames(AEC.combined), verbose = FALSE)
# Run the standard workflow for visualization and clustering
AEC.combined <- RunPCA(object = AEC.combined, npcs = 50, verbose = FALSE)

#Determine the ??dimensionality?? of the dataset
#Method 1 Jackstrawp
AEC.combined <- JackStraw(AEC.combined, num.replicate = 100)
AEC.combined <- ScoreJackStraw(AEC.combined, dims = 1:20)
JackStrawPlot(AEC.combined, dims = 1:20)
# Method 2
ElbowPlot(AEC.combined, ndims = 30)

#After the above step, set dims to 15
#t-SNE and Clustering
AEC.combined<- FindNeighbors(object = AEC.combined, reduction = "pca", dims = 1:15)
AEC.combined <- FindClusters(AEC.combined, resolution = 0.2)
AEC.combined <- RunUMAP(object = AEC.combined, reduction = "pca", dims = 1:15)


DefaultAssay(AEC.combined) <- "RNA"

# Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
AEC.combined <- RunUMAP(AEC.combined, dims = 1:10)
AEC.combined <- RunTSNE(AEC.combined, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(AEC.combined, reduction = "umap")
DimPlot(AEC.combined, reduction = "tsne")
# save AEC.combined to return a complete object following Seurat tutorial instructions
saveRDS(AEC.combined, file = "AEC.combined.rds")
# find markers for every cluster compared to all remaining cells, report only the positive ones
AEC.combined.markers <- FindAllMarkers(AEC.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(AEC.combined.markers, file = 'AEC.combined.markers_GSE131800.csv',row.names = FALSE)
top_AEC.combined.markers <- AEC.combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)
top20_AEC.combined.markers <- AEC.combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
# find all markers of cluster 12
cluster12.markers <- FindMarkers(AEC.combined, ident.1 = 12, min.pct = 0.25)
head(cluster12.markers, n = 5)
# Visualize marker expression
VlnPlot(AEC.combined, features = c("Rgcc", "Col3a1", "Col1a2", "Col1a1"))
# you can plot raw counts as well
VlnPlot(AEC.combined, features = c("Rgcc"), slot = "counts", log = TRUE)
# feature plot
FeaturePlot(AEC.combined, features = c("Rgcc", "Col3a1", "Col1a2", "Col1a1"))
# DimPlot and Feature Plots from afternoon 01JUL2021
DimPlot(AEC.combined, reduction = "tsne", label = TRUE)
FeaturePlot(AEC.combined, features = c("Acta2", 'Col1a1','Cdh1', 'Pecam1', 'Rgcc'),reduction = 'tsne',label = TRUE)
FeaturePlot(AEC.combined, features = c("Acta2", 'Col1a1','Cdh1', 'Pecam1', 'Rgcc', 'Il33'),reduction = 'tsne',label = TRUE)
# read in cluster annotation from SCSA to write to csv

# Figures for MS
# Panel A
tsnePalette <- CustomPalette(low = "green", high = "red", mid = 'blue', k = 15)
panelA <- DimPlot(AEC.combined, cols = tsnePalette, reduction = "tsne", label = TRUE, label.size = 8)
# Panel B
panelB <- FeaturePlot(AEC.combined, features = "Rgcc", reduction = 'tsne',label = TRUE, split.by = "group")


# Marker Gene for clusters
# Annotation of epithelial1 cells - Cluster 0
epithelial1 <- c("Cdh1", "Epcam", "Muc1", 'Sftpa1',"Sftpb","Sftpc",'Sftpd')
FeaturePlot(AEC.combined, features = epithelial1, reduction = 'tsne',label = TRUE)
# Annotation of epithelial2 cells - Cluster 9
epithelial2 <- c("Cdh1", "Epcam","Muc1", 'Hopx', 'Sfn', "Krt7", "Krt19")
FeaturePlot(AEC.combined, features = epithelial2, reduction = 'tsne',label = TRUE)
# Annotation of endothelial1 cells - Cluster 8
endothelial1 <- c("Pecam1", "Cdh5", "Cldn5", "Calcrl", "Ly6a",'Ly6c1')
FeaturePlot(AEC.combined, features = endothelial1, reduction = 'tsne',label = TRUE)
# Annotation of endothelial2 cells - Cluster 11
endothelial2 <- c("Pecam1", "Cdh5", "Cldn5", "Mmrn1")
FeaturePlot(AEC.combined, features = endothelial2, reduction = 'tsne',label = TRUE)
# Annotation of myofibroblasts - cluster 13
smoothMuscleCells <- c('Myh11','Actc1',"Cnn1",'Actg2','Acta2')
FeaturePlot(AEC.combined, features = smoothMuscleCells, reduction = 'tsne',label = TRUE)
# Annotation of fibroblasts1 - cluster 1
fibroblasts1 <- c('Col1a1', 'Col1a2', 'Col3a1','Col13a1','Cxcl14')
FeaturePlot(AEC.combined, features = fibroblasts1, reduction = 'tsne',label = TRUE)
# Annotation of fibroblasts2 - cluster 6
fibroblasts2 <- c('Col1a1', 'Col1a2', 'Col3a1','Col14a1', 'Ccl11')
FeaturePlot(AEC.combined, features = fibroblasts2, reduction = 'tsne',label = TRUE)
# Annotation of T Lymphocytes - Cluster 5
T_lymphocytes <- c('Cd3d', 'Cd3e', 'Cd3g', 'Il7r', 'Icos')
FeaturePlot(AEC.combined, features = T_lymphocytes, reduction = 'tsne',label = TRUE)
# Annotation of B Lymphocytes - Cluster 7
B_lymphocytes <- c('Ms4a1', 'Cd19','Cd79a','Mzb1','H2-DMb2')
FeaturePlot(AEC.combined, features = B_lymphocytes, reduction = 'tsne',label = TRUE)
# Annotation of Mesotheliocytes - Cluster 12
Mesothelial <- c('Msln')
FeaturePlot(AEC.combined, features = Mesothelial, reduction = 'tsne',label = TRUE)
# Annotation of Dendritic Cells - Cluster 4
Dendritic <- c('H2-DMb2', 'Ccl22', 'Cd80', 'Ccl17')
FeaturePlot(AEC.combined, features = Dendritic, reduction = 'tsne',label = TRUE)
# Annotation of Macrophages1 - Cluster 3
Macrophages1 <- c('Cd14','Cd86', 'Tlr2', 'Cd68', 'Ifitm6')
FeaturePlot(AEC.combined, features = Macrophages1, reduction = 'tsne',label = TRUE)
# Annotation of Macrophages2 - Cluster 2
Macrophages2 <- c('Cd14', 'Cd68', 'Tlr2', 'Spp1','Mrc1')
FeaturePlot(AEC.combined, features = Macrophages2, reduction = 'tsne',label = TRUE)
# Annotation of Dividing cells - Cluster 10
Dividing <- c('Pclaf', 'Top2a','Mki67','Cdca3', 'Ccnb2')
FeaturePlot(AEC.combined, features = Dividing, reduction = 'tsne',label = TRUE)
# Annotation of Ciliated Cells - Cluster 14
Ciliated <- c('Dynlrb2','Rsph1','Riiad1')
FeaturePlot(AEC.combined, features = Ciliated, reduction = 'tsne',label = TRUE)

# DimPlot of all cluseters
GSE131800_DimPlot <-
  DimPlot(
    AEC.combined,
    cols = rainbow(n = 15),
    reduction = "tsne",
    label = TRUE,
    label.size = 8
  )
# DimPlot without labels
GSE131800_DimPlot_nolabel <-
  DimPlot(
    AEC.combined,
    cols = rainbow(n = 15),
    reduction = "tsne",
    label = FALSE,
    label.size = 8
  )

# save dimplot 600dpi
ggsave(filename = "GSE131800_DimPlot.png",
       plot = GSE131800_DimPlot,
       dpi = 600,
       width = 4.25, 
       height = 3.65
)
ggsave(filename = "GSE131800_DimPlot_nolabel.png",
       plot = GSE131800_DimPlot_nolabel,
       dpi = 600,
       width = 4.25, 
       height = 3.65
)
# Split featureplot RGCC
RGCC_FeaturePlot <- FeaturePlot(AEC.combined, features = "Rgcc", cols = CustomPalette(low = "#BBBBBB",high = "red"), reduction = 'tsne',label = FALSE, split.by = 'group')
# save split RGCC feature plots
ggsave(filename = "GSE131800_Rgcc-SplitFeaturePlot.png",
       plot = RGCC_FeaturePlot,
       dpi = 600,
       width = 7.25, 
       height = 3.65
)

RgccVlnPlot <- VlnPlot(object = AEC.combined, features = "Rgcc", idents = 14:0, split.by = 'group', cols = c("#ffb5b5", "#d1f5ff"))
ggsave(filename = "GSE131800_RGCC-ViolinPlot.png",
       plot = RgccVlnPlot,
       dpi = 600
)
# Split.plot Rgcc Violin Plot
SplitRgccVlnPlot <- VlnPlot(AEC.combined, features = "Rgcc", pt.size = 0.1, idents = 0:14, split.by = "group", split.plot = TRUE, y.max = 7)
ggsave(filename = "GSE131800_Split-RGCC-ViolinPlot.png",
       plot = SplitRgccVlnPlot,
       dpi = 600
)
# print violin plots with significance test
for(i in clusters){
  junk <- VlnPlot(AEC.combined, features = "Rgcc", 
                  pt.size = 0.1, idents = i, 
                  group.by = "group") + 
    stat_compare_means(comparisons = comparisons, label = "p.signif")
  p <- VlnPlot(AEC.combined, features = "Rgcc", 
               pt.size = 0.1, idents = i, 
               group.by = "group", y.max = 1 + max(junk[[1]]$data[["Rgcc"]])) + 
    stat_compare_means(comparisons = comparisons, label = "p.signif") 
  ggsave(paste0("SigVlnPlot_Cluster", i, ".png"), p, dpi = 600)
}