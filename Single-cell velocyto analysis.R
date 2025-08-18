library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(dplyr)
library(patchwork)

DCA<-list()
ldat1 <- ReadVelocity(file = "~/Desktop/DCA/CellRanger210615/1SA/velocyto/1SA.loom")
DCA$bm1 <- as.Seurat(x = ldat1, min.cells=3)
DCA$bm1 <- SCTransform(object = DCA$bm1, assay = "spliced")
DCA$bm1 <- AddMetaData(DCA$bm1  , "1SA", col.name = "stim")

ldat2 <- ReadVelocity(file = "~/Desktop/DCA/CellRanger210615/2ACS/velocyto/2ACS.loom")
DCA$bm2 <- as.Seurat(x = ldat2)
DCA$bm2 <- SCTransform(object = DCA$bm2, assay = "spliced")
DCA$bm2 <- CreateSeuratObject(counts = DCA$bm1,min.cells = 3)
DCA$bm2 <- AddMetaData(DCA$bm2  , "2ACS", col.name = "stim")

DCA$bm1[["percent.mt"]] <- PercentageFeatureSet(DCA$bm1, pattern = "MT")
VlnPlot(DCA$bm1, features = c("nFeature_spliced", "nCount_spliced", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(DCA$bm1, feature1 = "nCount_spliced", feature2 = "percent.mt")
plot2 <- FeatureScatter(DCA$bm1, feature1 = "nCount_spliced", feature2 = "nFeature_spliced")
plot1 + plot2

DCA$bm1<- subset(DCA$bm1, subset = nFeature_spliced > 500 & nFeature_spliced < 5000 & percent.mt < 8)

DCA$bm2[["percent.mt"]] <- PercentageFeatureSet(DCA$bm2, pattern = "MT")
VlnPlot(DCA$bm2, features = c("nFeature_spliced", "nCount_spliced", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(DCA$bm2, feature1 = "nCount_spliced", feature2 = "percent.mt")
plot2 <- FeatureScatter(DCA$bm2, feature1 = "nCount_spliced", feature2 = "nFeature_spliced")
plot1 + plot2

DCA$bm2<- subset(DCA$bm2, subset = nFeature_spliced > 500 & nFeature_spliced < 5000 & percent.mt < 8)


# normalize and identify variable features for each dataset independently
all.genes.1 <- rownames(DCA$bm1)
DCA$bm1 <- ScaleData(DCA$bm1, features = all.genes.1)

all.genes.2 <- rownames(DCA$bm2)
DCA$bm2 <- ScaleData(DCA$bm2, features = all.genes.2)


all.genes <- rownames(DCA["1SA"]) +  rownames(DCA["2ACS"])

DCA <- lapply(X = DCA, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- ScaleData(x, features = all.genes)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = DCA)

# Integration
DCA.v.anchors <- FindIntegrationAnchors(object.list = DCA, anchor.features = features)

# this command creates an 'integrated' data assay
DCA.v.combined <- IntegrateData(anchorset =DCA.v.anchors)



DefaultAssay(DCA.v.combined) <- "integrated"

DCA.v.combined <- ScaleData(DCA.v.combined, verbose = FALSE)
DCA.v.combined <- RunPCA(DCA.v.combined, npcs = 30, verbose = FALSE)
DCA.v.combined <- RunUMAP(DCA.v.combined, reduction = "pca", dims = 1:30)
DCA.v.combined <- FindNeighbors(DCA.v.combined, reduction = "pca", dims = 1:30)
DCA.v.combined <- FindClusters(DCA.v.combined, resolution = 0.5)

p1 <- DimPlot(DCA.v.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(DCA.v.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 
p2

p1 <- DimPlot(DCA.v.combined, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(DCA.v.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2


# Myeloid----
VlnPlot(DCA.v.combined, features = c("CD68","CD14","CSF1R","IL1B","S100A8", "S100A9") )
DCA.v.combined_1 <- subset(DCA.v.combined, idents = c(2,5,8,9))
DCA.v.combined_1 <- ScaleData(DCA.v.combined_1, verbose = FALSE)
DCA.v.combined_1<- RunPCA(DCA.v.combined_1, npcs = 30, verbose = FALSE)
DCA.v.combined_1 <- RunUMAP(DCA.v.combined_1 , reduction = "pca", dims = 1:30)
DCA.v.combined_1 <- FindNeighbors(DCA.v.combined_1 , reduction = "pca", dims = 1:30)
DCA.v.combined_1 <- FindClusters(DCA.v.combined_1, resolution = 0.5)
p1 <- DimPlot(DCA.v.combined_1, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(DCA.v.combined_1, reduction = "umap", label = TRUE, repel = TRUE)
p2
VlnPlot(DCA.v.combined_1, features = c("CD68","CD14","CSF1R","IL1B","S100A8", "S100A9") )
VlnPlot(DCA.v.combined_1, features = c("HDC","KIT","TPSAB1","GZMB","CCL5","TPSB2") )
VlnPlot(DCA.v.combined_1, features = c("TET2","DNMT3A","ASXL1","JAK2","SF3B1","SRSF2") )
VlnPlot(DCA.v.combined_1, features = c("HLA-DQB1","HLA-DQB1","HLA-DPA1","TREM2","CLEC10A","FCER1A") )
VlnPlot(DCA.v.combined_1, features = c("CD9","MRC1","FCER2","CCL22","AREG","EREG") )
VlnPlot(DCA.v.combined_1, features = c("TNF","IL1B","CXCL8","CXCL3","CCL3","CXCL10") )
VlnPlot(DCA.v.combined_1, features = c("TLR4","TLR2","TLR7","TNF","KLF4", "NFkB") )
VlnPlot(DCA.v.combined_1, features = c("ABCG1","ABCA1","IL18","CD9", "OLR1","TREM2") )
FeaturePlot(DCA.v.combined_1,features = c("TLR4"))

DCA.v.v.combined_1 <- subset(DCA.v.combined_1, idents = c(0,1,2,3))
DCA.v.v.combined_1 <- ScaleData(DCA.v.v.combined_1, verbose = FALSE)
DCA.v.v.combined_1<- RunPCA(DCA.v.v.combined_1, npcs = 30, verbose = FALSE)
DCA.v.v.combined_1 <- RunUMAP(DCA.v.v.combined_1 , reduction = "pca", dims = 1:30)
DCA.v.v.combined_1 <- FindNeighbors(DCA.v.v.combined_1 , reduction = "pca", dims = 1:30)
DCA.v.v.combined_1 <- FindClusters(DCA.v.v.combined_1, resolution = 0.5)
p1 <- DimPlot(DCA.v.v.combined_1, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(DCA.v.v.combined_1, reduction = "umap", label = False, repel = TRUE)
p2

VlnPlot(DCA.v.v.combined_1, features = c("CD68","CD14","CSF1R","IL1B","TNF", "TREM2") )
VlnPlot(DCA.v.v.combined_1, features = c("HDC","KIT","TPSAB1","GZMB","CCL5","TPSB2") )
VlnPlot(DCA.v.v.combined_1, features = c("TET2","DNMT3A","ASXL1","JAK2","SF3B1","SRSF2") )
VlnPlot(DCA.v.v.combined_1, features = c("HLA-DQA1","HLA-DQB1","HLA-DPA1","TREM2","CLEC10A","FCER1A") )
VlnPlot(DCA.v.v.combined_1, features = c("CD9","MRC1","FCER2","CCL22","AREG","EREG") )
VlnPlot(DCA.v.v.combined_1, features = c("IL1B","CXCL8","CXCL3","CCL3","CXCL10") )
VlnPlot(DCA.v.v.combined_1, features = c("TLR4","TLR2","TLR7","TNF","KLF4", "NFkB") )
VlnPlot(DCA.v.v.combined_1, features = c("APOE","ABCA1","CEBPB","CD9", "OLR1","TREM2") )
VlnPlot(DCA.v.v.combined_1, features = c("C1QA","C1QB","C1QC","TYMS","NFKB1","STAT4") )


DCA.v.v.combined_1$celltype <- Idents(DCA.v.v.combined_1)
DCA.v.v.combined_1$celltype.stim <- paste(Idents(DCA.v.v.combined_1), DCA.v.v.combined_1$stim, sep="_")
Idents(DCA.v.v.combined_1) <- "celltype.stim" 
levels(DCA.v.v.combined_1) 


DCA.v.v.combined_Myeloid_1SA<- subset(DCA.v.v.combined_1, idents = c("0_1SA","1_1SA","2_1SA","3_1SA"))
DCA.v.v.combined_Myeloid_2ACS<- subset(DCA.v.v.combined_1, idents = c("0_2ACS","1_2ACS","2_2ACS","3_2ACS"))





# velocity
DCA.v.v.combined_Myeloid_1SA <- RunVelocity(object = DCA.v.v.combined_Myeloid_1SA, deltaT = 1, kCells = 25, fit.quantile = 0.02)
# deltaT: amount of time to project the cell forward
# kCells: number of k nearest neighbors (NN) to use in slope calculation smoothing
# fit.quantile: perform gamma fit on a top/bottom quantiles of expression magnitudes

ident.colors <- c("#F8766D", "#ABA300", "#00A9FF","#0CB702")
names(x = ident.colors) <- levels(x = DCA.v.v.combined_Myeloid_1SA)
cell.colors <- ident.colors[Idents(object = DCA.v.v.combined_Myeloid_1SA)]
names(x = cell.colors) <- colnames(x = DCA.v.v.combined_Myeloid_1SA)


show.velocity.on.embedding.cor(emb = Embeddings(object = DCA.v.v.combined_Myeloid_1SA, reduction = "umap"),
                               vel = Tool(object = DCA.v.v.combined_Myeloid_1SA, slot = "RunVelocity"),
                               n = 200, scale = "sqrt",
                               cell.colors = ac(x = cell.colors, alpha = 0.5),
                               cex = 0.8, arrow.scale = 2, show.grid.flow = TRUE,
                               min.grid.cell.mass = 0.5, grid.n = 15, arrow.lwd = 0.5,
                               do.par = FALSE, cell.border.alpha = 0.1, n.cores = 1, xlim=c(-8, 8), ylim=c(-7, 7))



DCA.v.v.combined_Myeloid_2ACS <- RunVelocity(object = DCA.v.v.combined_Myeloid_2ACS, deltaT = 1, kCells = 25, fit.quantile = 0.02)
# deltaT: amount of time to project the cell forward
# kCells: number of k nearest neighbors (NN) to use in slope calculation smoothing
# fit.quantile: perform gamma fit on a top/bottom quantiles of expression magnitudes

ident.colors <- c("#00A9FF","#0CB702","#ABA300","#F8766D")
names(x = ident.colors) <- levels(x = DCA.v.v.combined_Myeloid_2ACS)
cell.colors <- ident.colors[Idents(object = DCA.v.v.combined_Myeloid_2ACS)]
names(x = cell.colors) <- colnames(x = DCA.v.v.combined_Myeloid_2ACS)

show.velocity.on.embedding.cor(emb = Embeddings(object = DCA.v.v.combined_Myeloid_2ACS, reduction = "umap"), 
                               vel = Tool(object = DCA.v.v.combined_Myeloid_2ACS, slot = "RunVelocity"), 
                               n = 200, scale = "sqrt", 
                               cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 2, show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, grid.n = 15, arrow.lwd = 0.5, 
                               do.par = FALSE, cell.border.alpha = 0.1, n.cores = 1, xlim=c(-8, 8), ylim=c(-7, 7) )


