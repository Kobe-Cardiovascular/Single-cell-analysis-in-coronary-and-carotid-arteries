library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
install.packages("ggplot2")
install.packages("ggpie")
library("ggpie")

# split the dataset into a list of two seurat objects (stim and CTRL)----
CARO.list<-list()
data.sym7<- Read10X(data.dir = "~/Desktop/Human_carotid_public/Cellranger/sample7_Symp")
data.asym8<- Read10X(data.dir = "~/Desktop/Human_carotid_public/Cellranger/sample8_Asymp")
data.sym9<- Read10X(data.dir = "~/Desktop/Human_carotid_public/Cellranger/sample9_Symp")
data.sym10<- Read10X(data.dir = "~/Desktop/Human_carotid_public/Cellranger/sample10_Symp")
data.sym11<- Read10X(data.dir = "~/Desktop/Human_carotid_public/Cellranger/sample11_Symp")
data.asym12<- Read10X(data.dir = "~/Desktop/Human_carotid_public/Cellranger/sample12_Asymp")
data.asym13<- Read10X(data.dir = "~/Desktop/Human_carotid_public/Cellranger/sample13_Asymp")
data.asym14<- Read10X(data.dir = "~/Desktop/Human_carotid_public/Cellranger/sample14_Asymp")
data.sym15<- Read10X(data.dir = "~/Desktop/Human_carotid_public/Cellranger/sample15_Symp")
data.asym16<- Read10X(data.dir = "~/Desktop/Human_carotid_public/Cellranger/sample16_Asymp")
data.asym17<- Read10X(data.dir = "~/Desktop/Human_carotid_public/Cellranger/sample17_Asymp")
data.sym18<- Read10X(data.dir = "~/Desktop/Human_carotid_public/Cellranger/sample18_Symp")
data.sym19<- Read10X(data.dir = "~/Desktop/Human_carotid_public/Cellranger/sample19_Symp")
data.asym20<- Read10X(data.dir = "~/Desktop/Human_carotid_public/Cellranger/sample20_Asymp")
data.sym21<- Read10X(data.dir = "~/Desktop/Human_carotid_public/Cellranger/sample21_Symp")
data.sym22<- Read10X(data.dir = "~/Desktop/Human_carotid_public/Cellranger/sample22_Symp")
data.sym23<- Read10X(data.dir = "~/Desktop/Human_carotid_public/Cellranger/sample23_Symp")
data.sym24<- Read10X(data.dir = "~/Desktop/Human_carotid_public/Cellranger/sample24_Symp")


data.sym19 <- data.sym19$`Gene Expression`
data.asym20 <- data.asym20$`Gene Expression`
data.sym21 <- data.sym21$`Gene Expression`
data.sym22 <- data.sym22$`Gene Expression`
data.sym23 <- data.sym23$`Gene Expression`
data.sym24 <- data.sym24$`Gene Expression`

CARO.list$sym7 <- CreateSeuratObject(counts = data.sym7, project = "carotid",
                                     min.cells = 3, min.features = 200)
CARO.list$sym7 <- AddMetaData(CARO.list$sym7 , "sym", col.name = "stim")
CARO.list$sym7 <- AddMetaData(CARO.list$sym7 , "sym7", col.name = "sample")

CARO.list$asym8 <- CreateSeuratObject(counts = data.asym8, project = "carotid",
                                      min.cells = 3, min.features = 200)
CARO.list$asym8 <- AddMetaData(CARO.list$asym8 , "asym", col.name = "stim")
CARO.list$asym8 <- AddMetaData(CARO.list$asym8 , "asym8", col.name = "sample")

CARO.list$sym9 <- CreateSeuratObject(counts = data.sym9, project = "carotid",
                                     min.cells = 3, min.features = 200)
CARO.list$sym9 <- AddMetaData(CARO.list$sym9 , "sym", col.name = "stim")
CARO.list$sym9 <- AddMetaData(CARO.list$sym9 , "sym9", col.name = "sample")

CARO.list$sym10 <- CreateSeuratObject(counts = data.sym10, project = "carotid",
                                      min.cells = 3, min.features = 200)
CARO.list$sym10 <- AddMetaData(CARO.list$sym10 , "sym", col.name = "stim")
CARO.list$sym10 <- AddMetaData(CARO.list$sym10 , "sym10", col.name = "sample")

CARO.list$sym11 <- CreateSeuratObject(counts = data.sym11, project = "carotid",
                                      min.cells = 3, min.features = 200)
CARO.list$sym11 <- AddMetaData(CARO.list$sym11 , "sym", col.name = "stim")
CARO.list$sym11 <- AddMetaData(CARO.list$sym11 , "sym11", col.name = "sample")

CARO.list$asym12 <- CreateSeuratObject(counts = data.asym12, project = "carotid",
                                       min.cells = 3, min.features = 200)
CARO.list$asym12 <- AddMetaData(CARO.list$asym12 , "asym", col.name = "stim")
CARO.list$asym12 <- AddMetaData(CARO.list$asym12 , "asym12", col.name = "sample")

CARO.list$asym13 <- CreateSeuratObject(counts = data.asym13, project = "carotid",
                                       min.cells = 3, min.features = 200)
CARO.list$asym13 <- AddMetaData(CARO.list$asym13 , "asym", col.name = "stim")
CARO.list$asym13 <- AddMetaData(CARO.list$asym13 , "asym13", col.name = "sample")

CARO.list$asym14 <- CreateSeuratObject(counts = data.asym14, project = "carotid",
                                       min.cells = 3, min.features = 200)
CARO.list$asym14 <- AddMetaData(CARO.list$asym14 , "asym", col.name = "stim")
CARO.list$asym14 <- AddMetaData(CARO.list$asym14 , "asym14", col.name = "sample")

CARO.list$sym15 <- CreateSeuratObject(counts = data.sym15, project = "carotid",
                                      min.cells = 3, min.features = 200)
CARO.list$sym15 <- AddMetaData(CARO.list$sym15 , "sym", col.name = "stim")
CARO.list$sym15 <- AddMetaData(CARO.list$sym15 , "sym15", col.name = "sample")

CARO.list$asym16 <- CreateSeuratObject(counts = data.asym16, project = "carotid",
                                       min.cells = 3, min.features = 200)
CARO.list$asym16 <- AddMetaData(CARO.list$asym16 , "asym", col.name = "stim")
CARO.list$asym16 <- AddMetaData(CARO.list$asym16 , "asym16", col.name = "sample")

CARO.list$asym17 <- CreateSeuratObject(counts = data.asym17, project = "carotid",
                                       min.cells = 3, min.features = 200)
CARO.list$asym17 <- AddMetaData(CARO.list$asym17 , "asym", col.name = "stim")
CARO.list$asym17 <- AddMetaData(CARO.list$asym17 , "asym17", col.name = "sample")

CARO.list$sym18 <- CreateSeuratObject(counts = data.sym18, project = "carotid",
                                      min.cells = 3, min.features = 200)
CARO.list$sym18 <- AddMetaData(CARO.list$sym18 , "sym", col.name = "stim")
CARO.list$sym18 <- AddMetaData(CARO.list$sym18 , "sym18", col.name = "sample")

CARO.list$sym19 <- CreateSeuratObject(counts = data.sym19, project = "carotid",
                                      min.cells = 3, min.features = 200)
CARO.list$sym19 <- AddMetaData(CARO.list$sym19 , "sym", col.name = "stim")
CARO.list$sym19 <- AddMetaData(CARO.list$sym19 , "sym19", col.name = "sample")

CARO.list$asym20 <- CreateSeuratObject(counts = data.asym20, project = "carotid",
                                       min.cells = 3, min.features = 200)
CARO.list$asym20 <- AddMetaData(CARO.list$asym20 , "asym", col.name = "stim")
CARO.list$asym20 <- AddMetaData(CARO.list$asym20 , "asym20", col.name = "sample")

CARO.list$sym21 <- CreateSeuratObject(counts = data.sym21, project = "carotid",
                                      min.cells = 3, min.features = 200)
CARO.list$sym21 <- AddMetaData(CARO.list$sym21 , "sym", col.name = "stim")
CARO.list$sym21 <- AddMetaData(CARO.list$sym21 , "sym21", col.name = "sample")

CARO.list$sym22 <- CreateSeuratObject(counts = data.sym22, project = "carotid",
                                      min.cells = 3, min.features = 200)
CARO.list$sym22 <- AddMetaData(CARO.list$sym22 , "sym", col.name = "stim")
CARO.list$sym22 <- AddMetaData(CARO.list$sym22 , "sym22", col.name = "sample")

CARO.list$sym23 <- CreateSeuratObject(counts = data.sym23, project = "carotid",
                                      min.cells = 3, min.features = 200)
CARO.list$sym23 <- AddMetaData(CARO.list$sym23 , "sym", col.name = "stim")
CARO.list$sym23 <- AddMetaData(CARO.list$sym23 , "sym23", col.name = "sample")

CARO.list$sym24 <- CreateSeuratObject(counts = data.sym24, project = "carotid",
                                      min.cells = 3, min.features = 200)
CARO.list$sym24 <- AddMetaData(CARO.list$sym24 , "sym", col.name = "stim")
CARO.list$sym24 <- AddMetaData(CARO.list$sym24 , "sym24", col.name = "sample")

########
CARO.list$sym7 <- CreateSeuratObject(counts = data.sym7, project = "carotid",
                                     min.cells = 0, min.features = 0)
CARO.list$sym7 <- AddMetaData(CARO.list$sym7 , "sym", col.name = "stim")
CARO.list$sym7 <- AddMetaData(CARO.list$sym7 , "sym7", col.name = "sample")

CARO.list$asym8 <- CreateSeuratObject(counts = data.asym8, project = "carotid",
                                      min.cells = 0, min.features = 0)
CARO.list$asym8 <- AddMetaData(CARO.list$asym8 , "asym", col.name = "stim")
CARO.list$asym8 <- AddMetaData(CARO.list$asym8 , "asym8", col.name = "sample")

CARO.list$sym9 <- CreateSeuratObject(counts = data.sym9, project = "carotid",
                                     min.cells = 0, min.features = 0)
CARO.list$sym9 <- AddMetaData(CARO.list$sym9 , "sym", col.name = "stim")
CARO.list$sym9 <- AddMetaData(CARO.list$sym9 , "sym9", col.name = "sample")

CARO.list$sym10 <- CreateSeuratObject(counts = data.sym10, project = "carotid",
                                      min.cells = 0, min.features = 0)
CARO.list$sym10 <- AddMetaData(CARO.list$sym10 , "sym", col.name = "stim")
CARO.list$sym10 <- AddMetaData(CARO.list$sym10 , "sym10", col.name = "sample")

CARO.list$sym11 <- CreateSeuratObject(counts = data.sym11, project = "carotid",
                                      min.cells = 0, min.features = 0)
CARO.list$sym11 <- AddMetaData(CARO.list$sym11 , "sym", col.name = "stim")
CARO.list$sym11 <- AddMetaData(CARO.list$sym11 , "sym11", col.name = "sample")

CARO.list$asym12 <- CreateSeuratObject(counts = data.asym12, project = "carotid",
                                       min.cells = 0, min.features = 0)
CARO.list$asym12 <- AddMetaData(CARO.list$asym12 , "asym", col.name = "stim")
CARO.list$asym12 <- AddMetaData(CARO.list$asym12 , "asym12", col.name = "sample")

CARO.list$asym13 <- CreateSeuratObject(counts = data.asym13, project = "carotid",
                                       min.cells = 0, min.features = 0)
CARO.list$asym13 <- AddMetaData(CARO.list$asym13 , "asym", col.name = "stim")
CARO.list$asym13 <- AddMetaData(CARO.list$asym13 , "asym13", col.name = "sample")

CARO.list$asym14 <- CreateSeuratObject(counts = data.asym14, project = "carotid",
                                       min.cells = 0, min.features = 0)
CARO.list$asym14 <- AddMetaData(CARO.list$asym14 , "asym", col.name = "stim")
CARO.list$asym14 <- AddMetaData(CARO.list$asym14 , "asym14", col.name = "sample")

CARO.list$sym15 <- CreateSeuratObject(counts = data.sym15, project = "carotid",
                                      min.cells = 0, min.features = 0)
CARO.list$sym15 <- AddMetaData(CARO.list$sym15 , "sym", col.name = "stim")
CARO.list$sym15 <- AddMetaData(CARO.list$sym15 , "sym15", col.name = "sample")

CARO.list$asym16 <- CreateSeuratObject(counts = data.asym16, project = "carotid",
                                       min.cells = 0, min.features = 0)
CARO.list$asym16 <- AddMetaData(CARO.list$asym16 , "asym", col.name = "stim")
CARO.list$asym16 <- AddMetaData(CARO.list$asym16 , "asym16", col.name = "sample")

CARO.list$asym17 <- CreateSeuratObject(counts = data.asym17, project = "carotid",
                                       min.cells = 0, min.features = 0)
CARO.list$asym17 <- AddMetaData(CARO.list$asym17 , "asym", col.name = "stim")
CARO.list$asym17 <- AddMetaData(CARO.list$asym17 , "asym17", col.name = "sample")

CARO.list$sym18 <- CreateSeuratObject(counts = data.sym18, project = "carotid",
                                      min.cells = 0, min.features = 0)
CARO.list$sym18 <- AddMetaData(CARO.list$sym18 , "sym", col.name = "stim")
CARO.list$sym18 <- AddMetaData(CARO.list$sym18 , "sym18", col.name = "sample")

CARO.list$sym19 <- CreateSeuratObject(counts = data.sym19, project = "carotid",
                                      min.cells = 0, min.features = 0)
CARO.list$sym19 <- AddMetaData(CARO.list$sym19 , "sym", col.name = "stim")
CARO.list$sym19 <- AddMetaData(CARO.list$sym19 , "sym19", col.name = "sample")

CARO.list$asym20 <- CreateSeuratObject(counts = data.asym20, project = "carotid",
                                       min.cells = 0, min.features = 0)
CARO.list$asym20 <- AddMetaData(CARO.list$asym20 , "asym", col.name = "stim")
CARO.list$asym20 <- AddMetaData(CARO.list$asym20 , "asym20", col.name = "sample")

CARO.list$sym21 <- CreateSeuratObject(counts = data.sym21, project = "carotid",
                                      min.cells = 0, min.features = 0)
CARO.list$sym21 <- AddMetaData(CARO.list$sym21 , "sym", col.name = "stim")
CARO.list$sym21 <- AddMetaData(CARO.list$sym21 , "sym21", col.name = "sample")

CARO.list$sym22 <- CreateSeuratObject(counts = data.sym22, project = "carotid",
                                      min.cells = 0, min.features = 0)
CARO.list$sym22 <- AddMetaData(CARO.list$sym22 , "sym", col.name = "stim")
CARO.list$sym22 <- AddMetaData(CARO.list$sym22 , "sym22", col.name = "sample")

CARO.list$sym23 <- CreateSeuratObject(counts = data.sym23, project = "carotid",
                                      min.cells = 0, min.features = 0)
CARO.list$sym23 <- AddMetaData(CARO.list$sym23 , "sym", col.name = "stim")
CARO.list$sym23 <- AddMetaData(CARO.list$sym23 , "sym23", col.name = "sample")

CARO.list$sym24 <- CreateSeuratObject(counts = data.sym24, project = "carotid",
                                      min.cells = 0, min.features = 0)
CARO.list$sym24 <- AddMetaData(CARO.list$sym24 , "sym", col.name = "stim")
CARO.list$sym24 <- AddMetaData(CARO.list$sym24 , "sym24", col.name = "sample")


# low-quality cell-----------------------------------------------------

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
CARO.list$sym7[["percent.mt"]] <- PercentageFeatureSet(CARO.list$sym7, pattern = "^MT-")
CARO.list$asym8[["percent.mt"]] <- PercentageFeatureSet(CARO.list$asym8, pattern = "^MT-")
CARO.list$sym9[["percent.mt"]] <- PercentageFeatureSet(CARO.list$sym9, pattern = "^MT-")
CARO.list$sym10[["percent.mt"]] <- PercentageFeatureSet(CARO.list$sym10, pattern = "^MT-")
CARO.list$sym11[["percent.mt"]] <- PercentageFeatureSet(CARO.list$sym11, pattern = "^MT-")
CARO.list$asym12[["percent.mt"]] <- PercentageFeatureSet(CARO.list$asym12, pattern = "^MT-")
CARO.list$asym13[["percent.mt"]] <- PercentageFeatureSet(CARO.list$asym13, pattern = "^MT-")
CARO.list$asym14[["percent.mt"]] <- PercentageFeatureSet(CARO.list$asym14, pattern = "^MT-")
CARO.list$sym15[["percent.mt"]] <- PercentageFeatureSet(CARO.list$sym15, pattern = "^MT-")
CARO.list$asym16[["percent.mt"]] <- PercentageFeatureSet(CARO.list$asym16, pattern = "^MT-")
CARO.list$asym17[["percent.mt"]] <- PercentageFeatureSet(CARO.list$asym17, pattern = "^MT-")
CARO.list$sym18[["percent.mt"]] <- PercentageFeatureSet(CARO.list$sym18, pattern = "^MT-")
CARO.list$sym19[["percent.mt"]] <- PercentageFeatureSet(CARO.list$sym19, pattern = "^MT-")
CARO.list$asym20[["percent.mt"]] <- PercentageFeatureSet(CARO.list$asym20, pattern = "^MT-")
CARO.list$sym21[["percent.mt"]] <- PercentageFeatureSet(CARO.list$sym21, pattern = "^MT-")
CARO.list$sym22[["percent.mt"]] <- PercentageFeatureSet(CARO.list$sym22, pattern = "^MT-")
CARO.list$sym23[["percent.mt"]] <- PercentageFeatureSet(CARO.list$sym23, pattern = "^MT-")
CARO.list$sym24[["percent.mt"]] <- PercentageFeatureSet(CARO.list$sym24, pattern = "^MT-")

VlnPlot(CARO.list$sym7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CARO.list$asym8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CARO.list$sym9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CARO.list$sym10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CARO.list$sym11, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CARO.list$asym12, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CARO.list$asym13, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CARO.list$asym14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CARO.list$sym15, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CARO.list$asym16, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CARO.list$asym17, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CARO.list$sym18, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CARO.list$sym19, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CARO.list$asym20, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CARO.list$sym21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CARO.list$sym22, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CARO.list$sym23, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CARO.list$sym24, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# low-quality cell-----------------------------------------------------
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(CARO.list$sym1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CARO.list$sym1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(CARO.list$asym2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CARO.list$asym2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(CARO.list$asym3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CARO.list$asym3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(CARO.list$asym3A, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CARO.list$asym3A, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(CARO.list$asym4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CARO.list$asym4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(CARO.list$asym5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CARO.list$asym5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(CARO.list$sym6, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CARO.list$sym6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# Filtering low-quality cell-----------------------------------------------
CARO.list$sym7 <- subset(CARO.list$sym7, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)
CARO.list$asym8 <- subset(CARO.list$asym8, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)
CARO.list$sym9 <- subset(CARO.list$sym9, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)
CARO.list$sym10 <- subset(CARO.list$sym10, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)
CARO.list$sym11 <- subset(CARO.list$sym11, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)
CARO.list$asym12 <- subset(CARO.list$asym12, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)
CARO.list$asym13 <- subset(CARO.list$asym13, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)
CARO.list$asym14 <- subset(CARO.list$asym14, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)
CARO.list$sym15 <- subset(CARO.list$sym15, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)
CARO.list$asym16 <- subset(CARO.list$asym16, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)
CARO.list$asym17 <- subset(CARO.list$asym17, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)
CARO.list$sym18 <- subset(CARO.list$sym18, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)
CARO.list$sym19 <- subset(CARO.list$sym19, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)
CARO.list$asym20 <- subset(CARO.list$asym20, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)
CARO.list$sym21 <- subset(CARO.list$sym21, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)
CARO.list$sym22 <- subset(CARO.list$sym22, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)
CARO.list$sym23 <- subset(CARO.list$sym23, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)
CARO.list$sym24 <- subset(CARO.list$sym24, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30)


plot1 <- FeatureScatter(CARO.list$sym1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CARO.list$sym1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(CARO.list$asym2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CARO.list$asym2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(CARO.list$asym3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CARO.list$asym3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(CARO.list$asym3A, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CARO.list$asym3A, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(CARO.list$asym4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CARO.list$asym4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(CARO.list$asym5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CARO.list$asym5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(CARO.list$sym6, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CARO.list$sym6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# UMAP------------------------------------------------------------------------------

# normalize and identify variable features for each dataset independently
all.genes <- rownames(CARO.list["asym"]) +  rownames(CARO.list["sym"])

all.genes <- append(rownames(CARO.list$sym7),rownames(CARO.list$asym8))
all.genes <- append(all.genes,rownames(CARO.list$sym9))
all.genes <- append(all.genes,rownames(CARO.list$sym10))
all.genes <- append(all.genes,rownames(CARO.list$sym11))
all.genes <- append(all.genes,rownames(CARO.list$asym12))
all.genes <- append(all.genes,rownames(CARO.list$asym13))
all.genes <- append(all.genes,rownames(CARO.list$asym14))
all.genes <- append(all.genes,rownames(CARO.list$sym15))
all.genes <- append(all.genes,rownames(CARO.list$asym16))
all.genes <- append(all.genes,rownames(CARO.list$asym17))
all.genes <- append(all.genes,rownames(CARO.list$sym18))
all.genes <- append(all.genes,rownames(CARO.list$sym19))
all.genes <- append(all.genes,rownames(CARO.list$asym20))
all.genes <- append(all.genes,rownames(CARO.list$sym21))
all.genes <- append(all.genes,rownames(CARO.list$sym22))
all.genes <- append(all.genes,rownames(CARO.list$sym23))
all.genes <- append(all.genes,rownames(CARO.list$sym24))



CARO.list <- lapply(X = CARO.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- ScaleData(x, features = all.genes)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- RunPCA(x, features = VariableFeatures(object = x), verbose = FALSE)
})


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = CARO.list)

# Integration
CARO.anchors <- FindIntegrationAnchors(object.list = CARO.list, anchor.features = features, reduction = "rpca",k.filter = 35)

# this command creates an 'integrated' data assay
CARO.combined <- IntegrateData(anchorset =CARO.anchors, k.weight = 35)

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(CARO.combined) <- "integrated"

CARO.combined <- ScaleData(CARO.combined, verbose = FALSE)
CARO.combined <- RunPCA(CARO.combined, npcs = 50, verbose = FALSE)
ElbowPlot(CARO.combined, ndims = 50)
CARO.combined <- RunUMAP(CARO.combined, reduction = "pca", dims = 1:40)
CARO.combined <- FindNeighbors(CARO.combined, reduction = "pca", dims = 1:40)
CARO.combined <- FindClusters(CARO.combined, resolution = 0.3)

DimPlot(CARO.combined, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(CARO.combined, reduction = "umap", split.by = "stim2")
DimPlot(CARO.combined, reduction = "umap", group.by = "stim")


# Visualization ------
DimPlot(CARO.combined, reduction = "umap", split.by = "stim")
DimPlot(CARO.combined, reduction = "umap", repel = TRUE)

saveRDS(CARO.combined,"~/Desktop/Human_carotid_public/RDS/Human_carotid_sample2_all_240510.rds")
CARO.combined1 <- readRDS("~/Desktop/Human_carotid_public/RDS/Human_carotid_sample2_all_240510.rds")

CARO.combined1 <- RunUMAP(CARO.combined1, reduction = "pca", dims = 1:40)
CARO.combined1 <- FindNeighbors(CARO.combined1, reduction = "pca", dims = 1:40)
CARO.combined1 <- FindClusters(CARO.combined1, resolution = 0.3)

DimPlot(CARO.combined1, reduction = "umap", label = TRUE, repel = TRUE)

# cluster change name----
library(ggplot2)

levels(CARO.combined) # [1] "0" "1" "2" "3" "4" "5" "6" "7 "8" "9"
levels(CARO.combined) <- c("0", "1", "5", "7", "2", "3", "4", "9", "6", "8", "10", "11")

new.cluster.ids <- c("0","1","2","3","4", "5", "6","7","8","9","10","11")
names(new.cluster.ids) <- levels(CARO.combined)
CARO.combined<- RenameIdents(CARO.combined, new.cluster.ids)
DimPlot(CARO.combined, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

DefaultAssay(CARO.combined) <- "RNA"


VlnPlot(CARO.combined, assay = "RNA",features = c("MMP19") ,pt.size = 0)
FeaturePlot(CARO.combined, features=c('MMP19'), min.cutoff=0, max.cutoff='q90')

# Myeloid----
VlnPlot(CARO.combined, assay = "RNA",features = c("CD68","CD14","CSF3R","CSF1R","IL1B","S100A8", "S100A9","MKI67","FCGR3A") ,pt.size = 0)
VlnPlot(CARO.combined, assay = "RNA",features = c("CD163","CD14","PTPRC") ,pt.size = 0)

# Bcells-----
VlnPlot(CARO.combined1, features = c("CD79A","CD79B","FCER2","CD22","MYH11","CD34") ,pt.size = 0)

VlnPlot(CARO.combined1, features = c("CDH5","ACTA2","VWF","MKI67","MYH11","CD34","IGHG4") ,pt.size = 0)

#SMC
VlnPlot(CARO.combined, assay = "RNA", features = c("HLA-DQA1","HLA-DQB1","HLA-DRB5","HLA-DPB1") ,pt.size = 0)
VlnPlot(CARO.combined, assay = "RNA", features = c("ACTA2","HAND2","ERBB4","SERTAD4","DLX5","SPARCL1","TAGLN","MYOCD") ,pt.size = 0)


# Tcells-----
VlnPlot(CARO.combined1, features = c("CD3E","CD4","CD8A","IL7R","LEF1","GZMK") ,pt.size = 0)

VlnPlot(CARO.combined, features = c("GLS"), split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))

# Heatmap------------------------------------------------------------------------------
CARO.combined.markers <- FindAllMarkers(CARO.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CARO.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- CARO.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

CARO.combined01234 <- subset(CARO.combined, idents = c(0,1,2,3,4))
CARO.combined56789 <- subset(CARO.combined, idents = c(5,6,7,8,9))
CARO.combined1011121314 <- subset(CARO.combined, idents = c(10,11,12,13,14))
CARO.combined1516171819 <- subset(CARO.combined, idents = c(15,16,17,18,19))
CARO.combined202122 <- subset(CARO.combined, idents = c(20,21,22))

DoHeatmap(CARO.combined, features = top10$gene) + NoLegend()
DoHeatmap(CARO.combined01234, features = top10$gene) + NoLegend()
DoHeatmap(CARO.combined56789, features = top10$gene) + NoLegend()
DoHeatmap(CARO.combined1011121314, features = top10$gene) + NoLegend()
DoHeatmap(CARO.combined1516171819, features = top10$gene) + NoLegend()
DoHeatmap(CARO.combined202122, features = top10$gene) + NoLegend()

write.table(top10, "table.txt", quote=F, col.names=F, append=T)

#Vinplot----
VlnPlot(CARO.combined, features = c("SLC25A44"), split.by = "stim" )
VlnPlot(CARO.combined, features = c("FCGR3A") )
VlnPlot(CARO.combined, features = c("Th17A") )
VlnPlot(CARO.combined, features = c("CD68") , levels = new_order)

#myeloid-------
DefaultAssay(CARO.combined_Mye) <- "integrated"
CARO.combined_Mye <- subset(CARO.combined, idents = c(0, 6, 7, 8, 9, 12, 18))
CARO.combined_Mye <- ScaleData(CARO.combined_Mye, verbose = FALSE)
CARO.combined_Mye <- RunPCA(CARO.combined_Mye, npcs = 50, verbose = FALSE)
CARO.combined_Mye <- RunUMAP(CARO.combined_Mye , reduction = "pca", dims = 1:10)
CARO.combined_Mye <- FindNeighbors(CARO.combined_Mye , reduction = "pca", dims = 1:10)
CARO.combined_Mye <- FindClusters(CARO.combined_Mye, resolution = 0.8)

DimPlot(CARO.combined_Mye, reduction = "umap", label = TRUE, repel = TRUE) 
DimPlot(CARO.combined_Mye, reduction = "umap", split.by = "stim2")


CARO.combined_Mye.markers <- FindAllMarkers(CARO.combined_Mye, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CARO.combined_Mye.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- CARO.combined_Mye.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(CARO.combined_Mye, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)

saveRDS(CARO.combined_Mye,"~/Desktop/Human_carotid_public/RDS/Human_carotid_Mye_240426.rds")
CARO.combined_Mye <- readRDS("~/Desktop/Human_carotid_public/RDS/Human_carotid_Mye_240426.rds")
CARO.combined_Mye2 <- readRDS("~/Desktop/Human_carotid_public/RDS/Human_carotid_Mye_240426.rds")

DimPlot(CARO.combined_Mye, reduction = "umap", label = TRUE, repel = TRUE) 

DefaultAssay(CARO.combined_Mye) <- "RNA"

VlnPlot(CARO.combined_Mye, assay = "RNA",features = c("CD68") ,pt.size = 0)
VlnPlot(CARO.combined_Mye, assay = "RNA",features = c("TNF") ,pt.size = 0)
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("S100A8","FCGR3A","IL1B","LYZ","VCAN","FCN1") ,pt.size = 0)
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("IL1B","S100A12","NRG1","EREG","S100A8") ,pt.size = 0)
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("CEP295NL","PLCG2","NEB","CATSPERG") ,pt.size = 0)
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("CCL18","HS3ST2","APOE","RARRES1") ,pt.size = 0)
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("SPP1","SPOCD1","LPL","PAQR5","TREM2") ,pt.size = 0)
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("CLNK","TACSTD2","IDO1","SLCO5A1") ,pt.size = 0)
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("GZMB","PACSIN1","EPHB1","JCHAIN","IFNA") ,pt.size = 0)

VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("MAFB","STAT1","STAT5B","TCF7L2","RELA") ,pt.size = 0)
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("MLXIPL","MYCN","NUPR1","BHLHE40","ECSIT") ,pt.size = 0)

#apoptotic marker
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("CCNE2", "CDK1", "PCNA","CASP3", "CASP9", "BAX") ,pt.size = 0)

#tcell
VlnPlot(CARO.combined_Mye, assay = "RNA",features = c("CD3E","CD4","CD8A","IL7R","LEF1","GZMK") ,pt.size = 0)

#c1q-hi
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("C1QA","C1QB","SELENOP","GIPC2","C1QC","CD68","FOLR2","PLTP") ,pt.size = 0)
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("FOLR2","LYVE1") ,pt.size = 0)

#proliferative
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("MIF","FTL","TUBA1B","LGALS1","LGALS3") ,pt.size = 0)
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("SPC25","MKI67","TOP2A","BIRC5","FABP5") ,pt.size = 0)

#DC
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("HLA-DQA1","HLA-DQB1","HLA-DRB5","HLA-DPB1") ,pt.size = 0)
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("CLEC10A","CLEC9A") ,pt.size = 0)#CLEC10Aがcdc1マーカー、9Aがcdc2マーカー

#ACTA2 
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("ACTA2","HAND2","ERBB4","SERTAD4","DLX5","SPARCL1","TAGLN","MYOCD") ,pt.size = 0)

#SMC
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("HLA-DQA1","HLA-DQB1","HLA-DRB5","HLA-DPB1") ,pt.size = 0)

FeaturePlot(CARO.combined_Mye, features=c('CD3E'), min.cutoff=0, max.cutoff='q90')


#15-19-20-Tcell
#myeloid-------
DefaultAssay(CARO.combined_Mye) <- "integrated"
CARO.combined_Mye2 <- subset(CARO.combined_Mye, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 21, 22))
CARO.combined_Mye2 <- ScaleData(CARO.combined_Mye2, verbose = FALSE)
CARO.combined_Mye2 <- FindVariableFeatures(CARO.combined_Mye2, selection.method = "vst", nfeatures = 2000)
CARO.combined_Mye2 <- RunPCA(CARO.combined_Mye2, npcs = 50, verbose = FALSE)
CARO.combined_Mye2 <- RunUMAP(CARO.combined_Mye2 , reduction = "pca", dims = 1:30)#9
CARO.combined_Mye2 <- FindNeighbors(CARO.combined_Mye2 , reduction = "pca", dims = 1:30)
CARO.combined_Mye2 <- FindClusters(CARO.combined_Mye2, resolution = 0.35)#0.35

DimPlot(CARO.combined_Mye2, reduction = "umap", label = TRUE, repel = TRUE) 
DimPlot(CARO.combined_Mye2, reduction = "umap", split.by = "stim")
DimPlot(CARO.combined_Mye2, reduction = "umap", split.by = "sample")
DimPlot(CARO.combined_Mye2, reduction = "umap", split.by = "stim", label = TRUE, repel = TRUE, label.size = 8) 

CARO.combined_Mye2_sym <- subset(CARO.combined_Mye2, subset = stim == "sym")
CARO.combined_Mye2_asym <- subset(CARO.combined_Mye2, subset = stim == "asym")

DimPlot(CARO.combined_Mye2_sym, reduction = "umap", split.by = "sample")
DimPlot(CARO.combined_Mye2_asym, reduction = "umap", split.by = "sample")

saveRDS(CARO.combined_Mye2,"~/Desktop/Human_carotid_public/RDS/Human_carotid_sample2_Mye_240824.rds")
CARO.combined_Mye2 <- readRDS("~/Desktop/Human_carotid_public/RDS/Human_carotid_sample2_Mye_240824.rds")
DefaultAssay(CARO.combined_Mye2) <- "RNA"

###random sampling
cluster_0 <- subset(CARO.combined_Mye2, idents = c("0"))
cluster_1 <- subset(CARO.combined_Mye2, idents = c("1"))
cluster_2 <- subset(CARO.combined_Mye2, idents = c("2"))
cluster_3 <- subset(CARO.combined_Mye2, idents = c("3"))
cluster_4 <- subset(CARO.combined_Mye2, idents = c("4"))
cluster_5 <- subset(CARO.combined_Mye2, idents = c("5"))
cluster_6 <- subset(CARO.combined_Mye2, idents = c("6"))
cluster_7 <- subset(CARO.combined_Mye2, idents = c("7"))
cluster_8 <- subset(CARO.combined_Mye2, idents = c("8"))
cluster_9 <- subset(CARO.combined_Mye2, idents = c("9"))
cluster_10 <- subset(CARO.combined_Mye2, idents = c("10"))
cluster_11 <- subset(CARO.combined_Mye2, idents = c("11"))

object_cluster_0 <- subset(CARO.combined_Mye2, cells = Cells(cluster_0), downsample = (10000*ncol(cluster_0)/ncol(CARO.combined_Mye2)))
object_cluster_1 <- subset(CARO.combined_Mye2, cells = Cells(cluster_1), downsample = (10000*ncol(cluster_1)/ncol(CARO.combined_Mye2)))
object_cluster_2 <- subset(CARO.combined_Mye2, cells = Cells(cluster_2), downsample = (10000*ncol(cluster_2)/ncol(CARO.combined_Mye2)))
object_cluster_3 <- subset(CARO.combined_Mye2, cells = Cells(cluster_3), downsample = (10000*ncol(cluster_3)/ncol(CARO.combined_Mye2)))
object_cluster_4 <- subset(CARO.combined_Mye2, cells = Cells(cluster_4), downsample = (10000*ncol(cluster_4)/ncol(CARO.combined_Mye2)))
object_cluster_5 <- subset(CARO.combined_Mye2, cells = Cells(cluster_5), downsample = (10000*ncol(cluster_5)/ncol(CARO.combined_Mye2)))
object_cluster_6 <- subset(CARO.combined_Mye2, cells = Cells(cluster_6), downsample = (10000*ncol(cluster_6)/ncol(CARO.combined_Mye2)))
object_cluster_7 <- subset(CARO.combined_Mye2, cells = Cells(cluster_7), downsample = (10000*ncol(cluster_7)/ncol(CARO.combined_Mye2)))
object_cluster_8 <- subset(CARO.combined_Mye2, cells = Cells(cluster_8), downsample = (10000*ncol(cluster_8)/ncol(CARO.combined_Mye2)))
object_cluster_9 <- subset(CARO.combined_Mye2, cells = Cells(cluster_9), downsample = (10000*ncol(cluster_9)/ncol(CARO.combined_Mye2)))
object_cluster_10 <- subset(CARO.combined_Mye2, cells = Cells(cluster_10), downsample = (10000*ncol(cluster_10)/ncol(CARO.combined_Mye2)))
object_cluster_11 <- subset(CARO.combined_Mye2, cells = Cells(cluster_11), downsample = (10000*ncol(cluster_11)/ncol(CARO.combined_Mye2)))

object_subset_2 <- merge(object_cluster_0, y = c(object_cluster_1, object_cluster_2,object_cluster_3,object_cluster_4,object_cluster_5,object_cluster_6,object_cluster_7,object_cluster_8,object_cluster_9,object_cluster_10,object_cluster_11))

object_subset <- subset(CARO.combined_Mye2, cells = Cells(object_subset_2))

DimPlot(object_subset, reduction = "umap", label = TRUE, repel = TRUE) 

DefaultAssay(object_subset) <- "integrated"
object_subset <- subset(object_subset, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 21, 22))
object_subset <- ScaleData(object_subset, verbose = FALSE)
object_subset <- FindVariableFeatures(object_subset, selection.method = "vst", nfeatures = 2000)
object_subset <- RunPCA(object_subset, npcs = 50, verbose = FALSE)
object_subset <- RunUMAP(object_subset , reduction = "pca", dims = 1:30)#9
object_subset <- FindNeighbors(object_subset , reduction = "pca", dims = 1:30)
object_subset <- FindClusters(object_subset, resolution = 0.35)#0.35

DimPlot(CARO.combined_Mye2, reduction = "umap", label = TRUE, repel = TRUE) 


install.packages("ggpie")
library(ggpie)
Idents(CARO.combined_Mye2_sym)
CARO.combined_Mye2_sym$seurat_clusters <- Idents(CARO.combined_Mye2_sym)

df <- CARO.combined_Mye2_sym@meta.data
ggpie(data = df, 
      group_key = "seurat_clusters",
      count_type = "full",
      label_type = "none")


Idents(CARO.combined_Mye2_asym)
CARO.combined_Mye2_asym$seurat_clusters <- Idents(CARO.combined_Mye2_asym)

df <- CARO.combined_Mye2_asym@meta.data
ggpie(data = df, 
      group_key = "seurat_clusters",
      count_type = "full",
      label_type = "none")

DimPlot(CARO.combined_Mye2, reduction = "umap", label = TRUE, repel = TRUE) 
DimPlot(CARO.combined_Mye2_sym, reduction = "umap", label = TRUE, repel = TRUE) 
DimPlot(CARO.combined_Mye2_asym, reduction = "umap", label = TRUE, repel = TRUE) 


#Dotplot
CARO.combined_Mye2_trans <- SCTransform (CARO.combined_Mye2)

cd_genes <- c("LYVE1","FOLR2","IGF1","C1QA","C1QB","TNF","CCL18","APOE","RARRES1","TREM2","SPP1","LPL","PLCG2","S100A12","S100A8","VCAN","IL1B","EREG","AREG","CXCL2","CXCL3","HLA-DQA1","HLA-DQB1","CLEC10A","CD1C","CSF3R","FCGR3B","NAMPT","SPC25","MKI67","TOP2A","ACTA2","SPARCL1","TAGLN")
DotPlot(CARO.combined_Mye2_trans,features = cd_genes)+RotatedAxis()+coord_flip()
#final
cd_genes <- c("LYVE1","FOLR2","IGF1",　"TREM2","SPP1","LPL",　"C1QA","C1QB","TNF",　"C1QC",　"PLCG2", "CCL18","APOE","RARRES1","HLA-DQA1","HLA-DQB1","CLEC10A","CD1C","IL1B","EREG","AREG","CXCL2","CXCL3","S100A12","S100A8","VCAN","CSF3R","FCGR3B","NAMPT","SPC25","MKI67","TOP2A","ACTA2","SPARCL1","TAGLN")
DotPlot(CARO.combined_Mye2_trans,features = cd_genes)+RotatedAxis()+coord_flip()


VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("CD1C","CCL4","RAB27A","MARCKSL1","NFIL3","EHD1","CD55") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("CD52","CD68","CD14") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("ISG15","IFIT2","ISG20","IFIT1","CXCL8","CCL8","CCL2","TNFSF10","IL1RN","STAT1","IRF7","TNFSF13B") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("S100A8","FCGR3A","IL1B","LYZ","VCAN","FCN1") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("IL1B","S100A12","NRG1","EREG","S100A8") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("CEP295NL","PLCG2","NEB","CATSPERG") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("CCL18","HS3ST2","APOE","RARRES1") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("SPP1","SPOCD1","LPL","PAQR5","TREM2") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("CLNK","TACSTD2","IDO1","SLCO5A1") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("GZMB","PACSIN1","EPHB1","JCHAIN","IFNA") ,pt.size = 0)

VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("MAFB","STAT1","STAT5B","TCF7L2","RELA") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("MLXIPL","MYCN","NUPR1","BHLHE40","ECSIT") ,pt.size = 0)

#apoptotic marker
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("CCNE2", "CDK1", "PCNA","CASP3", "CASP9", "BAX") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("PLCG2") ,pt.size = 0)

#tcell
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("CD3E","CD4","CD8A","IL7R","LEF1","GZMK") ,pt.size = 0)

#bcell
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("CD19","FCRLA","FCRL2","PAX5") ,pt.size = 0)

#plasma cell
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("CD138","CD78") ,pt.size = 0)

#NK
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("NKG7") ,pt.size = 0)

#IL1B+
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("IL1B","EREG","AREG","HBEGF","CXCL2","VEGFA","CD44","IL1RN","CXCL3") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("IL1B","EREG","AREG","HBEGF","CXCL2","CXCL3") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("TNF","CXCL2") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("IL1B","AREG") ,pt.size = 0)

#TREM2+
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("TREM2","SPP1","LPL","APOE","ABCG1") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("TREM2","SPP1") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("APOE","C1QC") ,pt.size = 0)

#LYVE1+
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("LYVE1","FOLR2","IGF1","C1QA") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("LYVE1","FOLR2") ,pt.size = 0)

#mast 
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("HDC","KIT") ,pt.size = 0)

#Mreg DC　
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("FSCN1","CCR7","HLA-DRA","HLA-DPA1") ,pt.size = 0)

#clasical mono
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("LYZ","S100A12","S100A8","CD14","VCAN","FOS") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("S100A12","VCAN") ,pt.size = 0)

#nonclasical mono
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("CSF1R","FCGR3A","CXCL10","CCL4","CXCL8") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("CSF1R","FCGR3A") ,pt.size = 0)

#Neutrophils
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("CSF3R","FCGR3B","S100A8","CXCR2","NAMPT") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("FCGR3B","CSF3R","CXCR2","SLC27A2") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA",features = c("FCGR3B","CSF3R") ,pt.size = 0)

#c1q-hi
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("C1QA","C1QB","SELENOP","GIPC2","C1QC","CD68","FOLR2","PLTP") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("FOLR2","LYVE1","TREM2","SPP1") ,pt.size = 0)

#C1Q+ resident like
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("FTL","CD163","HMOX1","MAF","MRC1","SIGLEC1") ,pt.size = 0) #CD169=SIGLEC1,CD206=MRC1
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("APOC1","APOE","ALOX5AP","INHBA","GPNMB") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("C3AR1","CFD","FTH1","ACP5") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("FTL","PLCG2") ,pt.size = 0) #CD169=SIGLEC1,CD206=MRC1


#proliferative
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("MIF","FTL","TUBA1B","LGALS1","LGALS3") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("SPC25","MKI67","TOP2A","BIRC5","FABP5") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("MKI67","TOP2A") ,pt.size = 0)

#DC
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("HLA-DQA1","HLA-DQB1","HLA-DRB5","HLA-DPB1","CLEC10A","CLEC9A","CD1E","CD1C","CD74") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("CLEC10A","CLEC9A") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("GZMB","PACSIN1","EPHB1","JCHAIN") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("HLA-DQA1","HLA-DQB1","HLA-DPB1","CLEC10A","CD1E","CD1C") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("CD1C","CLEC10A") ,pt.size = 0)

#ACTA2 
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("ACTA2","HAND2","ERBB4","SERTAD4","DLX5","SPARCL1","TAGLN","MYOCD") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("ACTA2","SPARCL1","TAGLN","MYOCD") ,pt.size = 0)
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("ACTA2","TAGLN") ,pt.size = 0)

#SMC
VlnPlot(CARO.combined_Mye2, assay = "RNA", features = c("HLA-DQA1","HLA-DQB1","HLA-DRB5","HLA-DPB1") ,pt.size = 0)

CARO.combined_Mye2.markers <- FindAllMarkers(CARO.combined_Mye2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CARO.combined_Mye2.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- CARO.combined_Mye2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 <- CARO.combined_Mye2.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

DoHeatmap(CARO.combined_Mye2, features = top5$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


new.cluster.ids <- c("0","1","2","3","4","5","6","7")
new.cluster.ids <- c("0","1","2","3","4","5","6","2")

names(new.cluster.ids) <- levels(CARO.combined_Mye2)
CARO.combined_Mye2<- RenameIdents(CARO.combined_Mye2, new.cluster.ids)
DimPlot(CARO.combined_Mye2, reduction = "umap", label = TRUE, repel = TRUE) 

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = CARO.combined_1)))
names(x = ident.colors) <- levels(x = CARO.combined_1)
cell.colors <- ident.colors[Idents(object = CARO.combined_1)]
names(x = cell.colors) <- colnames(x = CARO.combined_1)
cell.colors <- c("#F8766D","#ABA300","#0CB702","#00A9FF","00C19A","FF61CC")

CARO.combined_1.markers <- FindAllMarkers(CARO.combined_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CARO.combined_1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- CARO.combined_1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(CARO.combined_1, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)

FeaturePlot(CARO.combined_Mye2, features=c('ASMT'), min.cutoff=0, max.cutoff='q90')
FeaturePlot(CARO.combined_Mye2, features=c('APOE'), min.cutoff=0, max.cutoff='q90')

FeaturePlot(CARO.combined_Mye2, features=c('ISG15'), min.cutoff=0, max.cutoff='q90')
FeaturePlot(CARO.combined_Mye2, features=c('TREM2'), min.cutoff=1.0, max.cutoff='q90')
FeaturePlot(CARO.combined_Mye2, features=c('TNF'), min.cutoff=0, max.cutoff='q90')#炎症マーカー
FeaturePlot(CARO.combined_Mye2, features=c('IL1B'), min.cutoff=0, max.cutoff='q90')#炎症マーカー
FeaturePlot(CARO.combined_Mye2, features=c('CCL4'), min.cutoff=0, max.cutoff='q90')#炎症マーカー
FeaturePlot(CARO.combined_Mye2, features=c('FOLR2'), min.cutoff=1.0, max.cutoff='q90')
FeaturePlot(CARO.combined_Mye2, features=c('SPP1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(CARO.combined_Mye2, features=c('ISG15'), min.cutoff=1.9, max.cutoff='q90')
FeaturePlot(CARO.combined_Mye2 ,features=c('FCGR3A'), min.cutoff=0, max.cutoff='q90')
FeaturePlot(CARO.combined_Mye2 ,features=c('MALAT1'), min.cutoff=6, max.cutoff='q90')
FeaturePlot(CARO.combined_Mye2 ,features=c('NEAT1'), min.cutoff=0, max.cutoff='q90')
FeaturePlot(CARO.combined_Mye2 ,features="percent.mt", min.cutoff=0, max.cutoff='q90')


DefaultAssay(CARO.combined_Mye2) <- "RNA"

FeaturePlot(CARO.combined_Mye2, features=c('MMP14'), min.cutoff=0, max.cutoff='q90')
2,7,9,12,14,19,25

CARO.combined_Mye2_sym <- subset(CARO.combined_Mye2, subset = stim == "sym")
DefaultAssay(CARO.combined_Mye2_sym) <- "RNA"
FeaturePlot(CARO.combined_Mye2_sym, features=c('MMP14'), min.cutoff=0, max.cutoff='q90')& 
  scale_colour_gradientn(colours = brewer.pal(n = 10, name = "YlOrRd"), 
                         limits = c(0,1.65))
2,7,9,12,14,19,25

CARO.combined_Mye2_asym <- subset(CARO.combined_Mye2, subset = stim == "asym")
DefaultAssay(CARO.combined_Mye2_asym) <- "RNA"
FeaturePlot(CARO.combined_Mye2_asym, features=c('MMP14'), min.cutoff=0, max.cutoff='q90')& 
  scale_colour_gradientn(colours = brewer.pal(n = 10, name = "YlOrRd"), 
                         limits = c(0,1.65))
2,7,9,12,14,19,25

CARO.combined_Mye2 <- AddMetaData(CARO.combined_Mye2 , CARO.combined_Mye2@active.ident, col.name = "seurat_clusters")

saveRDS(CARO.combined_Mye2,"~/Desktop/Human_carotid_public/RDS/Human_carotid_sample2_Mye_240824.rds")
CARO.combined_Mye2 <- readRDS("~/Desktop/Human_carotid_public/RDS/Human_carotid_sample2_Mye_240824.rds")






#Macrophages

DefaultAssay(CARO.combined_Mye2) <- "integrated"
CARO.combined_Mac <- subset(CARO.combined_Mye2, idents = c(0, 1, 2, 3, 4, 5, 7, 10, 11))
CARO.combined_Mac <- ScaleData(CARO.combined_Mac, verbose = FALSE)
CARO.combined_Mac <- FindVariableFeatures(CARO.combined_Mac, selection.method = "vst", nfeatures = 2000)
CARO.combined_Mac <- RunPCA(CARO.combined_Mac, npcs = 50, verbose = FALSE)
CARO.combined_Mac <- RunUMAP(CARO.combined_Mac , reduction = "pca", dims = 1:19)
CARO.combined_Mac <- FindNeighbors(CARO.combined_Mac , reduction = "pca", dims = 1:19)
ElbowPlot(CARO.combined_Mac)

DimPlot(CARO.combined_Mac, reduction = "umap", label = TRUE, repel = TRUE) 
DimPlot(CARO.combined_Mac, reduction = "umap", split.by = "stim")
DimPlot(CARO.combined_Mac, reduction = "umap", split.by = "sample")

CARO.combined_Mac_sym <- subset(CARO.combined_Mac, subset = stim == "sym")
CARO.combined_Mac_asym <- subset(CARO.combined_Mac, subset = stim == "asym")

DimPlot(CARO.combined_Mac_sym, reduction = "umap", split.by = "sample")
DimPlot(CARO.combined_Mac_asym, reduction = "umap", split.by = "sample")

saveRDS(CARO.combined_Mac,"~/Desktop/Human_carotid_public/RDS/Human_carotid_sample2_Mac_240918.rds")
CARO.combined_Mac <- readRDS("~/Desktop/Human_carotid_public/RDS/Human_carotid_sample2_Mac_240918.rds")




VlnPlot(CARO.combined_1, features = c("FCGR3A","FCGR2A","ITGB2","CD16","CEACAM8","CD55","CD44","ELANE") )
VlnPlot(CARO.combined_1, features = c("PTPRC","CD14","MPO","LYZ","ITGB2","CD16","FCGR3A","FCGR3B","CD55","CD44","ELANE","VCAN", "CD52","FCN1","S100A12","IL1B","CSF1R","S100A8","S100A9","CAMP","DEFB4A","CTSS","ITGAM","CD68","CYBB") )
VlnPlot(CARO.combined_1, features = c("DEFB1","DEFA1","DEFB103B") )
VlnPlot(CARO.combined_1, features = c("S100A6","FOS","S100A12","S100A9","HMGB2","KLF6"))


VlnPlot(CARO.combined_1, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(CARO.combined_1, features = c("TLR4","TLR2","TLR7","TLR3","TLR12", "TLR9") )
VlnPlot(CARO.combined_1, features = c("HDC","KIT","TPSAB1","CCL5") )
VlnPlot(CARO.combined_1, features = c("IL1B","TNF","CXCL8","NR1H3", "CEBPA","CEBPB") )
VlnPlot(CARO.combined_1, features = c("CD68","CSF1R","ITGAM","S100A8","S100A9","S100A12") )
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("S100A8","FCGR3A","IL1B","LYZ","VCAN","FCN1") ,pt.size = 0)
VlnPlot(CARO.combined_1, features = c("CASP1","CASP4","DEFB4A","CTSS","CYBB","IL1B","VEGFA") )
VlnPlot(CARO.combined_1, features = c("CCL2","CCL20","CCL4","HMGB2","KLF6","VEGFA") )
VlnPlot(CARO.combined_Mye, assay = "RNA", features = c("TREM2","C1QA1","CLEC10A","SPP1") ,pt.size = 0)


VlnPlot(CARO.combined_1, features = c("IL1B","TNF","CXCL8","CXCL3", "ABCG1","ABCA1") )
VlnPlot(CARO.combined_1, features = c("CD9","TREM2","MRC1","FCER2","C1QA", "C1QB") )
VlnPlot(CARO.combined_1, features = c("C1QC","CD59","APOE","APOC1","NUPR1", "C1QB") )

VlnPlot(CARO.combined_1, features = c("CEACAM8","MPO","IL6","CD63","S100A8","S100A9") )
VlnPlot(CARO.combined_1, features = c("MKI67","ASP175","CSF1R","IL1B","S100A8","S100A9") )
VlnPlot(CARO.combined_1, features = c("CD80","CD86","CD40","CD209","IL10","CD83") )
VlnPlot(CARO.combined_1, features = c("MT2A","MALAT1","HNRNPU","CD3E","KLF4","CCL5") )
VlnPlot(CARO.combined_1, features = c("CD80","CD86","CD40","CD209","IL10","CD83") )
VlnPlot(CARO.combined_1, features = c("CD36","CCR2","CD14","CD61","CX3CR1","CCL5","PTPRC") )
VlnPlot(CARO.combined_1, features = c("HNRNPU","MALAT1","CD14","ICAM1","CX3CR1","CCL5","PTPRC") )
VlnPlot(CARO.combined_1, features = c("HDC","KIT","TPSAB1","GZMB","CCL5","TPSB2") )
VlnPlot(CARO.combined_1, features = c("ENPP3","MITF","CD63","ITGA2","CD33","ITGAM") )
VlnPlot(CARO.combined_1, features = c("TNF","TNFSF13","CXCL8","CXCL3","CCL3") )
VlnPlot(CARO.combined_1, features = c("CD9","MRC1","FCER2","CCL22","AREG","EREG") )
VlnPlot(CARO.combined_1, features = c("MKI67","TUBB","STMN1","TYMS","NFKB1","STAT4") )
VlnPlot(CARO.combined_1, features = c("C1QA","C1QB","C1QC","TYMS","NFKB1","STAT4") )
VlnPlot(CARO.combined_1, features = c("TET2","DNMT3A","ASXL1","PPM1D","NLRP3","STING1","CGAS" ))
VlnPlot(CARO.combined_1, features = c("OLR1","SCARB1","SCARB2","MSR1","MARCO","SRSF2") )
VlnPlot(CARO.combined_1, features = c("NR1H3","CTSD","CTSL","SPP1","MARCO","FABP4") )
VlnPlot(CARO.combined_1, features = c("TET2","TET1","TET3","SPP1","MARCO","SLC25A44") )
VlnPlot(CARO.combined_1, features = c("TNF","TNFSF13"),split.by = "stim",  cols = c("sa" ="deepskyblue", "sym"="red"))
VlnPlot(CARO.combined_1, features = c("CD163","F13A1","MS4A4A","SPP1","MARCO","SLC25A44") )

#Tcells--------
CARO.combined_2 <- subset(CARO.combined, idents = c(0,1,2,3,4))
CARO.combined_2 <- ScaleData(CARO.combined_2, verbose = FALSE)
CARO.combined_2<- RunPCA(CARO.combined_2, npcs = 30, verbose = FALSE)
CARO.combined_2 <- RunUMAP(CARO.combined_2 , reduction = "pca", dims = 1:30)
CARO.combined_2 <- FindNeighbors(CARO.combined_2, reduction = "pca", dims = 1:30)
CARO.combined_2 <- FindClusters(CARO.combined_2, resolution = 0.5)
p1 <- DimPlot(CARO.combined_2, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(CARO.combined_2, reduction = "umap", label = FALSE, repel = TRUE)
p2

levels(CARO.combined_2) # [1] "0" "1" "2" "3" "4" 
levels(CARO.combined_2) <- c("0", "2", "1", "3", "4")

new.cluster.ids <- c("0","1","2","3","4")
names(new.cluster.ids) <- levels(CARO.combined_2)
CARO.combined_2<- RenameIdents(CARO.combined_2, new.cluster.ids)


VlnPlot(CARO.combined_2, features = c("CD3E","CD4","CD8A","IL7R","GZMB","TYROBP") )
VlnPlot(CARO.combined_2, features = c("GZMA","GZMK","PRF1","LEF1","LEF1","SELL") )
VlnPlot(CARO.combined_2, features = c("RORA","GATA3","CD40LG","IL10","IL4","IFNG") )
VlnPlot(CARO.combined_2, features = c("IL10","IL4","IL6","IFNG","IL17A") )
VlnPlot(CARO.combined_2, features = c("CCR7","GATA3","RORC","FOXP3") )
VlnPlot(CARO.combined_2, features = c("GZMB","TBX21","NKG7","GNLY","CX3CR1","CD69") )

VlnPlot(CARO.combined_2, features = c("RORA","GATA3","PD1","CD40LG","IL10","IL4","LILRB1R") ,split.by = "stim")

CARO.combined.markers_2 <- FindAllMarkers(CARO.combined_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CARO.combined.markers_2 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- CARO.combined.markers_2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(CARO.combined_2, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)

#CD8----

CARO.combined_2_1 <- subset(CARO.combined_2, idents = c(0,1))
CARO.combined_2_1 <- ScaleData(CARO.combined_2_1, verbose = FALSE)
CARO.combined_2_1<- RunPCA(CARO.combined_2_1, npcs = 30, verbose = FALSE)
CARO.combined_2_1 <- RunUMAP(CARO.combined_2_1, reduction = "pca", dims = 1:30)
CARO.combined_2_1 <- FindNeighbors(CARO.combined_2_1, reduction = "pca", dims = 1:30)
CARO.combined_2_1 <- FindClusters(CARO.combined_2_1, resolution = 0.5)
p1 <- DimPlot(CARO.combined_2_1, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(CARO.combined_2_1, reduction = "umap", label = FALSE, repel = TRUE)
p2

VlnPlot(CARO.combined_2_1, features = c("GZMA","GZMK","GZMB","TBX21","CX3CR1","NKG7") )
VlnPlot(CARO.combined_2_1, features = c("CD69","IL7R","IFNG","CD103","SELL","CD27") )
VlnPlot(CARO.combined_2_1, features = c("TNF","IL2","IFNG","EOMES","SELL","CD27") )
VlnPlot(CARO.combined_2_1, features = c("TNF","CD27","IFNG","MKI67","CD74","CXCR6") )


CARO.combined.markers_2_1 <- FindAllMarkers(CARO.combined_2_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CARO.combined.markers_2_1 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- CARO.combined.markers_2_1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(CARO.combined_2_1, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)



#CD4----

CARO.combined_2_2 <- subset(CARO.combined_2, idents = c(2,3))
CARO.combined_2_2 <- ScaleData(CARO.combined_2_2, verbose = FALSE)
CARO.combined_2_2<- RunPCA(CARO.combined_2_2, npcs = 30, verbose = FALSE)
CARO.combined_2_2 <- RunUMAP(CARO.combined_2_2, reduction = "pca", dims = 1:30)
CARO.combined_2_2 <- FindNeighbors(CARO.combined_2_2, reduction = "pca", dims = 1:30)
CARO.combined_2_2 <- FindClusters(CARO.combined_2_2, resolution = 0.7)

p1 <- DimPlot(CARO.combined_2_2, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(CARO.combined_2_2, reduction = "umap", label = FALSE, repel = TRUE)
p2

VlnPlot(CARO.combined_2_2, features = c("GZMA","GZMK","PRF1","LEF1","IL7R","SELL") )
VlnPlot(CARO.combined_2_2, features = c("TBX21","GATA3","CCR7","RORC","FOXP3","CTLA4") )
VlnPlot(CARO.combined_2_2, features = c("IL17A","IL23","IL2","IL10","IL4","IFNG") )
VlnPlot(CARO.combined_2_2, features = c("GZMB","CD28","PRF1","CD69","IL7R","SELL") )
VlnPlot(CARO.combined_2_2, features = c("CXCR3","CCR4","CCR6","TER1"), split.by ="stim",cols = c("sa" ="deepskyblue", "sym"="red") )

VlnPlot(CARO.combined_2_2, features = c("CCL4","CCL5","TNFAIP3","RORC","FOXP3","CTLA4") )


FeaturePlot(CARO.combined_2_2, features=c('CD69'), min.cutoff=2.0, max.cutoff='q90',split.by = "stim" )
FeaturePlot(CARO.combined_2_2, features=c('CD69'), min.cutoff=2.0, max.cutoff='q90')



levels(CARO.combined_2_2) # [1] "0" "1" "2" "3" "4" 
levels(CARO.combined_2_2) <- c("0", "1", "3", "2", "4")

new.cluster.ids <- c("0","1","2","3","4")
names(new.cluster.ids) <- levels(CARO.combined_2_2)
CARO.combined_2_2<- RenameIdents(CARO.combined_2_2, new.cluster.ids)


CARO.combined.markers_2_2 <- FindAllMarkers(CARO.combined_2_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CARO.combined.markers_2_2 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- CARO.combined.markers_2_2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(CARO.combined_2_2, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)



#Bcells---------
CARO.combined_3 <- subset(CARO.combined, idents = c(8,9))
CARO.combined_3 <- ScaleData(CARO.combined_3, verbose = FALSE)
CARO.combined_3 <- RunPCA(CARO.combined_3, npcs = 30, verbose = FALSE)
CARO.combined_3 <- RunUMAP(CARO.combined_3 , reduction = "pca", dims = 1:30)
CARO.combined_3 <- FindNeighbors(CARO.combined_3 , reduction = "pca", dims = 1:30)
CARO.combined_3 <- FindClusters(CARO.combined_3, resolution = 0.5)
p1 <- DimPlot(CARO.combined_3, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(CARO.combined_3, reduction = "umap", label = FALSE, repel = TRUE)
p2 

VlnPlot(CARO.combined_3, features = c("CD79A","CD79B","CD27","FCER2","CD22") )
VlnPlot(CARO.combined_3, features = c("CD38","PTPRC","CD40","CD24","CD79B","FCER2") )
VlnPlot(CARO.combined_3, features = c("IGHG1","IGHG3","IGHG","IGHM","IGHD","IGHE"))

VlnPlot(CARO.combined_3, features = c("IGHG1","IGHG3","IGHG","IGHM","IGHD","IGHE"),split.by = "stim",  cols = c("sa" ="deepskyblue", "sym"="red"))
VlnPlot(CARO.combined_3, features = c("CD40","IL4R","IL21R","PAX5","BACH2","ICOSL","CD80","CD86","HLA-DQB1","HLA-DQB1","HLA-DPA1"), split.by = "stim",  cols = c("sa" ="deepskyblue", "sym"="red"))
VlnPlot(CARO.combined_3, features = c("IGHG1","IGHG2","IGHG3","IGLG1"),split.by = "stim",  cols = c("sa" ="deepskyblue", "sym"="red"))


CARO.combined_3_1 <- subset(CARO.combined_3, idents = c(1))
CARO.combined_3_1 <- ScaleData(CARO.combined_3_1, verbose = FALSE)
CARO.combined_3_1<- RunPCA(CARO.combined_3_1, npcs = 30, verbose = FALSE)
CARO.combined_3_1 <- RunUMAP(CARO.combined_3_1 , reduction = "pca", dims = 1:30)
CARO.combined_3_1 <- FindNeighbors(CARO.combined_3_1 , reduction = "pca", dims = 1:30)
CARO.combined_3_1 <- FindClusters(CARO.combined_3_1, resolution = 0.5)
p1 <- DimPlot(CARO.combined_3_1, reduction = "umap", split.by = "stim")
p1

# heatmap------------------------------------------------------------------------------
CARO.combined_1.markers <- FindAllMarkers(CARO.combined_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CARO.combined_1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- CARO.combined_1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(CARO.combined_1, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)

CARO.combined_1 <- readRDS (file =  "/home/cardiovascular/Desktop/CARO_human/RDS_Myeloid_cells/CARO_Myeoloid_cluster.rds")



# VInplot----
VlnPlot(CARO.combined_1, features = c("CXCL10","S100A8","HMGB2","HLA-DPA1","HLA-DQA1","HLA-DPB1") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))
VlnPlot(CARO.combined_1, features = c("APOE") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))
VlnPlot(CARO.combined_1, features = c("CXCL10") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))

CARO.combined_1_DC <- subset(CARO.combined_1, idents = c("5_sa", "5_sym"))
VlnPlot(CARO.combined_1_DC, features = c("CD40","CD80","CD86","HLA-DPA1","HLA-DQA1","HLA-DPB1") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))
VlnPlot(CARO.combined_1_DC, features = c("CLEC10A","FCER1A","CD1C","IL3RA","RELB") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))
VlnPlot(CARO.combined_1_DC, features = c("NECTIN2","PVR","TNFSF4","IL3RA","RELB") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))


CARO.combined_2_2_cmCD4 <- subset(CARO.combined_2_2, idents = c("3_sa", "3_sym"))
VlnPlot(CARO.combined_2_2_cmCD4, features = c("CD28","CD40LG","CTLA4","TIGHT","FOXP3","TNFRSF4") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))
VlnPlot(CARO.combined_2_2_cmCD4, features = c("TIGIT","CD226" ,"IL17A","IL4","IL10") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))

CARO.combined_2_2_eCD4 <- subset(CARO.combined_2_2, idents = c("2_sa", "2_sym"))
VlnPlot(CARO.combined_2_eCD4, features = c("CD28","CD40LG","PDCD1","LAG3","CTLA4","TNFRSF4") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))
VlnPlot(CARO.combined_2_eCD4, features = c("TIGIT","CD226","HAVCR2","PDCD1") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))

CARO.combined_2_2_nCD4 <- subset(CARO.combined_2_2, idents = c("0_sa","3_sa","1_sa", "0_sym", "1_sym", "3_sym"))
VlnPlot(CARO.combined_2_2_nCD4, features = c("CD28","CD40LG","CTLA4","LAG3","CTLA4","TNFRSF4") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))
VlnPlot(CARO.combined_2_2_nCD4, features = c("TIGIT","CD226","HAVCR2","PDCD1") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))

CARO.combined_2_2_e1CD4 <- subset(CARO.combined_2_2, idents = c("0_sa", "0_sym"))
VlnPlot(CARO.combined_2_2_e1CD4, features = c("CD28","CD40LG","PDCD1","LAG3","CTLA4","TNFRSF4") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))
VlnPlot(CARO.combined_2_2_e1CD4, features = c("TIGIT","CD226","HAVCR2","PDCD1") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))

CARO.combined_2_2_regCD4 <- subset(CARO.combined_2_2, idents = c("4_sa", "4_sym"))
VlnPlot(CARO.combined_2_2_regCD4, features = c("CD28","CD40LG","PDCD1","LAG3","CTLA4","TNFRSF4") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))
VlnPlot(CARO.combined_2_2_regCD4, features = c("TIGIT","CD226","HAVCR2","PDCD1") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))


CARO.combined_2.markers <- FindAllMarkers(CARO.combined_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CARO.combined_2.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- CARO.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(CARO.combined_2, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


CARO.combined_3.markers <- FindAllMarkers(CARO.combined_3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CARO.combined_3.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- CARO.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(CARO.combined_3, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)

# volcano all myeloid------

CARO.combined$celltype <- Idents(CARO.combined)
CARO.combined$celltype.stim <- paste(Idents(CARO.combined), CARO.combined$stim, sep="_")
Idents(CARO.combined) <- "celltype.stim" 
levels(CARO.combined) 


CARO.combined_Myeloid_sa<- subset(CARO.combined, idents = c("5_sa","6_sa","7_sa"))
CARO.combined_Myeloid_sym<- subset(CARO.combined, idents = c("5_sym","6_sym","7_sym"))

CARO.combined_Tcells_sa<- subset(CARO.combined, idents = c("0_sa","1_sa","2_sa","3_sa","4_sa"))
CARO.combined_Tcells_sym<- subset(CARO.combined, idents = c("0_sym","1_sym","2_sym","3_sym","4_sym"))

CARO.combined_Bcells_sa<- subset(CARO.combined, idents = c("8_sa","9_sa"))
CARO.combined_Bcells_sym<- subset(CARO.combined, idents = c("8_sym","9_sym"))



CARO.table_myeloid <- FindMarkers(CARO.combined, ident.1 = c("5_sa","6_sa","7_sa"), ident.2 =c("5_sym","6_sym","7_sym"), verbose = FALSE, logfc.threshold = 0)

CARO.table_myeloid$logp <- -log10(CARO.table_myeloid$p_val)

CARO.table_myeloid_filtered_left = subset(CARO.table_myeloid, logp>=5 & avg_log2FC <= -1.0)
CARO.table_myeloid_filtered_right = subset(CARO.table_myeloid, logp>=5 & avg_log2FC >= 1.0)

genes.to.label.left <- rownames(CARO.table_myeloid_filtered_left)
genes.to.label.right <- rownames(CARO.table_myeloid_filtered_right)

p1 <- ggplot(CARO.table_myeloid, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="blue", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="red", repel = TRUE, xnudge=0)
p1

# VOLCANO Trem2 macs----



# volcano all CD4------

CARO.combined_2$celltype <- Idents(CARO.combined_2)
CARO.combined_2$celltype.stim <- paste(Idents(CARO.combined_2), CARO.combined_2$stim, sep="_")
Idents(CARO.combined_2) <- "celltype.stim" 
levels(CARO.combined_2) 

CARO.table_CD4 <- FindMarkers(CARO.combined_2, ident.1 = c("2_sa"), ident.2 =c("2_sym"), verbose = FALSE, logfc.threshold = 0)

CARO.table_CD4$logp <- -log10(CARO.table_CD4$p_val)

CARO.table_CD4_filtered_left = subset(CARO.table_CD4, logp>=2 & avg_log2FC <= -0.5)
CARO.table_CD4_filtered_right = subset(CARO.table_CD4, logp>=2 & avg_log2FC >= 0.5)

genes.to.label.left <- rownames(CARO.table_CD4_filtered_left)
genes.to.label.right <- rownames(CARO.table_CD4_filtered_right)

p1 <- ggplot(CARO.table_CD4, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="blue", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="red", repel = TRUE, xnudge=0)
p1

write.table(CARO.table_CD4, "effecCD4.txt", quote=F, col.names=F, append=T)

# central memory CD4

CARO.combined_2_2$celltype <- Idents(CARO.combined_2_2)
CARO.combined_2_2$celltype.stim <- paste(Idents(CARO.combined_2_2), CARO.combined_2_2$stim, sep="_")
Idents(CARO.combined_2_2) <- "celltype.stim" 
levels(CARO.combined_2_2) 

CARO.table_cCD4 <- FindMarkers(CARO.combined_2_2, ident.1 = c("3_sa"), ident.2 =c("3_sym"), verbose = FALSE, logfc.threshold = 0)

CARO.table_cCD4$logp <- -log10(CARO.table_cCD4$p_val)

CARO.table_cCD4_filtered_left = subset(CARO.table_cCD4, logp>=2 & avg_log2FC <= -0.4)
CARO.table_cCD4_filtered_right = subset(CARO.table_cCD4, logp>=2 & avg_log2FC >= 0.4)

genes.to.label.left <- rownames(CARO.table_cCD4_filtered_left)
genes.to.label.right <- rownames(CARO.table_cCD4_filtered_right)

p1 <- ggplot(CARO.table_cCD4, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="blue", repel = TRUE, xnudge=-0.1)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="red", repel = TRUE, xnudge=-0.1)
p1

write.table(CARO.table1, "table.txt", quote=F, col.names=F, append=T)

# volcano myeloid subset
CARO.combined_1$celltype <- Idents(CARO.combined_1)
CARO.combined_1$celltype.stim <- paste(Idents(CARO.combined_1), CARO.combined_1$stim, sep="_")
Idents(CARO.combined_1) <- "celltype.stim" 
levels(CARO.combined_1)

CARO.table1 <- FindMarkers(CARO.combined_1, ident.1 = "5_sym", ident.2 = "5_sa", verbose = FALSE, logfc.threshold = 0)
CARO.table2 <- FindMarkers(CARO.combined_1, ident.1 = "0_sym", ident.2 = "0_sa", verbose = FALSE, logfc.threshold = 0)
CARO.table3 <- FindMarkers(CARO.combined_1, ident.1 = "1_sym", ident.2 = "1_sa", verbose = FALSE, logfc.threshold = 0)

CARO.table1$logp <- -log10(CARO.table1$p_val)
CARO.table2$logp <- -log10(CARO.table2$p_val)
CARO.table3$logp <- -log10(CARO.table3$p_val)
CARO.table4$logp <- -log10(CARO.table4$p_val)

# cluster1----
CARO.table1_filtered_left = subset(CARO.table1, logp>=1.30103& avg_log2FC <= -0.8)
CARO.table1_filtered_right = subset(CARO.table1, logp>=1.30103 & avg_log2FC >= 0.8)

genes.to.label.left <- rownames(CARO.table1_filtered_left)
genes.to.label.right <- rownames(CARO.table1_filtered_right)

p1 <- ggplot(CARO.table1, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="red", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="blue", repel = TRUE, xnudge=0)
p1

write.table(CARO.table1, "table.txt", quote=F, col.names=F, append=T)

VlnPlot(CARO.combined_1, ident.1 = "5_sa", ident.2 = "5_sym", features = c("CXCL10") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))

# cluster2----
CARO.table2_filtered_left = subset(CARO.table2, logp>=2 & avg_log2FC <= -0.8)
CARO.table2_filtered_right = subset(CARO.table2, logp>=2 & avg_log2FC >= 0.8)

genes.to.label.left <- rownames(CARO.table2_filtered_left)
genes.to.label.right <- rownames(CARO.table2_filtered_right)

p1 <- ggplot(CARO.table2, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="red", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="blue", repel = TRUE, xnudge=0)
p1

write.table(CARO.table2, "table.txt", quote=F, col.names=F, append=T)

# cluster3----

CARO.table3_filtered_left = subset(CARO.table3, logp>=2 & avg_log2FC <= -0.8)
CARO.table3_filtered_right = subset(CARO.table3, logp>=2 & avg_log2FC >= 0.8)

genes.to.label.left <- rownames(CARO.table3_filtered_left)
genes.to.label.right <- rownames(CARO.table3_filtered_right)

p1 <- ggplot(CARO.table3, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="red", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="blue", repel = TRUE, xnudge=0)
p1

write.table(CARO.table2, "table.txt", quote=F, col.names=F, append=T)


plot1 <- FeatureScatter(CARO.combined_3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CARO.combined_3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

CARO.combined_3

CARO.combined_Mye2

CARO.combined_Mye2$celltype <- Idents(CARO.combined_Mye2)
CARO.combined_Mye2$celltype.stim <- paste(Idents(CARO.combined_Mye2), CARO.combined_Mye2$stim, sep="_")
Idents(CARO.combined_Mye2) <- "celltype.stim" 
levels(CARO.combined_Mye2) 

options(future.globals.maxSize = 8000 * 1024^2)
CARO.table_Mye2 <- FindMarkers(CARO.combined_Mye2, ident.1 = c("0_sym","1_sym","2_sym","3_sym","4_sym","5_sym","6_sym","7_sym","8_sym","9_sym","10_sym","11_sym"), ident.2 =c("0_asym","1_asym","2_asym","3_asym","4_asym","5_asym","6_asym","7_asym","8_asym","9_asym","10_asym","11_asym"), verbose = FALSE, logfc.threshold = 0)

CARO.table_Mye2$logp <- -log10(CARO.table_Mye2$p_val)

CARO.table_Mye2_filtered_left = subset(CARO.table_Mye2, logp>=1 & avg_log2FC <= -0.1)
CARO.table_Mye2_filtered_right = subset(CARO.table_Mye2, logp>=1 & avg_log2FC >= 0.1)

genes.to.label.left <- rownames(CARO.table_Mye2_filtered_left)
genes.to.label.right <- rownames(CARO.table_Mye2_filtered_right)

p1 <- ggplot(CARO.table_Mye2, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="blue", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="red", repel = TRUE, xnudge=0)
p1

CARO.table_Mye2["MMP19",1:6]
CARO.table_Mye2["MMP9",1:6]

