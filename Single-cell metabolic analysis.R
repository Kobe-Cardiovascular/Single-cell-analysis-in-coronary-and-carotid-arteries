#library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scMetabolism)
library(rsvd)

seurat <- readRDS ("~/Desktop/Human_carotid_public/RDS/Human_carotid_sample2_Mye_240824.rds")

DimPlot(seurat, reduction = "umap", label = TRUE, repel = TRUE)


asym <- subset(seurat, stim == "asym")
sym <- subset(seurat, stim == "sym")

seurat <- subset(seurat, idents = c(10), invert = TRUE)

#8.6GB
DefaultAssay(seurat) <- "RNA"
countexp.Seurat<-sc.metabolism.Seurat(obj = seurat, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG")
#saveRDS(countexp.Seurat,"~/Desktop/Human_DCA/RDS/DCA_Myeloid_scmetab_240516.rds")

saveRDS(countexp.Seurat,"~/Desktop/Human_carotid_public/RDS/Human_carotid_sample2_Mye_invert10_scmetab_250402.RDS")

countexp.Seurat <- readRDS("~/Desktop/Human_carotid_public/RDS/Human_carotid_sample2_Mye_invert10_scmetab_250402.RDS")


input.pathway<-c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)","Folate biosynthesis",
                 "Fatty acid biosynthesis","Fatty acid elongation","Fatty acid degradation",
                 "Oxidative phosphorylation", )

input.pathway<-c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)","Folate biosynthesis",
                 "Fatty acid biosynthesis","Fatty acid elongation","Fatty acid degradation","Oxidative phosphorylation",
                 "Fructose and mannose metabolism", "Tyrosine metabolism", "Riboflavin metabolism", "Phosphonate and phosphinate metabolism",
                 "Porphyrin and chlorophyll metabolism","Folate biosynthesis", "Glycosaminoglycan degradation",
                 "Amino sugar and nucleotide sugar metabolism", "Ascorbate and aldarate metabolism",
                 "Other glycan degradation","Pentose and glucuronate interconversions","Phenylalanine metabolism")

input.pathway<-c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)","Folate biosynthesis",
                 "Fatty acid biosynthesis","Fatty acid elongation","Fatty acid degradation","Oxidative phosphorylation",
                 "Folate biosynthesis", "Glycosaminoglycan degradation")

input.pathway<-c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)","Folate biosynthesis",
                 "Fatty acid degradation","Oxidative phosphorylation",
                 "Folate biosynthesis", "Glycosaminoglycan degradation")


DimPlot(countexp.Seurat, reduction = "umap", label = TRUE, repel = TRUE,raster = FALSE)
DimPlot(countexp.Seurat, reduction = "umap", split.by = "sample")
DimPlot(countexp.Seurat, reduction = "umap", split.by = "stim")

DimPlot.metabolism(obj = countexp.Seurat, pathway = "Glycolysis / Gluconeogenesis", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 0.5)

limit_plot(obj = countexp.Seurat, pathway = "Citrate cycle (TCA cycle)", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 0.5, limits_scale = c(0,1.6))
limit_plot(obj = countexp.Seurat, pathway = "Glycolysis / Gluconeogenesis", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 0.5, limits_scale = c(0,2.2))



paste(countexp.Seurat$stim, countexp.Seurat@active.ident, sep="_")

new.cluster.ids <- c("00","01","02","03","04","05","06","07","08","09","10","11")
names(new.cluster.ids) <- levels(countexp.Seurat)
countexp.Seurat<- RenameIdents(countexp.Seurat, new.cluster.ids)

new.cluster.ids <- c("00","01","02","03","04","05","06","07","08","09","11")
names(new.cluster.ids) <- levels(countexp.Seurat)
countexp.Seurat<- RenameIdents(countexp.Seurat, new.cluster.ids)

countexp.Seurat<- AddMetaData(countexp.Seurat , countexp.Seurat@active.ident, col.name = "ident")

countexp.Seurat<- AddMetaData(countexp.Seurat , paste(countexp.Seurat$stim, countexp.Seurat@active.ident, sep=" Mye."), col.name = "stim_ident")

DotPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "sample", norm = "y")
DotPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "stim", norm = "y")
DotPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "ident", norm = "y")
DotPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "stim_ident", norm = "y")

saveRDS(countexp.Seurat,"~/Desktop/Sarcoidosis_project/RDS/Sar_ICM_Mye_asym_scmetab_0507.RDS")


limit_plot <- function (obj, pathway, dimention.reduction.type = "umap", dimention.reduction.run = T, 
          size = 1, limits_scale = c(0,3)) 
{
  if (dimention.reduction.type == "umap") {
    if (dimention.reduction.run == T) 
      obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = 1:40)
    umap.loc <- obj@reductions$umap@cell.embeddings
    row.names(umap.loc) <- colnames(obj)
    signature_exp <- obj@assays$METABOLISM$score
    input.pathway <- pathway
    signature_ggplot <- data.frame(umap.loc, t(signature_exp[input.pathway, 
    ]))
    library(wesanderson)
    pal <- wes_palette("Zissou1", 100, type = "continuous")
    library(ggplot2)
    plot <- ggplot(data = signature_ggplot, aes(x = UMAP_1, 
                                                y = UMAP_2, color = signature_ggplot[, 3])) + geom_point(size = size) + 
      scale_fill_gradientn(colours = pal,limits = limits_scale) + scale_color_gradientn(colours = pal,limits = limits_scale) + 
      labs(color = input.pathway) + xlab("UMAP 1") + ylab("UMAP 2") + 
      theme(aspect.ratio = 1) + theme(panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                      axis.line = element_line(colour = "black"))
  }
  if (dimention.reduction.type == "tsne") {
    if (dimention.reduction.run == T) 
      obj <- Seurat::RunTSNE(obj, reduction = "pca", dims = 1:40)
    tsne.loc <- obj@reductions$tsne@cell.embeddings
    row.names(tsne.loc) <- colnames(obj)
    signature_exp <- obj@assays$METABOLISM$score
    input.pathway <- pathway
    signature_ggplot <- data.frame(tsne.loc, t(signature_exp[input.pathway, 
    ]))
    pal <- wes_palette("Zissou1", 100, type = "continuous")
    library(ggplot2)
    plot <- ggplot(data = signature_ggplot, aes(x = tSNE_1, 
                                                y = tSNE_2, color = signature_ggplot[, 3])) + geom_point(size = size) + 
      scale_fill_gradientn(colours = pal) + scale_color_gradientn(colours = pal) + 
      labs(color = input.pathway) + xlab("tSNE 1") + ylab("tSNE 2") + 
      theme(aspect.ratio = 1) + theme(panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                      axis.line = element_line(colour = "black"))
  }
  plot
}




