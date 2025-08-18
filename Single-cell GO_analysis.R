#install.packages('Seurat')
#install.packages('Rcpp')
#update.packages()
library(Rcpp)
library(Seurat)
library(dplyr)

data <- readRDS (file =  "~/Desktop/Human_carotid_public/RDS/Human_carotid_sample2_Mac_240918.rds")

type0.markers <- FindMarkers(data, ident.1 = "0", ident.2 = c("1","2","3","4","5","7","10","11"))
type0.markers <- type0.markers[type0.markers$p_val_adj < 0.05&type0.markers$avg_log2FC>0,]
type1.markers <- FindMarkers(data, ident.1 = "1", ident.2 = c("0","2","3","4","5","7","10","11"))
type1.markers <- type1.markers[type1.markers$p_val_adj < 0.05&type1.markers$avg_log2FC>0,]
type2.markers <- FindMarkers(data, ident.1 = "2", ident.2 = c("0","1","2","3","4","5","7","10","11"))
type2.markers <- type2.markers[type2.markers$p_val_adj < 0.05&type2.markers$avg_log2FC>0,]
type3.markers <- FindMarkers(data, ident.1 = "3", ident.2 = c("0","1","2","4","5","7","10","11"))
type3.markers <- type3.markers[type3.markers$p_val_adj < 0.05&type3.markers$avg_log2FC>0,]
type4.markers <- FindMarkers(data, ident.1 = "4", ident.2 = c("0","1","2","3","5","7","10","11"))
type4.markers <- type4.markers[type4.markers$p_val_adj < 0.05&type4.markers$avg_log2FC>0,]
type5.markers <- FindMarkers(data, ident.1 = "5", ident.2 = c("0","1","2","3","4","7","10","11"))
type5.markers <- type5.markers[type5.markers$p_val_adj < 0.05&type5.markers$avg_log2FC>0,]
type7.markers <- FindMarkers(data, ident.1 = "7", ident.2 = c("0","1","2","3","4","5","10","11"))
type7.markers <- type7.markers[type7.markers$p_val_adj < 0.05&type7.markers$avg_log2FC>0,]
type10.markers <- FindMarkers(data, ident.1 = "10", ident.2 = c("0","1","2","3","4","5","7","11"))
type10.markers <- type10.markers[type10.markers$p_val_adj < 0.05&type10.markers$avg_log2FC>0,]
type11.markers <- FindMarkers(data, ident.1 = "11", ident.2 = c("0","1","2","3","4","5","7","10"))
type11.markers <- type11.markers[type11.markers$p_val_adj < 0.05&type11.markers$avg_log2FC>0,]

library(org.Hs.eg.db)
hs <- org.Hs.eg.db


library(clusterProfiler)
type0.gene_SYMBOLs <- rownames(type0.markers)
type0.gene_IDs <- AnnotationDbi::select(hs, keys=type0.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID

type1.gene_SYMBOLs <- rownames(type1.markers)
type1.gene_IDs <- AnnotationDbi::select(hs, keys=type1.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID

type2.gene_SYMBOLs <- rownames(type2.markers)
type2.gene_IDs <- AnnotationDbi::select(hs, keys=type2.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID

type3.gene_SYMBOLs <- rownames(type3.markers)
type3.gene_IDs <- AnnotationDbi::select(hs, keys=type3.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID

type4.gene_SYMBOLs <- rownames(type4.markers)
type4.gene_IDs <- AnnotationDbi::select(hs, keys=type4.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID

type5.gene_SYMBOLs <- rownames(type5.markers)
type5.gene_IDs <- AnnotationDbi::select(hs, keys=type5.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID

type7.gene_SYMBOLs <- rownames(type7.markers)
type7.gene_IDs <- AnnotationDbi::select(hs, keys=type7.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID

type10.gene_SYMBOLs <- rownames(type10.markers)
type10.gene_IDs <- AnnotationDbi::select(hs, keys=type10.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID

type11.gene_SYMBOLs <- rownames(type11.markers)
type11.gene_IDs <- AnnotationDbi::select(hs, keys=type11.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID


# compare cluster

type0.gene_SYMBOLs <- rownames(type0.markers)
type0.gene_IDs <- AnnotationDbi::select(hs, keys=type0.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID
type1.gene_SYMBOLs <- rownames(type1.markers)
type1.gene_IDs <- AnnotationDbi::select(hs, keys=type1.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID
type2.gene_SYMBOLs <- rownames(type2.markers)
type2.gene_IDs <- AnnotationDbi::select(hs, keys=type2.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID
type3.gene_SYMBOLs <- rownames(type3.markers)
type3.gene_IDs <- AnnotationDbi::select(hs, keys=type3.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID
type4.gene_SYMBOLs <- rownames(type4.markers)
type4.gene_IDs <- AnnotationDbi::select(hs, keys=type4.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID
type5.gene_SYMBOLs <- rownames(type5.markers)
type5.gene_IDs <- AnnotationDbi::select(hs, keys=type5.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID
type7.gene_SYMBOLs <- rownames(type7.markers)
type7.gene_IDs <- AnnotationDbi::select(hs, keys=type7.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID
type10.gene_SYMBOLs <- rownames(type10.markers)
type10.gene_IDs <- AnnotationDbi::select(hs, keys=type10.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID
type11.gene_SYMBOLs <- rownames(type11.markers)
type11.gene_IDs <- AnnotationDbi::select(hs, keys=type11.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID


genelist <- list(type0.gene_IDs, type1.gene_IDs, type2.gene_IDs, type3.gene_IDs, type4.gene_IDs, type5.gene_IDs, type7.gene_IDs, type10.gene_IDs, type11.gene_IDs)
names(genelist) <- c("type0", "type1", "type2", "type3", "type4", "type5", "type7", "type10", "type11")

# GO:BP
cgBP <- compareCluster(geneCluster = genelist, fun = enrichGO, ont="BP",OrgDb='org.Hs.eg.db',pvalueCutoff = 1)
dotplot(cgBP,showCategory = 2) + 
  theme(axis.text.y = element_text(size = 9, lineheight = 0.7),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

any("organelle fission" == cgBP@compareClusterResult[["Description"]])

# GO:MF

cgMF <- compareCluster(geneCluster = genelist, fun = enrichGO, ont="MF",OrgDb='org.Hs.eg.db')
dotplot(cgMF,showCategory = 2) + 
  theme(axis.text.y = element_text(size = 10, lineheight = 0.9),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# GO:CC
cgCC <- compareCluster(geneCluster = genelist, fun = enrichGO, ont="CC",OrgDb='org.Hs.eg.db')
dotplot(cgCC,showCategory = 2) + 
  theme(axis.text.y = element_text(size = 10, lineheight = 0.9),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

sets2plot <- c("retinol metabolic process","organelle fission", "chromosome segregation","mitotic nuclear division","DNA replication")

dotplot(cgBP,showCategory = set2plot)

# GO:BPで、type1 empty
p1 <- dotplot(cgBP)
df <- data.frame("type1\n(39)","GO:0071346","cellular response to interferon-gamma",NA,"1/18866",2.073631e-17,6.26029e-14,4.793361e-14,"3120/3119/3113/3123/3115/3122/4261/3127/3117/115361/3717/3118/2210/2934/2209/6352/6349/6355/6772/6348",20)
names(df) <- c("Cluster","ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count")
p1$data <- rbind(p1$data, df)
p1$data$Cluster <- factor(p1$data$Cluster, levels=c("type0\n(150)","type1\n(39)","type2\n(122)"))
p1


