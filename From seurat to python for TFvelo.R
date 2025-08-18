library(SeuratData)
library(MuDataSeurat)
library(SeuratDisk)
library(patchwork)

library(reticulate)

use_condaenv("r-reti2")
repl_python()

#seurat to python
seurat <- readRDS ("~/Desktop/Human_DCA/RDS/DCA_Myeloid.rds")

DefaultAssay(seurat) <- "RNA"

DimPlot(seurat, reduction = "umap", label = TRUE, repel = TRUE) 

seurat <- AddMetaData(seurat, seurat@active.ident, col.name = "seurat_clusters")


#color code
library(scales)
show_col(hue_pal()(6))

#seurat_disk
seurat <- AddMetaData(seurat, seurat@active.ident, col.name = "seurat_clusters")

SaveH5Seurat(seurat, filename = "~/Desktop/Human_DCA/anndata/DCA_Mac.h5seurat", overwrite = TRUE)
Convert("~/Desktop/Human_DCA/anndata/DCA_Mac.h5seurat", dest = "h5ad",assay="RNA", overwrite = TRUE)

