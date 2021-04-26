library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)

DefaultAssay(bm) <- "RNA"
SaveH5Seurat(bm, filename = "mouseBM.h5Seurat")
Convert("mouseBM.h5Seurat", dest = "h5ad")