#!/usr/bin/env Rscript

# Load packages
library(getopt)
library(Seurat)
library(future)
library(cowplot)
library(clustree)
library(gridExtra)
library(grid)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(scHelper)

spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    plot_path = "./output/NF-downstream_analysis_stacas/seurat/cell_state_classification/plots/"
    rds_path = "./output/NF-downstream_analysis_stacas/seurat/cell_state_classification/rds_files/"
    data_path = "./output/NF-downstream_analysis_stacas/seurat/6_contamination_filt/rds_files/"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multiprocess", workers = ncores)
    options(future.globals.maxSize = 16* 1024^3) # 32gb
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

seurat_data <- readRDS(list.files(data_path, full.names = TRUE))

# Set RNA to default assay for plotting expression data
DefaultAssay(seurat_data) <- "RNA"

########################################################################################################
#                                      Cell state classification                                    #
########################################################################################################

###################### Nerual crest and placode cell type identification ###################### 
plac_nc_genes <- c(
  # delaminating NC
  "ETS1", "SOX10", "SOX8", "LMO4",  "TFAP2B", "SOX9",
  # NPB
  "DRAXIN", "TFAP2A", "MSX1", "CSRNP1", "PAX7", "BMP5", "MSX2",
  # NC
  "WNT6",
  # Placodes
  "PITX1", "PITX2", "ZNF385C",  "SIX1", "EYA2", "DLX6", "HOMER2", 'SIX3', 'SHISA2',
  # Epidermis
  "KRT18", "KRT7", "GATA2"
)

ncol = ceiling((length(plac_nc_genes)+1)/8)+1
nrow = ceiling((length(plac_nc_genes)+1)/ncol)

# plot expression of NC and placodal genes
png(paste0(plot_path, 'multi_feature_plac_nc.png'), width = ncol*10, height = nrow*10, units = "cm", res = 200)
MultiFeaturePlot(seurat_data, plot_stage = TRUE, stage_col = "stage", gene_list = plac_nc_genes, n_col = ncol, label = '')
graphics.off()


# add neural crest and placodal cell types to metadata
plac_nc_clusters <- c("Delaminating NC" = 14, "NC Progenitors" = 12, "Epi/Plac Progenitors" = 0, "Placodes" = 16, "Epidermis" = 15)

seurat_data@meta.data$plac_nc_clusters <- unlist(lapply(seurat_data@meta.data$seurat_clusters, function(x){
  ifelse(any(plac_nc_clusters %in% x), names(plac_nc_clusters)[plac_nc_clusters %in% x], NA)
}))

# plot annotated neural crest and placodal clusters
png(paste0(plot_path, "plac_nc_clusters.png"), width = 13, height = 10, res = 200, units = "cm")
DimPlot(seurat_data, group.by = "plac_nc_clusters")  + ggtitle("Clusters") + theme(plot.title = element_text(hjust = 0.5))
graphics.off()

# plot dotplot for neural crest and placodal genes and clusters
seurat_data@meta.data$plac_nc_clusters <- factor(seurat_data@meta.data$plac_nc_clusters, levels = rev(c("Epi/Plac Progenitors", "Epidermis", "Placodes", "NC Progenitors", "Delaminating NC")))

png(paste0(plot_path, "plac_nc_dotplot.png"), width = 25, height = 10, res = 200, units = "cm")
DotPlot(seurat_data[, !is.na(seurat_data@meta.data$plac_nc_clusters)], group.by = "plac_nc_clusters", features = rev(plac_nc_genes)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
graphics.off()


###################### Nerual crest and placode cell type identification ###################### 
np_genes <- c(
  # HB
  "GBX2", "HOXA2", "HOXA3", "HOXB2", "KROX20",
  # MHB
  "WNT4", "PAX2", "FGF8", "WNT1",
  # Anterior
  "PAX6", "OTX2", "OLIG2" , "SIX3",
  # Ventral floor
  "SHH", "NKX2-2", "FOXA2"
)

ncol = ceiling((length(np_genes)+1)/8)+1
nrow = ceiling((length(np_genes)+1)/ncol)

# plot expression of np genes
png(paste0(plot_path, 'multi_feature_np.png'), width = ncol*10, height = nrow*10, units = "cm", res = 200)
MultiFeaturePlot(seurat_data, plot_stage = TRUE, stage_col = "stage", gene_list = np_genes, n_col = ncol, label = '')
graphics.off()


# add neural crest and placodal cell types to metadata
np_clusters <- c("Hindbrain" = 3, "Midbrain" = 7, "Forebrain" = 9, "Ventral Forebrain" = 5) 

seurat_data@meta.data$np_clusters <- unlist(lapply(seurat_data@meta.data$seurat_clusters, function(x){
  ifelse(any(np_clusters %in% x), names(np_clusters)[np_clusters %in% x], NA)
}))

# plot annotated np clusters
png(paste0(plot_path, "np_clusters.png"), width = 13, height = 10, res = 200, units = "cm")
DimPlot(seurat_data, group.by = "np_clusters")  + ggtitle("Clusters") + theme(plot.title = element_text(hjust = 0.5))
graphics.off()

# plot dotplot for NP genes and clusters
seurat_data@meta.data$np_clusters <- factor(seurat_data@meta.data$np_clusters, levels = c("Hindbrain", "Midbrain", "Forebrain", "Ventral Forebrain"))

png(paste0(plot_path, "np_dotplot.png"), width = 25, height = 10, res = 200, units = "cm")
DotPlot(seurat_data[, !is.na(seurat_data@meta.data$np_clusters)], group.by = "np_clusters", features = rev(np_genes)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
graphics.off()
