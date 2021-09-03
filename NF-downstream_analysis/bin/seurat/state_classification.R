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
# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/filtered_seurat/seurat/state_classification/rds_files/cell_state_classification.RDS')

# Set RNA to default assay for plotting expression data
DefaultAssay(seurat_data) <- "RNA"

########################################################################################################
#                                      Cell state classification                                    #
########################################################################################################



# cell_type_markers <- list(delaminating_nc = c("ETS1", "LMO4", "SOX10", "SOX8", "FOXD3"),
#                           nc = c("ETS1", "ENSGALG00000030902", "FOXD3", "TFAP2B", "TFAP2A"),
#                           p_npb = c("PAX7", "MSX1", "GBX2", "SIX1", "DLX5", "DLX6", "TFAP2B", "TFAP2A", "TFAP2C"), # check ZIC1 expression
#                           a_npb = c("SIX3", "PAX6", "OTX2", "SIX1", "DLX5", "DLX6", "TFAP2B", "TFAP2A", "TFAP2C"),
#                           npb = c("PAX7", "SIX1", "DLX5", "DLX6", "TFAP2B", "TFAP2A", "TFAP2C"),
#                           placodes = c("SIX1", "EYA1", "EYA2", "DLX3", "DLX5", "DLX6", "GATA2", "GATA3"),
#                           placodal_progenitors = c("DLX5", "DLX6", "GATA2", "GATA3"),
#                           # PPR = c("SIX1", "EYA2", "DLX5", "DLX6"),
#                           # aPPR = c("SIX1", "EYA2", "DLX5", "DLX6", "SIX3", "OTX2"),
#                           # pPPR = c("SIX1", "EYA2", "DLX5", "DLX6", "GBX2"),
#                           non_neural = c("EPAS1", "MSX2", "GATA2", "GATA3", "VGLL1", "GRHL3"),
#                           hindbrain = c("GBX2", "HOXA2", "HOXA3", "HOXB2", "KROX20", "SOX2", "SOX21"),
#                           midbrain = c("WNT4", "PAX2", "FGF8", "WNT1", "OTX2", "SOX2", "SOX21"),
#                           forebrain = c("PAX6" , "SIX3", "OTX2", "SOX2", "SOX21"),
#                           early_hindbrain = c("GBX2", "SOX2", "SOX3", "FRZB"),
#                           early_midbrain = c("WNT4", "OTX2", "SOX2", "SOX3", "FRZB"),
#                           early_forebrain = c("PAX6" , "SIX3", "OTX2", "SOX2", "SOX3", "FRZB"),
#                           neural_progrenitors = c("SOX3", "SOX2", "FRZB", "LMO1"))

# ventral_floor = c("SHH", "NKX2-2", "FOXA2"),


# ENSGALG00000030902 == SNAIL2
# HIF2A == EPAS1
# SIP1 == ZEB2
# FOXI3 not present

cell_type_markers <- list(
  # hh5/6
  early_neural = c("SOX2", "SOX3"),
  early_border = c("SOX2", "SOX3", "DLX5", "DLX6", "GATA2", "GATA3"),
  early_non_neural = c("DLX5", "DLX6", "GATA2", "GATA3"),
  
  # hh6/7
  non_neural = c("EPAS1", "MSX2", "GATA2", "GATA3", "GRHL3", "EPAS1", "MSX2"),#krt19 keith mclaren 2003
  NPB = c("SIX1", "EYA2", "DLX5", "DLX6", "TFAP2B", "TFAP2A", "TFAP2C", "PRDM1", "SOX2", "SOX3", "MSX2"),
  pNPB = c("PAX7", "MSX1", "GBX2", "SIX1", "EYA2", "DLX5", "DLX6", "TFAP2B", "TFAP2A", "TFAP2C", "PRDM1", "SOX2", "SOX3", "MSX2"), # check ZIC1 expression
  aNPB = c("SIX3", "PAX6", "OTX2", "SIX1", "EYA2", "DLX5", "DLX6", "TFAP2B", "TFAP2A", "TFAP2C", "PRDM1", "SOX2", "SOX3", "MSX2"),
  early_aPPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "PRDM1", "SIX3", "PAX6", "HESX1", "OTX2", "TFAP2A"), # PNOC, SSTR5
  early_pPPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "PRDM1", "GBX2", "TFAP2A"), # FOXI3
  early_PPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "PRDM1", "TFAP2A"), # TFAPs
  neural_progenitors = c("SOX2", "SOX21", "LMO1", "ZEB2", "SOX1", "SOX3", "FRZB"),
  a_neural_progenitors = c("OTX2", "SIX3", "HESX1", "SOX2", "SOX21", "LMO1", "ZEB2", "SOX3", "FRZB"),
  p_neural_progenitors = c("GBX2", "SOX2", "SOX21", "LMO1", "ZEB2", "SOX3", "FRZB"),
  early_hindbrain = c("GBX2", "HOXA2", "HOXA3", "HOXB2", "KROX20", "SOX2", "SOX21", "LMO1", "ZEB2", "SOX1", "SOX3", "FRZB"),
  early_midbrain = c("WNT4", "PAX2", "FGF8", "WNT1", "OTX2", "SOX2", "SOX21", "LMO1", "ZEB2", "SOX1", "SOX3", "FRZB"),
  early_forebrain = c("PAX6" , "SIX3", "OTX2", "SOX2", "SOX21", "LMO1", "ZEB2", "SOX1", "SOX3", "FRZB"), #TLX1
  
  # ss4
  NC = c("ETS1", "ENSGALG00000030902", "FOXD3", "TFAP2B", "TFAP2A"),
  aPPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "SIX3", "PAX6", "HESX1", "OTX2"), #TFAP2?
  pPPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "GBX2", "PAX2", "SOX8"), #FOXI3
  iPPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "Pax3"),
  hindbrain = c("GBX2", "HOXA2", "HOXA3", "HOXB2", "KROX20", "SOX2", "SOX21", "LMO1", "ZEB2", "GLI2", "ZNF423"), #ZNF423/GLI2 (Trevers 2021) #ZEB2
  midbrain = c("WNT4", "PAX2", "FGF8", "WNT1", "OTX2", "SOX2", "SOX21", "LMO1", "ZEB2", "GLI2", "ZNF423"),
  forebrain = c("PAX6" , "SIX3", "OTX2", "SOX2", "SOX21", "LMO1", "ZEB2", "GLI2", "ZNF423"), #TLX1
  
  # ss8
  delaminating_NC = c("ETS1", "LMO4", "SOX10", "SOX8", "FOXD3")
)

seurat_data <- ClusterClassification(seurat_obj = seurat_data, cell_type_markers = cell_type_markers, quantile = 0.8, force_assign = TRUE, plot_path = paste0(plot_path, "scHelper_log/"))

# Plot UMAP for clusters and developmental stage
png(paste0(plot_path, "scHelper_celltype_force_umap.png"), width=50, height=20, units = 'cm', res = 200)
ClustStagePlot(seurat_data, stage_col = "stage", cluster_col = "scHelper_cell_type", label_clusters = TRUE)
graphics.off()

seurat_data <- ClusterClassification(seurat_obj = seurat_data, cell_type_markers = cell_type_markers, quantile = 0.8, plot_path = paste0(plot_path, "scHelper_log/"))

# Plot UMAP for clusters and developmental stage
png(paste0(plot_path, "scHelper_celltype_umap.png"), width=50, height=20, units = 'cm', res = 200)
ClustStagePlot(seurat_data, stage_col = "stage", cluster_col = "scHelper_cell_type", label_clusters = TRUE)
graphics.off()

saveRDS(seurat_data, paste0(rds_path, "cell_state_classification.RDS"), compress = FALSE)

# Plot stacked violins for each of the cell type classes to check genes used are good markers
curr_plot_path = paste0(plot_path, "cell_type_violins/")
dir.create(curr_plot_path)

for(i in names(cell_type_markers)){
  png(paste0(curr_plot_path, i, ".png"), width = (length(cell_type_markers[[i]])+2)*3, height = 15, units = 'cm', res = 200)
  print(VlnPlot(seurat_data, cell_type_markers[[i]], group.by = "scHelper_cell_type", stack = TRUE))
  graphics.off()
}


# Plot multi-feature plots for each of the cell type classes
curr_plot_path = paste0(plot_path, "cell_type_feature_plots/")
dir.create(curr_plot_path)

for(i in names(cell_type_markers)){
  ncol = ceiling((length(cell_type_markers[[i]])+1)/8)+1
  nrow = ceiling((length(cell_type_markers[[i]])+1)/ncol)
  
  png(paste0(curr_plot_path, i, '.png'), width = ncol*10, height = nrow*10, units = "cm", res = 200)
  MultiFeaturePlot(seurat_data, plot_stage = TRUE, stage_col = "stage", plot_celltype = TRUE, celltype_col = "scHelper_cell_type",
                   gene_list = cell_type_markers[[i]], n_col = ncol, label = '')
  graphics.off()
}