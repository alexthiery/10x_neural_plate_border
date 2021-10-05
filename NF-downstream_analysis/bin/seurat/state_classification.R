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
    
    plot_path = "./test/state_classification/plots/"
    rds_path = "./test/state_classification/rds_files/"
    data_path = "./output/NF-downstream_analysis_stacas/filtered_seurat/seurat/state_classification/rds_files/"
    
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

# Retrieve seurat object label
label <- sub('_.*', '', list.files(data_path))

seurat_data <- readRDS(list.files(data_path, full.names = TRUE))
# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/hh4_splitstage_data/seurat/stage_cluster/rds_files/hh4_clustered_data.RDS')

########################################################################################################
#                                      Cell state classification                                    #
#######################################################################################

cell_states = list(
  extra_embryonic = c()

  early_NNE = c()
  NNE = c()
  prospective_epidermis = c()

  PPR = c()
  aPPR = c()
  pPPR = c()

  early_NPB = c()
  NPB = c()
  aNPB = c()
  pNPB = c()

  NC = c()
  delaminating_NC = c()

  early_neural = c()
  early_caudal_neural = c()

  NP = c()
  pNP = c()
  iNP = c()
  pNP = c()
  hindbrain = c()
  midbrain = c()
  forebrain = c()
  ventral_floorplate = c()

  node = c()
  streak = c()
)

hh4_cell_type_markers = list(  node = c('EOMES', 'ADMP', 'CHRD', 'SHH', 'GNOT2', 'CNOT1', 'FOXA2'),#OTX2
                               streak = c('EOMES', 'TBX6', 'TBXT', 'GSC'),
                               early_neural = c("SOX2", "SOX3", 'OTX2', 'EPCAM', 'MAFA', 'FRZB', "YEATS4", 'SOX11', 'ERN1'), # many from trevers 2021 / katherine thesis
                               early_border = c("SOX2", "SOX3", 'OTX2', 'EPCAM', 'MAFA', 'FRZB', "YEATS4", 'SOX11', 'ERN1', "DLX5", "DLX6", "GATA2", "GATA3"),
                               early_NNE = c("DLX5", "DLX6", "GATA2", "GATA3"),
                               extra_embryonic = c('VGLL1', 'GRHL3', 'GATA2', 'GATA3'))

hh5_cell_type_markers = list(early_caudal_neural = c('GBX2', 'SP5', 'HOXB1', 'CDX2', 'SOX2', 'SOX21', 'SOX3'), #TBXT...
                             early_neural_plate = c('OTX2', 'SOX2', 'SOX21', 'SOX3', 'FRZB'),
                             early_pNPB = c("PAX7", "MSX1", "GBX2", 'SP5', "DLX5", "DLX6", "TFAP2A", "TFAP2C", "PRDM1", "SOX2", 'SOX3', "SOX21", "ERN1", 'FRZB'),
                             early_aNPB = c("SIX3", "OTX2", "DLX5", "DLX6", "TFAP2A", "TFAP2C", "PRDM1", "SOX2", 'SOX3', "SOX21", "ERN1", 'FRZB'),
                             early_NPB = c("DLX5", "DLX6", "TFAP2A", "TFAP2C", "PRDM1", "SOX2", 'SOX3', "SOX21", "ERN1", 'FRZB'),
                             NNE = c('ASTL', "DLX5", "DLX6", 'TFAP2A', "TFAP2C", "GATA2", "GATA3", "EPAS1"))

hh6_cell_type_markers = list(  prospective_epidermis = c("MSX2", "EPAS1", "GATA2", "GATA3", "GRHL3"),#krt19 keith mclaren 2003 MSX2?
                               NPB = c("MYC", "ZIC1", "SIX1", "EYA2", "DLX5", "DLX6", "TFAP2B", "TFAP2A", "TFAP2C", "PRDM1", "SOX2", "SOX3", 'SOX21', "MSX2", 'FRZB', 'PAX7'),
                               pNPB = c('ZIC1', "PAX7", "MSX1", "GBX2", "SIX1", "EYA2", "DLX5", "DLX6", "TFAP2B", "TFAP2A", "TFAP2C", "PRDM1", "SOX2", "SOX3",'SOX21', "MSX2", 'FRZB'), # check ZIC1 expression
                               aNPB = c("MYC", "ZIC1", "SIX3", "PAX6", "OTX2", "SIX1", "EYA2", "DLX5", "DLX6", "TFAP2B", "TFAP2A", "TFAP2C", "PRDM1", "SOX2", "SOX3",'SOX21', "MSX2", 'FRZB'),
                               aPPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "SIX3", "PAX6", "OTX2"), #TFAP2?
                               pPPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "GBX2", "PAX2", "SOX8"), #, "FOXI3"
                              #  iPPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "Pax3"),
                               #  early_aPPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "PRDM1", "SIX3", "PAX6", "HESX1", "OTX2", "TFAP2A"), # PNOC, SSTR5
                               #  early_pPPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "PRDM1", "GBX2", "TFAP2A"), # FOXI3
                               PPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "PRDM1", "TFAP2A"), # TFAPs
                               pan_neural = c("SOX2", "SOX21", "LMO1", "ZEB2", "SOX1", "SOX3", "FRZB"),
                              #  a_neural_progenitors = c("OTX2", "SIX3", "HESX1", "SOX2", "SOX21", "LMO1", "ZEB2", "SOX3", "FRZB"),
                              #  p_neural_progenitors = c("GBX2", "SOX2", "SOX21", "LMO1", "ZEB2", "SOX3", "FRZB"),
                               posterior_neural = c("GBX2", "HOXA2", "HOXA3", "HOXB2", "KROX20", "SOX2", "SOX21", "LMO1", "ZEB2", "SOX1", "SOX3", "FRZB", "WNT8A"),
                               intermediate_neural = c("WNT4", "PAX2", "FGF8", "WNT1", "OTX2", "SOX2", "SOX21", "LMO1", "ZEB2", "SOX1", "SOX3", "FRZB"),
                               anterior_neural = c("SIX3", "OTX2", "SOX2", "SOX21", "LMO1", "ZEB2", "SOX1", "SOX3", "FRZB", 'SHISA2', 'HESX1')) #TLL1/NR2E1) # "PAX6"

hh7_cell_type_markers = list(  NC = c("PAX7", "MSX1", "MSX2", "ETS1", "ENSGALG00000030902", "FOXD3", "TFAP2B", "TFAP2A"), # MSX2?
                               hindbrain = c("GBX2", "HOXA2", "HOXA3", "HOXB2", "KROX20", "SOX2", "SOX21", 'FRZB', "LMO1", "ZEB2", "GLI2"), #ZNF423/GLI2 (Trevers 2021) #ZEB2
                               midbrain = c("WNT4", "PAX2", "FGF8", "WNT1", "OTX2", "SOX2", "SOX21", 'FRZB', "LMO1", "ZEB2", "GLI2"))

ss4_cell_type_markers = list( forebrain = c("MYC", "PAX6" , "SIX3", "OTX2", "SOX2", "SOX21", 'FRZB', "LMO1", "ZEB2", "GLI2", 'SHISA2', 'ZIC1', 'HESX1'))

ss8_cell_type_markers = list( delaminating_NC = c("ETS1", "LMO4", "SOX10", "SOX8", "FOXD3"),
                              ventral_floorplate = c('SHH', 'FOXA2', 'NKX2-2'))


cell_type_markers = list(hh4 = c(hh4_cell_type_markers, hh5_cell_type_markers),
                         hh5 = c(hh4_cell_type_markers, hh5_cell_type_markers, hh6_cell_type_markers),
                         hh6 = c(hh4_cell_type_markers, hh5_cell_type_markers, hh6_cell_type_markers),
                         hh7 = c(hh6_cell_type_markers, hh7_cell_type_markers, ss4_cell_type_markers, ss8_cell_type_markers),
                         ss4 = c(hh6_cell_type_markers, hh7_cell_type_markers, ss4_cell_type_markers, ss8_cell_type_markers),
                         ss8 = c(hh6_cell_type_markers, hh7_cell_type_markers, ss4_cell_type_markers, ss8_cell_type_markers))

# Run classification using different resolutions for different stages
stage = unique(seurat_data@meta.data$stage)

if(length(stage) == 1){
  cell_type_markers = cell_type_markers[[stage]]
  cluster_res = list(hh4 = 1, hh5 = 1, hh6 = 1, hh7 = 1.2, ss4 = 1.2, ss8 = 1.2)[[stage]]
} else {
  cell_type_markers = flatten(cell_type_markers)
  cell_type_markers = cell_type_markers[!duplicated(cell_type_markers)]
  
  # Set cluster res to 2 for run split subsets - 3 for integrated data
  cluster_res = ifelse(length(unique(seurat_data$run)) == 1, 2, 3)
}

DefaultAssay(seurat_data) <- "integrated"

seurat_data <- FindClusters(seurat_data, resolution = cluster_res)

# Set RNA to default assay for plotting expression data
DefaultAssay(seurat_data) <- "RNA"

seurat_data <- ClusterClassification(seurat_obj = seurat_data, cell_type_markers = cell_type_markers, force_assign = FALSE, quantile = 0.5, plot_path = paste0(plot_path, "scHelper_log/"))

# Plot UMAP for clusters and developmental stage
png(paste0(plot_path, "scHelper_celltype_umap.png"), width=40, height=20, units = 'cm', res = 200)
ClustStagePlot(seurat_data, stage_col = "stage", cluster_col = "scHelper_cell_type", label_clusters = TRUE)
graphics.off()

png(paste0(plot_path, "scHelper_celltype_umap2.png"), width=20, height=20, units = 'cm', res = 200)
DimPlot(seurat_data, group.by = 'scHelper_cell_type', label = TRUE, label.size = 3, label.box = TRUE, repel = TRUE) + ggplot2::theme(legend.position = "none")
graphics.off()

saveRDS(seurat_data, paste0(rds_path, label, "_cell_state_classification.RDS"), compress = FALSE)

# Plot stacked violins for each of the cell type classes to check genes used are good markers
curr_plot_path = paste0(plot_path, "cell_type_dotplots/")
dir.create(curr_plot_path)

for(i in names(cell_type_markers)){
  png(paste0(curr_plot_path, i, ".png"), width = (length(cell_type_markers[[i]])+2)*3, height = 15, units = 'cm', res = 200)
  print(DotPlot(seurat_data, features = cell_type_markers[[i]], group.by = "scHelper_cell_type"))
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