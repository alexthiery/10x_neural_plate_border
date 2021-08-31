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

# ENSGALG00000030902 == SNAIL2

cell_type_markers <- list(delaminating_nc = c("ETS1", "LMO4", "SOX10", "SOX8", "FOXD3"),
                          nc = c("ETS1", "ENSGALG00000030902", "FOXD3", "TFAP2B", "TFAP2A"),
                          p_npb = c("PAX7", "MSX1", "GBX2", "SIX1", "DLX5", "DLX6", "TFAP2B", "TFAP2A", "TFAP2C"), # check ZIC1 expression
                          a_npb = c("SIX3", "PAX6", "OTX2", "SIX1", "DLX5", "DLX6", "TFAP2B", "TFAP2A", "TFAP2C"),
                          npb = c("PAX7", "SIX1", "DLX5", "DLX6", "TFAP2B", "TFAP2A", "TFAP2C"),
                          placodes = c("SIX1", "EYA1", "EYA2", "DLX3", "DLX5", "DLX6", "GATA2", "GATA3"),
                          placodal_progenitors = c("DLX5", "DLX6", "GATA2", "GATA3"),
                          # PPR = c("SIX1", "EYA2", "DLX5", "DLX6"),
                          # aPPR = c("SIX1", "EYA2", "DLX5", "DLX6", "SIX3", "OTX2"),
                          # pPPR = c("SIX1", "EYA2", "DLX5", "DLX6", "GBX2"),
                          non_neural = c("EPAS1", "MSX2", "GATA2", "GATA3", "VGLL1", "GRHL3"),
                          hindbrain = c("GBX2", "HOXA2", "HOXA3", "HOXB2", "KROX20", "SOX2", "SOX21"),
                          midbrain = c("WNT4", "PAX2", "FGF8", "WNT1", "OTX2", "SOX2", "SOX21"),
                          forebrain = c("PAX6" , "SIX3", "OTX2", "SOX2", "SOX21"),
                          early_hindbrain = c("GBX2", "SOX2", "SOX3", "FRZB"),
                          early_midbrain = c("WNT4", "OTX2", "SOX2", "SOX3", "FRZB"),
                          early_forebrain = c("PAX6" , "SIX3", "OTX2", "SOX2", "SOX3", "FRZB"),
                          neural_progrenitors = c("SOX3", "SOX2", "FRZB", "LMO1"))

# hh5/6
neural = c("SOX2", "SOX3")
border = c("SOX2", "SOX3", "DLX5", "DLX6", "GATA2", "GATA3")
non_neural = c("DLX5", "DLX6", "GATA2", "GATA3")

# hh6/7
non_neural = c("EPAS1", "MSX2", "GATA2", "GATA3", "GRHL3", "HIF2A", "MSX2")#krt19 keith mclaren 2003
npb = c("SIX1", "EYA2", "DLX5", "DLX6", "TFAP2B", "TFAP2A", "TFAP2C", "PRDM1", "SOX2", "SOX3", "MSX2")
p_npb = c("PAX7", "MSX1", "GBX2", "SIX1", "EYA2", "DLX5", "DLX6", "TFAP2B", "TFAP2A", "TFAP2C", "PRDM1", "SOX2", "SOX3", "MSX2") # check ZIC1 expression
a_npb = c("SIX3", "PAX6", "OTX2", "SIX1", "EYA2", "DLX5", "DLX6", "TFAP2B", "TFAP2A", "TFAP2C", "PRDM1", "SOX2", "SOX3", "MSX2")
aPPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "PRDM1", "SIX3", "PAX6", "HESX1", "OTX2", "TFAP2A") # PNOC, SSTR5
pPPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "PRDM1", "GBX2", "FOXI3", "TFAP2A")
PPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "PRDM1", "TFAP2A") # TFAPs

neural_progenitors = c("SOX2", "SOX21", "LMO1", "SIP1", "SOX3", "FRZB")
a_neural_progenitors = c("OTX2", "SIX3", "HESX1", "SOX2", "SOX21", "LMO1", "SIP1", "SOX3", "FRZB")
p_neural_progenitors = c("GBX2", "SOX2", "SOX21", "LMO1", "SIP1", "SOX3", "FRZB")
pan_neural = c("SOX2", "SOX21", "LMO1", "SIP1", "SOX1", "SOX3", "FRZB")
hindbrain = c("GBX2", "HOXA2", "HOXA3", "HOXB2", "KROX20", "SOX2", "SOX21", "LMO1", "SIP1", "SOX1", "SOX3", "FRZB") #ZEB2
midbrain = c("WNT4", "PAX2", "FGF8", "WNT1", "OTX2", "SOX2", "SOX21", "LMO1", "SIP1", "SOX1", "SOX3", "FRZB")
forebrain = c("PAX6" , "SIX3", "OTX2", "SOX2", "SOX21", "LMO1", "SIP1", "SOX1", "SOX3", "FRZB") #TLX1

# ss4
early_NC = c("ETS1", "ENSGALG00000030902", "FOXD3", "TFAP2B", "TFAP2A")
aPPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "SIX3", "PAX6", "HESX1", "OTX2") #TFAP2?
pPPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "GBX2", "PAX2", "SOX8", "FOXI3")
iPPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "PAX3")
hindbrain = c("GBX2", "HOXA2", "HOXA3", "HOXB2", "KROX20", "SOX2", "SOX21", "LMO1", "SIP1", "GLI2", "ZNF423") #ZNF423/GLI2 (Trevers 2021) #ZEB2
midbrain = c("WNT4", "PAX2", "FGF8", "WNT1", "OTX2", "SOX2", "SOX21", "LMO1", "SIP1", "GLI2", "ZNF423")
forebrain = c("PAX6" , "SIX3", "OTX2", "SOX2", "SOX21", "LMO1", "SIP1", "GLI2", "ZNF423") #TLX1

# ss8
late_NC = c("ETS1", "LMO4", "SOX10", "SOX8", "FOXD3")






# ventral_floor = c("SHH", "NKX2-2", "FOXA2"),
# neural_progenitors = c("SOX21", "FRZB", "LMO1"))
# anterior_neural_progenitor = c("SOX2", "SOX21", "FRZB", "LMO1", "PAX6", "OTX2", "SIX3"),
# posterior_neural_progenitor = c("SOX2", "SOX21", "FRZB", "LMO1", "GBX2"))


seurat_data <- ClusterClassification(seurat_obj = seurat_data, cell_type_markers = cell_type_markers, quantile = 0.8, plot_path = paste0(plot_path, "scHelper_log/"))

# Plot UMAP for clusters and developmental stage
png(paste0(plot_path, "scHelper_celltype_umap.png"), width=50, height=20, units = 'cm', res = 200)
ClustStagePlot(seurat_data, stage_col = "stage", cluster_col = "scHelper_cell_type")
graphics.off()

saveRDS(seurat_data, paste0(rds_path, "cell_state_classification.RDS"), compress = FALSE)

# 
# 
# genes <- c("SOX10", "DLX2", "SIX1", "TEAD3", "TEAD4", "SOX21", "SOX2")
# 
# for(i in genes){
#   png(paste0(plot_path, "temp_", i, ".png"), width=25, height=20, units = 'cm', res = 200)
#   print(Seurat::VlnPlot(seurat_data, i, group.by = 'scHelper_cell_type'))
#   graphics.off()
# }
# 



# 
# # Plot multi feature plot
# nc_genes = cell_type_markers[c("early_nc", "delaminating_nc")]
# 
# ncol = ceiling((length(unlist(nc_genes))+1)/8)+1
# nrow = ceiling((length(unlist(nc_genes))+1)/ncol)
# 
# # plot expression of NC and placodal genes
# png(paste0(plot_path, 'multi_feature_nc.png'), width = ncol*10, height = nrow*10, units = "cm", res = 200)
# MultiFeaturePlot(seurat_data, plot_stage = TRUE, stage_col = "stage", gene_list = unlist(nc_genes), n_col = ncol, label = '')
# graphics.off()
# 
# 
# # plot annotated neural crest and placodal clusters
# png(paste0(plot_path, "plac_nc_clusters.png"), width = 13, height = 10, res = 200, units = "cm")
# DimPlot(seurat_data, group.by = "plac_nc_clusters")  + ggtitle("Clusters") + theme(plot.title = element_text(hjust = 0.5))
# graphics.off()
# 
# # plot dotplot for neural crest and placodal genes and clusters
# seurat_data@meta.data$plac_nc_clusters <- factor(seurat_data@meta.data$plac_nc_clusters, levels = rev(c("Epi/Plac Progenitors", "Epidermis", "Placodes", "NC Progenitors", "Delaminating NC")))
# 
# png(paste0(plot_path, "plac_nc_dotplot.png"), width = 25, height = 10, res = 200, units = "cm")
# DotPlot(seurat_data[, !is.na(seurat_data@meta.data$plac_nc_clusters)], group.by = "plac_nc_clusters", features = rev(plac_nc_genes)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# graphics.off()
# 
# 
# ###################### Neural crest and placode cell type identification ###################### 
# np_genes <- c(
#   # HB
#   "GBX2", "HOXA2", "HOXA3", "HOXB2", "KROX20",
#   # MHB
#   "WNT4", "PAX2", "FGF8", "WNT1",
#   # Anterior
#   "PAX6", "OTX2", "OLIG2" , "SIX3",
#   # Ventral floor
#   "SHH", "NKX2-2", "FOXA2"
# )
# 
# ncol = ceiling((length(np_genes)+1)/8)+1
# nrow = ceiling((length(np_genes)+1)/ncol)
# 
# # plot expression of np genes
# png(paste0(plot_path, 'multi_feature_np.png'), width = ncol*10, height = nrow*10, units = "cm", res = 200)
# MultiFeaturePlot(seurat_data, plot_stage = TRUE, stage_col = "stage", gene_list = np_genes, n_col = ncol, label = '')
# graphics.off()
# 
# # add neural crest and placodal cell types to metadata
# np_clusters <- c("Hindbrain" = 3, "Midbrain" = 7, "Forebrain" = 9, "Ventral Forebrain" = 5) 
# 
# seurat_data@meta.data$np_clusters <- unlist(lapply(seurat_data@meta.data$seurat_clusters, function(x){
#   ifelse(any(np_clusters %in% x), names(np_clusters)[np_clusters %in% x], NA)
# }))
# 
# # plot annotated np clusters
# png(paste0(plot_path, "np_clusters.png"), width = 13, height = 10, res = 200, units = "cm")
# DimPlot(seurat_data, group.by = "np_clusters")  + ggtitle("Clusters") + theme(plot.title = element_text(hjust = 0.5))
# graphics.off()
# 
# # plot dotplot for NP genes and clusters
# seurat_data@meta.data$np_clusters <- factor(seurat_data@meta.data$np_clusters, levels = c("Hindbrain", "Midbrain", "Forebrain", "Ventral Forebrain"))
# 
# png(paste0(plot_path, "np_dotplot.png"), width = 25, height = 10, res = 200, units = "cm")
# DotPlot(seurat_data[, !is.na(seurat_data@meta.data$np_clusters)], group.by = "np_clusters", features = rev(np_genes)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# graphics.off()
# 
# 
# ###################### progenitor cell type identification ###################### 
# progenitor_genes <- c("KRT7", "KRT18", "GATA2", "HOMER2", "DLX6", "DLX5", "EYA2", "SIX1", "ZNF385C", "TFAP2A", "MSX1", "CSRNP1", "PAX7", "DRAXIN", "SOX21", "SHISA2", "FRZB", "LMO1", "PAX6")
# 
# ncol = ceiling((length(progenitor_genes)+1)/8)+1
# nrow = ceiling((length(progenitor_genes)+1)/ncol)
# 
# # plot expression of np genes
# png(paste0(plot_path, 'multi_feature_progenitors.png'), width = ncol*10, height = nrow*10, units = "cm", res = 200)
# MultiFeaturePlot(seurat_data, plot_stage = TRUE, stage_col = "stage", gene_list = progenitor_genes, n_col = ncol, label = '')
# graphics.off()
# 
# 
# # add progenitor cell types to metadata
# progenitor_clusters <- c("Placodal progenitors?" = 0, "Neural/NC progenitors?" = 6, "Naive ectoderm?" = 2, "Forebrain progenitors?" = 9, "Neural progenitors?" = 8)
# 
# seurat_data@meta.data$progenitor_clusters <- unlist(lapply(seurat_data@meta.data$seurat_clusters, function(x){
#   ifelse(any(progenitor_clusters %in% x), names(progenitor_clusters)[progenitor_clusters %in% x], NA)
# }))
# 
# # plot annotated np clusters
# png(paste0(plot_path, "progenitor_clusters.png"), width = 13, height = 10, res = 200, units = "cm")
# DimPlot(seurat_data, group.by = "progenitor_clusters")  + ggtitle("Clusters") + theme(plot.title = element_text(hjust = 0.5))
# graphics.off()
# 
# # plot dotplot for NP genes and clusters
# seurat_data@meta.data$progenitor_clusters <- factor(seurat_data@meta.data$progenitor_clusters, levels = c("Placodal progenitors?", "Naive ectoderm?", "Neural/NC progenitors?", "Neural progenitors?", "Forebrain progenitors?"))
# 
# png(paste0(plot_path, "progenitors_dotplot.png"), width = 25, height = 10, res = 200, units = "cm")
# DotPlot(seurat_data[, !is.na(seurat_data@meta.data$progenitor_clusters)], group.by = "progenitor_clusters", features = rev(progenitor_genes)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# graphics.off()
# 


