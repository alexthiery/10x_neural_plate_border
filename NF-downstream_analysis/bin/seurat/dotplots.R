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
    data_path = "./input/"
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
label <- sub('_.*', '', list.files(data_path, pattern = '*.RDS'))

seurat_data <- readRDS(list.files(data_path, full.names = TRUE, pattern = '*.RDS'))
hh4_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/hh4_splitstage_data/seurat/stage_state_classification/rds_files/hh4_cell_state_classification.RDS')
hh5_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/hh5_splitstage_data/seurat/stage_state_classification/rds_files/hh5_cell_state_classification.RDS')
# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/hh6_splitstage_data/seurat/stage_cluster/rds_files/hh6_clustered_data.RDS')
# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/hh7_splitstage_data/seurat/stage_cluster/rds_files/hh7_clustered_data.RDS')
# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/ss4_splitstage_data/seurat/stage_cluster/rds_files/ss4_clustered_data.RDS')
# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/ss8_splitstage_data/seurat/stage_cluster/rds_files/ss8_clustered_data.RDS')


DimPlot(seurat_data, group.by = 'scHelper_cell_type')

meta_col = 'scHelper_cell_type'
if(meta_col == 'scHelper_cell_type'){
  scHelper_cell_type_order <- c('extra_embryonic', 'NNE', 'prospective_epidermis', 'PPR', 'aPPR', 'pPPR', 'early_NPB', 'NPB',
                                'aNPB', 'pNPB', 'NC', 'delaminating_NC', 'early_neural', 'early_caudal_neural', 'NP', 'pNP',
                                'hindbrain', 'iNP', 'midbrain', 'aNP', 'forebrain', 'ventral_forebrain', 'node', 'streak')
  
  scHelper_cell_type_order <- scHelper_cell_type_order[scHelper_cell_type_order %in% unique(seurat_data@meta.data[[meta_col]])]
  
  if(!all(seurat_data@meta.data$scHelper_cell_type %in% scHelper_cell_type_order)){
    stop('Check scHelper_cell_type_order. Cell types found within "seurat_data@meta.data$scHelper_cell_type" are missing from custom order vector')
  }
  
  seurat_data@meta.data$scHelper_cell_type <- factor(seurat_data@meta.data$scHelper_cell_type, levels = scHelper_cell_type_order)
}

Idents(seurat_data) <- meta_col


data_subset <- subset(seurat_data, cells = rownames(filter(seurat_data@meta.data, scHelper_cell_type %in% c('NNE', 'early_NPB', 'early_neural'))))

data_subset <- subset(seurat_data, cells = rownames(filter(seurat_data@meta.data, scHelper_cell_type %in% c('NNE', 'aPPR', 'pNPB', 'pNP', 'iNP', 'aNP'))))


goi <- list(hh4 = c('DLX5', 'GATA2', 'MSX1', 'BMP4', 'TFAP2A', 'TFAP2C', 'PAX7', 'SOX3', 'OTX2', 'MAFA', 'FRZB'),
            hh5 = c('SOX21', 'SIX3', 'OTX2', 'GBX2', 'BTG2', 'MSX1', 'PAX7', 'TFAP2B', 'SIX1', 'EYA2', 'TFAP2C', 'DLX5', 'BMP4'))


DotPlot(hh4_data, features = genes) + 
  theme(axis.text.x = element_text(angle = 90)) + coord_flip()




seurat_data@misc <- data.frame(temp = 1:10)


seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/transfer_labels/seurat/rds_files/seurat_label_transfer.RDS')





dotplot <- function(data, cells, genes, group_by, facet, group_by_levels=NULL, facet_dir = 'h'){
  
  if(is.null(group_by_levels)){
    group_by_levels <- unique(seurat_data@meta.data[cells, group_by])
  }
  # Extract normalised data for genes and cells of interest and scale
  scale_data <- t(GetAssayData(data, assay = 'RNA', slot = 'data')[genes, cells])
  
  # cbind facet and group data with expression data
  plot_data <- cbind(as.data.frame(scale_data), data@meta.data[rownames(scale_data),c(facet, group_by)])
  
  temp <- plot_data %>%
    pivot_longer(!c(!!sym(facet), !!sym(group_by)), names_to = 'gene', values_to = 'normalised_expression') %>%
    mutate(gene = factor(gene, levels = genes)) %>%
    mutate(!!sym(group_by) := factor(!!sym(group_by), levels = group_by_levels)) %>%
    group_by(gene) %>%
    mutate(scaled_expression = scale(normalised_expression)) %>%
    dplyr::group_by(!!sym(facet), !!sym(group_by), gene) %>%
    mutate(percent_expressed = 100*(sum(normalised_expression > 0, na.rm = TRUE)/n())) %>%
    mutate(average_expression = mean(scaled_expression)) %>%
    ungroup() %>%
    dplyr::select(c(!!sym(facet), !!sym(group_by), gene, percent_expressed, average_expression)) %>%
    distinct()
  
  temp[temp$average_expression > 1, 'average_expression'] <- 1
  
  
  plot <- ggplot(temp, aes(x = !!sym(group_by), y = gene, colour = average_expression)) +
    geom_point(aes(size = percent_expressed)) +
    scale_radius(limits = c(0, NA), range = c(0, 8)) +
    theme_bw() +
    theme(panel.grid.major = element_blank()) +
    scale_colour_gradient(low = 'gray90', high = 'blue') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  if(facet_dir == 'h'){
    plot <- plot + facet_wrap(~ stage, scales = "free_x", dir = 'h', nrow = 1)
  } else if(facet_dir == 'v'){
    plot <- plot + facet_wrap(~ stage, dir = 'v', ncol = 1)
  }
    
  return(plot)
}





genes <- c('SOX21', 'SIX3', 'OTX2', 'GBX2', 'BTG2', 'MSX1', 'PAX7', 'TFAP2B', 'SIX1', 'EYA2', 'TFAP2C', 'DLX5', 'BMP4')
cell_type_order <- c('extra_embryonic', 'NNE', 'prospective_epidermis', 'PPR', 'aPPR', 'pPPR', 'early_NPB', 'NPB',
                     'aNPB', 'pNPB', 'NC', 'delaminating_NC', 'early_neural', 'early_caudal_neural', 'NP', 'pNP',
                     'hindbrain', 'iNP', 'midbrain', 'aNP', 'forebrain', 'ventral_forebrain', 'node', 'streak')



cells <- rownames(filter(seurat_data@meta.data, stage %in% c('hh4', 'hh5', 'hh6')))

dotplot(seurat_data, cells, genes, group_by = 'scHelper_cell_type', facet = 'stage', group_by_levels = cell_type_order)


cells <- rownames(filter(seurat_data@meta.data,scHelper_cell_type %in% c('NP', 'pNP','hindbrain', 'iNP', 'midbrain', 'aNP', 'forebrain', 'early_caudal_neural')))
genes <- c('SIX3', 'OTX2', 'GBX2')

dotplot(seurat_data, cells=cells, genes, group_by = 'scHelper_cell_type', facet = 'stage', group_by_levels = cell_type_order, facet_dir = 'h')








