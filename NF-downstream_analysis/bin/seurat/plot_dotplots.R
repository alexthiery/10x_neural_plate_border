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
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
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

seurat_data <- readRDS(list.files(data_path, full.names = TRUE, pattern = '*.RDS'))

dotplot <- function(data, cells, genes, group_by, facet, group_by_levels=NULL, facet_dir = 'h', limits = c(0, NA), range = c(0,8)){
  
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
    scale_radius(limits = limits, range = range) +
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


cell_type_order <- c('extra_embryonic', 'NNE', 'prospective_epidermis', 'PPR', 'aPPR', 'pPPR', 'early_NPB', 'NPB',
                     'aNPB', 'pNPB', 'NC', 'delaminating_NC', 'early_neural', 'early_caudal_neural', 'NP', 'pNP',
                     'hindbrain', 'iNP', 'midbrain', 'aNP', 'forebrain', 'ventral_forebrain', 'node', 'streak')


genes <- rev(c('EPAS1', 'GATA3', 'SIX1', 'EYA2', 'DLX5', 'BMP4', 'MSX1', 'TFAP2A', 'TFAP2B', 'PAX7', 'SOX2', 'OTX2', 'YEATS4', 'SOX21', 'GBX2', 'SIX3', 'ADMP', 'EOMES'))
cells <- rownames(filter(seurat_data@meta.data, stage %in% c('hh4', 'hh5', 'hh6')))

png(paste0(plot_path, 'hh4-hh6_dotplot.png'), width = 23, height = 13, units = 'cm', res = 400)
dotplot(seurat_data, cells, genes, group_by = 'scHelper_cell_type', facet = 'stage', group_by_levels = cell_type_order, limits = c(5, NA), range = c(0,6)) +
  guides(size=guide_legend(title="Percent cells expressing"), colour = guide_colourbar(title="Average expression")) +
  theme(strip.text.x = element_text(size=15),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=11))
graphics.off()


genes <- rev(c('GATA3', 'DLX5', 'SIX1', 'EYA2', 'MSX1', 'TFAP2A', 'TFAP2B', 'Pax3', 'PAX7', 'CSRNP1', 'ENSGALG00000030902', 'SOX10', 'SOX2', 'SOX21', 'GBX2', 'PAX2', 'SIX3'))
cells <- rownames(filter(seurat_data@meta.data, stage %in% c('hh7', 'ss4', 'ss8')))

png(paste0(plot_path, 'hh7-ss8_dotplot.png'), width = 23, height = 13, units = 'cm', res = 400)
dotplot(seurat_data, cells, genes, group_by = 'scHelper_cell_type', facet = 'stage', group_by_levels = cell_type_order, limits = c(5, NA), range = c(0,6)) +
  guides(size=guide_legend(title="Percent cells expressing"), colour = guide_colourbar(title="Average expression")) +
  theme(strip.text.x = element_text(size=15),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=11))
graphics.off()
