#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(scHelper)
library(phateR)
library(Rmagic)

spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    plot_path = "./output/NF-downstream_analysis/8_phate/plots/"
    data_path = "./output/NF-downstream_analysis/6_contamination_filt/rds_files/"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
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
}

seurat_data <- readRDS(paste0(data_path, 'contamination_filt_data.RDS'))

DefaultAssay(seurat_data) <- 'integrated'

seurat_sub <- subset(seurat_data, cells = rownames(filter(seurat_data@meta.data, !grepl("hh4", orig.ident))))

phate_dat <- t(GetAssayData(seurat_sub, slot = "scale.data"))

phate_out <- phate(phate_dat, knn = 10, decay = 100, t = 50)

# png("phateR.output.pdf", width = 10, height = 10)
ggplot(phate_out, aes(PHATE1, PHATE2, colour = seurat_sub$orig.ident)) +
  geom_point() +
  labs(colour='Cell Type')
# dev.off()

