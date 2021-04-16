#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
library(Seurat)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library(scHelper)

library(phateR)
library(readr)
library(viridis)
library(Rmagic)
library(Seurat)



spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')

    plot_path = "./output/NF-downstream_analysis/1_gms_all_cells/plots/"
    rds_path = "./output/NF-downstream_analysis/1_gms_all_cells/rds_files/"
    antler_path = "./output/NF-downstream_analysis/1_gms_all_cells/antler_data/"
    data_path = "./output/NF-downstream_analysis/6_contamination_filt/rds_files/"
    
    ncores = 8

  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')

    plot_path = "./plots/"
    rds_path = "./rds_files/"
    antler_path = "./antler_data/"
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
  dir.create(antler_path, recursive = T)
}

seurat_data <- readRDS(paste0(data_path, 'contamination_filt_data.RDS'))

seurat_data <- DietSeurat(seurat_data, counts = TRUE, data = TRUE, scale.data = TRUE, assays = 'RNA')


phate_dat <- t(GetAssayData(seurat_data, slot = "counts"))


phate_dat <- library.size.normalize(phate_dat)
phate_dat <- sqrt(phate_dat)
phate_out <- phate(phate_dat)


png("phateR.output.pdf", width = 10, height = 10)
ggplot(phate_out, aes(PHATE1, PHATE2, color=temp$celltype1)) +
  geom_point() +
  labs(colour='Cell Type')
dev.off()

