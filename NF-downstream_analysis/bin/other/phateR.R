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
    
    plot_path = "./output/NF-downstream_analysis/phateR/plots/"
    data_path = "./output/NF-downstream_analysis/seurat/6_contamination_filt/rds_files/"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multiprocess", workers = ncores)
    options(future.globals.maxSize = 16* 1024^3) # 16gb
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
}

seurat_data <- readRDS(paste0(data_path, 'contamination_filt_data_uncomp.RDS'))

DefaultAssay(seurat_data) <- 'integrated'

phate_dat <- t(GetAssayData(seurat_data, slot = "scale.data"))

phate_out <- phate(phate_dat)

png(paste0(plot_path, "phate_hh4_default.png"), width=20, height=20, units = 'cm', res = 200)
ggplot(phate_out, aes(PHATE1, PHATE2, colour = seurat_data$orig.ident)) +
  geom_point() +
  labs(colour='Cell Type')
graphics.off()

phate_out <- phate(phate_dat, knn = 20)

png(paste0(plot_path, "phate_knn20.png"), width=20, height=20, units = 'cm', res = 200)
ggplot(phate_out, aes(PHATE1, PHATE2, colour = seurat_data$orig.ident)) +
  geom_point() +
  labs(colour='Cell Type')
graphics.off()

phate_out <- phate(phate_dat, t = 50)

png(paste0(plot_path, "phate_t50.png"), width=20, height=20, units = 'cm', res = 200)
ggplot(phate_out, aes(PHATE1, PHATE2, colour = seurat_data$orig.ident)) +
  geom_point() +
  labs(colour='Cell Type')
graphics.off()

phate_out <- phate(phate_dat, gamma = 0)

png(paste0(plot_path, "phate_gamma0.png"), width=20, height=20, units = 'cm', res = 200)
ggplot(phate_out, aes(PHATE1, PHATE2, colour = seurat_data$orig.ident)) +
  geom_point() +
  labs(colour='Cell Type')
graphics.off()

phate_out <- phate(phate_dat, knn = 20, t = 50)

png(paste0(plot_path, "phate_knn20_t50.png"), width=20, height=20, units = 'cm', res = 200)
ggplot(phate_out, aes(PHATE1, PHATE2, colour = seurat_data$orig.ident)) +
  geom_point() +
  labs(colour='Cell Type')
graphics.off()

phate_out <- phate(phate_dat, knn = 20, t = 50, gamma = 0)

png(paste0(plot_path, "phate_knn20_t50_gamma0.png"), width=20, height=20, units = 'cm', res = 200)
ggplot(phate_out, aes(PHATE1, PHATE2, colour = seurat_data$orig.ident)) +
  geom_point() +
  labs(colour='Cell Type')
graphics.off()






seurat_data <- subset(seurat_data, cells = rownames(filter(seurat_data@meta.data, !grepl("hh4", orig.ident))))

phate_dat <- t(GetAssayData(seurat_data, slot = "scale.data"))

phate_out <- phate(phate_dat)

png(paste0(plot_path, "phate_hh4_default.png"), width=20, height=20, units = 'cm', res = 200)
ggplot(phate_out, aes(PHATE1, PHATE2, colour = seurat_data$orig.ident)) +
  geom_point() +
  labs(colour='Cell Type')
graphics.off()

phate_out <- phate(phate_dat, knn = 20)

png(paste0(plot_path, "phate_hh4_knn20.png"), width=20, height=20, units = 'cm', res = 200)
ggplot(phate_out, aes(PHATE1, PHATE2, colour = seurat_data$orig.ident)) +
  geom_point() +
  labs(colour='Cell Type')
graphics.off()

phate_out <- phate(phate_dat, t = 50)

png(paste0(plot_path, "phate_hh4_t50.png"), width=20, height=20, units = 'cm', res = 200)
ggplot(phate_out, aes(PHATE1, PHATE2, colour = seurat_data$orig.ident)) +
  geom_point() +
  labs(colour='Cell Type')
graphics.off()

phate_out <- phate(phate_dat, gamma = 0)

png(paste0(plot_path, "phate_hh4_gamma0.png"), width=20, height=20, units = 'cm', res = 200)
ggplot(phate_out, aes(PHATE1, PHATE2, colour = seurat_data$orig.ident)) +
  geom_point() +
  labs(colour='Cell Type')
graphics.off()

phate_out <- phate(phate_dat, knn = 20, t = 50)

png(paste0(plot_path, "phate_hh4_knn20_t50.png"), width=20, height=20, units = 'cm', res = 200)
ggplot(phate_out, aes(PHATE1, PHATE2, colour = seurat_data$orig.ident)) +
  geom_point() +
  labs(colour='Cell Type')
graphics.off()

phate_out <- phate(phate_dat, knn = 20, t = 50, gamma = 0)

png(paste0(plot_path, "phate_hh4_knn20_t50_gamma0.png"), width=20, height=20, units = 'cm', res = 200)
ggplot(phate_out, aes(PHATE1, PHATE2, colour = seurat_data$orig.ident)) +
  geom_point() +
  labs(colour='Cell Type')
graphics.off()