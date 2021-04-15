#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
reticulate::use_python('/usr/bin/python3.7')
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


pdf("phateR.output.pdf", width = 10, height = 10)
ggplot(phate_out, aes(PHATE1, PHATE2)) +
  geom_point() +
  labs(colour='Cell Type')
dev.off()




# try running phate with just later two stages
phate_dat <- t(GetAssayData(norm.data.clustfilt.cc, slot = "counts"))

phate_dat <- phate_dat[!grepl("hh4", rownames(phate_dat)) & !grepl("hh6", rownames(phate_dat)),]
temp <- norm.data.clustfilt.cc@meta.data[!grepl("hh4", rownames(norm.data.clustfilt.cc@meta.data)) & !grepl("hh6", rownames(norm.data.clustfilt.cc@meta.data)),]


phate_dat <- library.size.normalize(phate_dat)
phate_dat <- sqrt(phate_dat)

phate_out <- phate(phate_dat)


pdf("phateR.output2.pdf", width = 10, height = 10)
ggplot(phate_out, aes(PHATE1, PHATE2, color=temp$celltype1)) +
  geom_point() +
  labs(colour='Cell Type')
dev.off()


temp <- subset(norm.data.clustfilt.cc, cells = rownames(norm.data.clustfilt.cc@meta.data[!grepl("hh4", rownames(norm.data.clustfilt.cc@meta.data)),]))
temp <- RunPCA(object = temp, verbose = FALSE)
temp <- FindNeighbors(temp, dims = 1:15, verbose = FALSE)
temp <- RunUMAP(temp, dims = 1:15, verbose = FALSE)
temp <- FindClusters(temp, resolution = 0.5)

DimPlot(temp, group.by = "celltype1")






