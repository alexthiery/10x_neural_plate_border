#!/usr/bin/env Rscript

# Load packages
library(getopt)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)

spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores', 'c', 2, "integer",
  'assay', 'a', 2, "character",
  'outfile', 'o', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
cat('pipeline running through Nextflow\n')
data_path = "./input/rds_files/"

seurat_object <- readRDS(list.files(data_path, full.names = TRUE))

DefaultAssay(seurat_object) <- opt$assay
seurat_object <- DietSeurat(seurat_object, counts = TRUE, assays = opt$assay, dimreducs = c('pca', 'umap'))

SaveH5Seurat(seurat_object, filename = paste0(opt$outfile, '.h5Seurat'))
Convert(paste0(opt$outfile, '.h5Seurat'), dest = "h5ad")