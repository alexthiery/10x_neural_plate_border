#!/usr/bin/env Rscript

# Load packages
library(getopt)
library(Seurat)
library(SeuratDisk)

spec = matrix(c(
  'assay', 'a', 2, "character",
  'outfile', 'o', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
cat('pipeline running through Nextflow\n')
data_path = "./input/rds_files/"

seurat_object <- readRDS(list.files(data_path, full.names = TRUE, recursive = TRUE))

DefaultAssay(seurat_object) <- opt$assay
seurat_object <- DietSeurat(seurat_object, counts = TRUE, assays = opt$assay, dimreducs = c('pca', 'umap'))

# SaveH5Seurat sometimes encounters a recursion error. File is already written by this point so error can be ignored with try().
try(SaveH5Seurat(seurat_object, filename = paste0(opt$outfile, '.h5Seurat')), silent = TRUE)
Convert(paste0(opt$outfile, '.h5Seurat'), dest = "h5ad")

# Remove intermediate h5Seurat file
file.remove(paste0(opt$outfile, '.h5Seurat'))