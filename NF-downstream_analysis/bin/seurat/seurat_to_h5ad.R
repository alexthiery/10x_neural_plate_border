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

input_rds <- readRDS(list.files(data_path, full.names = TRUE))

DefaultAssay(input_rds) <- opt$assay
SaveH5Seurat(input_rds, filename = opt$outfile)
Convert(opt$outfile, dest = "h5ad")