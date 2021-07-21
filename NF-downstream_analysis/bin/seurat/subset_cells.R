#!/usr/bin/env Rscript

# Load packages
library(getopt)
library(Seurat)

spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores', 'c', 2, "integer",
  'groups', 's', 2, "character",
  'metadata_col', 'm', 2, "character"
), byrow=TRUE, ncol=4)

opt = getopt(spec)

# Set paths and load data
plot_path = "./plots/"
rds_path = "./rds_files/"
data_path = "./input/rds_files/"

dir.create(plot_path, recursive = T)
dir.create(rds_path, recursive = T)

seurat_data <- readRDS(list.files(data_path, full.names = TRUE))

seurat_subset <- subset(seurat_data, subset = scHelper_cell_type %in% celltypes)

# save RDS object for each stage/run
saveRDS(seurat_subset, paste0(rds_path, "seurat_subset.RDS"), compress = FALSE)
