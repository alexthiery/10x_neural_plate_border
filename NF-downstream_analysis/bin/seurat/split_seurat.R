#!/usr/bin/env Rscript

# Load packages
library(getopt)
library(Seurat)

spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores', 'c', 2, "integer",
  'split', 's', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

if(!opt$split %in% c('stage', 'run')){
    stop("'--split' must be one of either 'run' or 'stage'")
}

# Set paths and load data
plot_path = "./plots/"
rds_path = "./rds_files/"
data_path = "./input/rds_files/"

dir.create(plot_path, recursive = T)
dir.create(rds_path, recursive = T)

split_data <- readRDS(list.files(data_path, full.names = TRUE))

split_data <- SplitObject(split_data, split.by = opt$split)

# save RDS object for each stage/run
for(split in names(split_data)){
  saveRDS(split_data[[split]], paste0(rds_path, split, "_split", opt$split, "_data.RDS"), compress = FALSE)
}