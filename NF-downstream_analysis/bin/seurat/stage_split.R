#!/usr/bin/env Rscript

# Load packages
library(getopt)
library(Seurat)

spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores', 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
plot_path = "./plots/"
rds_path = "./rds_files/"
data_path = "./input/rds_files/"

dir.create(plot_path, recursive = T)
dir.create(rds_path, recursive = T)

split_stage_data <- readRDS(list.files(data_path, full.names = TRUE))

split_stage_data <- SplitObject(split_stage_data, split.by = 'stage')

# save RDS object for each stage
for(stage in names(split_stage_data)){
  saveRDS(split_stage_data[[stage]], paste0(rds_path, stage, "_split_stage_data.RDS"), compress = FALSE)
}