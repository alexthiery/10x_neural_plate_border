#!/usr/bin/env Rscript

# Load packages
library(getopt)
library(Seurat)
library(optparse)

# spec = matrix(c(
#   'runtype', 'l', 2, "character",
#   'cores', 'c', 2, "integer",
#   'groups', 's', 2, "character", # classification of cells to subset from dataset
#   'meta_col', 'm', 2, "character" # name of metadata column which relates to 'groups' above
# ), byrow=TRUE, ncol=4)

# opt = getopt(spec)

# Read in command line opts
option_list <- list(
    make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
    make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
    make_option(c("-m", "--meta_col"), action = "store", type = "character", help = "Name of metadata column containing groups to subset"),
    make_option(c("-g", "--groups"), action = "store", type = "character", help = "Classifications of cells (within meta_col) to subset from dataset. \
    If multiple classifications are used to subest, must be provided as a comma separated list i.e. --groups celltype1,celltype2"),
    make_option(c("", "--verbose"), action = "store_true", type = "logical", help = "Verbose", default = FALSE)
    )

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

# Split group opt by 
groups <- strsplit(opt$groups, ",")[[1]])

# Set paths and load data
plot_path = "./plots/"
rds_path = "./rds_files/"
data_path = "./input/rds_files/"

dir.create(plot_path, recursive = T)
dir.create(rds_path, recursive = T)

seurat_data <- readRDS(list.files(data_path, full.names = TRUE))
seurat_subset <- subset(seurat_data, subset = opt$meta_col %in% groups)

# save RDS object for each stage/run
saveRDS(seurat_subset, paste0(rds_path, "seurat_subset.RDS"), compress = FALSE)










