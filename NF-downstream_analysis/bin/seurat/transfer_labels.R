#!/usr/bin/env Rscript

# Load packages
library(optparse)
library(Seurat)
library(tidyverse)
library(scHelper)


# Read in command line opts
option_list <- list(
    make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
    make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
    make_option(c("-g", "--group_by"), action = "store", type = "character", help = "Column to group metadata by", default = 'scHelper_cell_type'),
    make_option(c("", "--verbose"), action = "store_true", type = "logical", help = "Verbose", default = FALSE)
    )

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

# Set paths and load data
plot_path = "./plots/"
rds_path = "./rds_files/"
data_path = "./input/rds_files/"
ncores = opt$cores

# Multi-core when running from command line
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 16* 1024^3) # 32gb

cat(paste0("script ran with ", ncores, " cores\n"))
dir.create(plot_path, recursive = T)
dir.create(rds_path, recursive = T)
  
# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/hh4_splitstage_data/seurat/stage_cluster/rds_files/seurat_data.RDS')

files <- list.files(data_path, recursive = TRUE, full.names = TRUE)

cell_states <- lapply(files, readRDS) %>% lapply(., function(x) x@meta.data[opt$group_by]) %>% do.call('rbind', .)


# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/seurat_filtering/6_contamination_filt/rds_files/contamination_filt_data.RDS')

seurat_data@meta.data[[opt$group_by]] <- cell_states[match(rownames(seurat_data@meta.data), rownames(cell_states)), ]
