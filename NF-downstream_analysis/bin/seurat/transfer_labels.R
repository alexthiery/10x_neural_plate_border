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
    make_option(c("-d", "--full_data"), action = "store", type = "character", help = "Name of full dataset to transfer labels to", default = 'contamination_filt_data.RDS'),
    make_option(c("", "--verbose"), action = "store_true", type = "logical", help = "Verbose", default = FALSE)
    )

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

# Set paths and load data
plot_path = "./plots/"
rds_path = "./rds_files/"
data_path = "./input"
ncores = opt$cores

cat(paste0("script ran with ", ncores, " cores\n"))
dir.create(plot_path, recursive = T)
dir.create(rds_path, recursive = T)

files <- list.files(data_path, full.names = TRUE, pattern = '*.RDS')
stage_data <- grep(opt$full_data, files, invert = T, value = TRUE)
full_data <- grep(opt$full_data, files, invert = F, value = TRUE)

cell_states <- lapply(stage_data, readRDS) %>% lapply(., function(x) x@meta.data[opt$group_by]) %>% do.call('rbind', .)

seurat_data <- readRDS(full_data)

seurat_data@meta.data[[opt$group_by]] <- cell_states[match(rownames(seurat_data@meta.data), rownames(cell_states)), ]

# Plot QC for each cluster
png(paste0(plot_path, "UMAP.png"), width=30, height=25, units = 'cm', res = 200)
DimPlot(seurat_data, group.by = opt$group_by)
graphics.off()

saveRDS(seurat_data, paste0(rds_path, 'seurat_label_transfer.RDS'))