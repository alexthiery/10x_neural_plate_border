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
png(paste0(plot_path, "scHelper_celltype_umap2.png"), width=20, height=20, units = 'cm', res = 200)
DimPlot(seurat_data, group.by = opt$group_by, label = TRUE, label.size = 3, label.box = TRUE) + ggplot2::theme(legend.position = "none")
graphics.off()

# read in scHelper_cell_type order and colours for every script
scHelper_cell_type_order <- c('extra_embryonic', 'NNE', 'prospective_epidermis', 'PPR', 'aPPR', 'pPPR',
                              'early_NPB', 'NPB', 'aNPB', 'pNPB','NC', 'delaminating_NC',
                              'early_neural', 'early_caudal_neural', 'NP', 'pNP', 'hindbrain', 'iNP', 'midbrain', 
                              'aNP', 'forebrain', 'ventral_forebrain', 'node', 'streak')

scHelper_ann_colours <- c("#676060", "#AD2828", "#551616", "#FF0000", "#DE4D00", "#FF8300",
                          "#C8E81E", "#A5E702", "#6EE702", "#16973F", "#19B4A1", "#10E0E8",
                          "#BA3CA5", "#8A4FC5", "#0A0075", "#3B0075", "#8000FF", "#D800FF",
                          "#FF00D4", "#F16DDB", "#FFBAF3", "#B672AA", "#BBBEBE", "#787878")
names(scHelper_ann_colours) <- scHelper_cell_type_order
##################################
# set levels and extract appropriate colours depending on seurat obj
cell_type_order <- scHelper_cell_type_order[scHelper_cell_type_order %in% unique(seurat_data@meta.data[["scHelper_cell_type"]])]
seurat_data@meta.data$scHelper_cell_type <- factor(seurat_data@meta.data$scHelper_cell_type, levels = cell_type_order)
cols = scHelper_ann_colours[levels(seurat_data@meta.data$scHelper_cell_type)]
names(cols) <- NULL

# Plot pretty DimPlot
DimPlot(seurat_data, group.by = 'scHelper_cell_type', label = TRUE, label.size = 5, 
        label.box = TRUE, repel = TRUE,
        pt.size = 2, cols = cols) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none", 
                   plot.title = element_blank())

saveRDS(seurat_data, paste0(rds_path, 'seurat_label_transfer.RDS'))