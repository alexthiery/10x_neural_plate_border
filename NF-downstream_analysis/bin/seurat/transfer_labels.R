#!/usr/bin/env Rscript

# Load packages
library(optparse)
library(Seurat)
library(tidyverse)
library(scHelper)

########################       CELL STATE COLOURS    ########################################
scHelper_all_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR',
                        'eNPB', 'NPB', 'aNPB', 'pNPB','NC', 'dNC',
                        'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                        'aNP', 'FB', 'vFB', 'node', 'streak')

scHelper_all_colours <- c("#676060", "#AD2828", "#551616", "#FF0000", "#DE4D00", "#FF8300",
                          "#C8E81E", "#A5E702", "#6EE702", "#16973F", "#19B4A1", "#10E0E8",
                          "#BA3CA5", "#8A4FC5", "#0A0075", "#3B0075", "#8000FF", "#D800FF",
                          "#FF00D4", "#F16DDB", "#FFBAF3", "#B672AA", "#BBBEBE", "#787878")
names(scHelper_all_colours) <- scHelper_all_order
########################       STAGE COLOURS     ###########################################
stage_all_order <- c("hh4", "hh5", "hh6", "hh7", "ss4", "ss8")

stage_all_colours = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F")
names(stage_all_colours) <- stage_all_order
############################################################################################

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


#####################   Set levels and colours for scHelper_cell_type  ###########################################
scHelper_order <- scHelper_all_order[scHelper_all_order %in% unique(seurat_data@meta.data[["scHelper_cell_type"]])]
seurat_data@meta.data$scHelper_cell_type <- factor(seurat_data@meta.data$scHelper_cell_type, levels = scHelper_all_order)

scHelper_cols = scHelper_all_colours[levels(seurat_data@meta.data$scHelper_cell_type)]
names(scHelper_cols) <- NULL

#####################   Set levels and colours for stage   ###########################################
stage_order <- stage_all_order[stage_all_order %in% unique(seurat_data@meta.data[["stage"]])]
seurat_data@meta.data$stage <- factor(seurat_data@meta.data$stage, levels = stage_all_order)

stage_cols = stage_all_colours[levels(seurat_data@meta.data$stage)]
names(stage_cols) <- NULL

################ DimPlot of scHelper_cell_types
png(paste0(plot_path, "scHelper_celltype_umap.png"), width=20, height=20, units = 'cm', res = 200)
DimPlot(seurat_data, group.by = 'scHelper_cell_type', label = TRUE, label.size = 3, 
        label.box = TRUE, repel = TRUE,
        pt.size = 1, cols = scHelper_cols) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none", 
                   plot.title = element_blank())
graphics.off()

################ DimPlot of stages
png(paste0(plot_path, "stage_umap.png"), width=20, height=20, units = 'cm', res = 200)
DimPlot(seurat_data, group.by = 'stage', label = TRUE, label.size = 5, 
        label.box = TRUE, repel = TRUE,
        pt.size = 1, cols = stage_cols) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none", 
                   plot.title = element_blank())
graphics.off()

saveRDS(seurat_data, paste0(rds_path, 'seurat_label_transfer.RDS'))