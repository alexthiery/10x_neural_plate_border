#!/usr/bin/env Rscript

# Load packages
library(optparse)
library(Seurat)
library(tidyverse)
library(scHelper)

########################       CELL STATE COLOURS    ########################################
scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR',
                              'eNPB', 'NPB', 'aNPB', 'pNPB','NC', 'dNC',
                              'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                              'aNP', 'FB', 'vFB', 'node', 'streak')

scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", "#53A651", "#6D8470",
                                "#87638F", "#A5548D", "#C96555", "#ED761C", "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30",
                                "#CC9F2C", "#AD6428", "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3")

names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 'MB', 
                                       'vFB', 'aNP', 'node', 'FB', 'pEpi')
########################       STAGE COLOURS     ###########################################
stage_order <- c("HH4", "HH5", "HH6", "HH7", "ss4", "ss8")

stage_colours = c("#E78AC3", "#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order
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

#####################   Set levels
seurat_data@meta.data$scHelper_cell_type <- factor(seurat_data@meta.data$scHelper_cell_type, levels = scHelper_cell_type_order)
seurat_data@meta.data$stage <- factor(seurat_data@meta.data$stage, levels = stage_order)

#####################   Set colours
scHelper_cols <- scHelper_cell_type_colours[levels(droplevels(seurat_data@meta.data$scHelper_cell_type))]
stage_cols <- stage_colours[levels(droplevels(seurat_data@meta.data$stage))]

################ DimPlot of scHelper_cell_types
png(paste0(plot_path, "scHelper_celltype_umap.png"), width=26, height=20, units = 'cm', res = 200)
DimPlot(seurat_data, group.by = 'scHelper_cell_type', label = FALSE, 
        pt.size = 0.9, cols = scHelper_cols, shuffle = TRUE) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = element_blank(),
                   legend.key.size = unit(1.5, 'cm'), #change legend key size
                   legend.key.height = unit(1, 'cm'), #change legend key height
                   legend.key.width = unit(1, 'cm'), #change legend key width
                   legend.text = element_text(size=20))
graphics.off()

################ DimPlot of stages
png(paste0(plot_path, "stage_umap.png"), width=20, height=20, units = 'cm', res = 200)
DimPlot(seurat_data, group.by = 'stage', label = TRUE, label.size = 8, 
        label.box = TRUE, repel = TRUE,
        pt.size = 0.9, cols = stage_cols, shuffle = TRUE) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none", 
                   plot.title = element_blank())
graphics.off()

saveRDS(seurat_data, paste0(rds_path, 'seurat_label_transfer.RDS'))