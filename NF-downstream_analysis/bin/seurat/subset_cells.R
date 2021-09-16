#!/usr/bin/env Rscript

# Load packages
library(Seurat)
library(optparse)
library(ggplot2)
library(scHelper)
library(grid)
library(gridExtra)

# Read in command line opts
option_list <- list(
    make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
    make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
    make_option(c("-m", "--meta_col"), action = "store", type = "character", help = "Name of metadata column containing groups to subset", default = NULL),
    make_option(c("-o", "--output"), action = "store", type = "character", help = "Name of output RDS file", default = 'seurat_subset'),
    make_option(c("-g", "--groups"), action = "store", type = "character", help = "Classifications of cells (within meta_col) to subset from dataset. \
    If multiple classifications are used to subest, must be provided as a comma separated list i.e. --groups celltype1,celltype2", default = NULL),
    make_option(c("-i", "--invert"), action = "store", type = "logical", help = "Boolean for whether to invert group selection", default = FALSE),
    make_option(c("", "--verbose"), action = "store_true", type = "logical", help = "Verbose", default = FALSE))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

# opt$groups = 'early_non_neural,early_aPPR,early_pPPR,aNPB,pPPR,non_neural'
# opt$meta_col = 'scHelper_cell_type'

opt$groups = strsplit(opt$groups, ',')[[1]]

if(is.na(opt$meta_col)){
    stop("meta_col parameter must be provided. See script usage (--help)")
}

if(is.null(opt$groups)){
    stop("groups parameter must be provided. See script usage (--help)")
}

# Set paths and load data
plot_path = "./plots/"
rds_path = "./rds_files/"
data_path = "./input/rds_files/"

dir.create(plot_path, recursive = T)
dir.create(rds_path, recursive = T)

seurat_data <- readRDS(list.files(data_path, full.names = TRUE))

# If invert is true, then subset the inverted groups from the seurat object
if(opt$invert == TRUE){
    opt$groups <- as.character(unique(seurat_data@meta.data[[opt$meta_col]])[!unique(seurat_data@meta.data[[opt$meta_col]]) %in% opt$groups])
}

# Plot DimPlot of subset
seurat_data@meta.data[[opt$meta_col]] <- as.factor(seurat_data@meta.data[[opt$meta_col]])
colours = ggPlotColours(length(levels(seurat_data@meta.data[[opt$meta_col]])))

max_x = round(max(seurat_data@reductions$umap@cell.embeddings[,1]))
min_x = round(min(seurat_data@reductions$umap@cell.embeddings[,1]))
max_y = round(max(seurat_data@reductions$umap@cell.embeddings[,2]))
min_y = round(min(seurat_data@reductions$umap@cell.embeddings[,2]))

full_plot = DimPlot(seurat_data, group.by = opt$meta_col, pt.size = 0.3, cols = colours, label = TRUE, label.size = 3, label.box = TRUE) +
                ylim(min_y, max_y) +
                xlim(min_x, max_x) +
                theme(legend.position = "none") +
                ggtitle("")

subset_plot = DimPlot(seurat_data, group.by = opt$meta_col, pt.size = 0.3, cols = colours[levels(seurat_data@meta.data[[opt$meta_col]]) %in% opt$groups],
        cells = rownames(seurat_data@meta.data)[seurat_data@meta.data[[opt$meta_col]] %in% opt$groups]) +
                ylim(min_y, max_y) +
                xlim(min_x, max_x) +
                theme(legend.position = "none") +
                ggtitle("")

png(paste0(plot_path, 'cell_subset_umap.png'), width = 40, height = 20, units = 'cm', res = 200)
grid.arrange(full_plot, subset_plot, nrow = 1)
graphics.off()

# Subset cells
seurat_subset <- subset(seurat_data, subset = !!as.symbol(opt$meta_col) %in% opt$groups)

# save RDS object for each stage/run
saveRDS(seurat_subset, paste0(rds_path, opt$output, ".RDS"), compress = FALSE)