#!/usr/bin/env Rscript

# Load packages
library(Seurat)
library(SeuratDisk)
library(scHelper)
library(optparse)

# Read in command line opts
option_list <- list(
  make_option(c("-a", "--assay"), action = "store", type = "character", help = "Assay to export from seurat object ('integrated' or 'RNA')", default = 'integrated'),
  make_option(c("-o", "--outfile"), action = "store", type = "character", help = "Name of outfile"),
  make_option(c("-g", "--group_by"), action = "store", type = "character", help = "Name of metadata column containing groups to colour by"),
  make_option(c("", "--verbose"), action = "store_true", type = "logical", help = "Verbose", default = FALSE)
)

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

# Set paths and load data
cat('pipeline running through Nextflow\n')
data_path = "./input/rds_files/"

# # Options for testing
# opt <- list()
# opt$assay = 'integrated'
# opt$outfile = 'seurat'
# opt$group_by = 'seurat_clusters'
# seurat_object <- readRDS("~/output/NF-downstream_analysis_stacas/seurat/6_contamination_filt/rds_files/contamination_filt_data.RDS")

seurat_object <- readRDS(list.files(data_path, full.names = TRUE, recursive = TRUE))

DefaultAssay(seurat_object) <- opt$assay
seurat_object <- DietSeurat(seurat_object, counts = TRUE, assays = opt$assay, dimreducs = c('pca', 'umap'))

# remove anything from misc slot as depending on the contents this can cause recursion errors
seurat_object@misc <- list()

# if --group_by is specified, generate cell colours gor group_by column
if(!is.null(opt[['group_by']]) && opt[['group_by']] %in% colnames(seurat_object@meta.data)){
  # Get ggcolours for cell states
  colours = ggPlotColours(length(unique(seurat_object@meta.data[[opt$group_by]])))
  cell_state = as.factor(seurat_object@meta.data[[opt$group_by]])
  names(colours) <- levels(cell_state)
  
  # Add colours to new metadata colummn
  seurat_object@meta.data[['cell_colours']] <- unname(colours[cell_state])
  # If any na values are present in seurat identities, set colour to grey
  seurat_object@meta.data[['cell_colours']][is.na(seurat_object@meta.data[['cell_colours']])] <- '#AAAAAA'
  
}else if(!is.null(opt[['group_by']])){
  stop('--group_by missing from seurat@meta.data')
}

# SaveH5Seurat sometimes encounters a recursion error. File is already written by this point so error can be ignored with try().
# try(SaveH5Seurat(seurat_object, filename = paste0(opt$outfile, '.h5Seurat')), silent = TRUE)
SaveH5Seurat(seurat_object, filename = paste0(opt$outfile, '.h5Seurat'))
Convert(paste0(opt$outfile, '.h5Seurat'), dest = "h5ad")

# Remove intermediate h5Seurat file
file.remove(paste0(opt$outfile, '.h5Seurat'))