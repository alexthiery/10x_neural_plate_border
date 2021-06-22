#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
library(future)
library(Seurat)
library(pheatmap)
library(tidyverse)
library(Antler)
library(RColorBrewer)
library(scHelper)

spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    plot_path = "./output/NF-downstream_analysis_stacas/gene_modules/plots/"
    rds_path = "./output/NF-downstream_analysis_stacas/gene_modules/rds_files/"
    antler_path = "./output/NF-downstream_analysis_stacas/gene_modules/antler_data/"
    data_path = "./output/NF-downstream_analysis_stacas/seurat/6_contamination_filt/rds_files/"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    antler_path = "./antler_data/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multiprocess", workers = ncores)
    options(future.globals.maxSize = 16* 1024^3) # 32gb
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
  dir.create(antler_path, recursive = T)
}

seurat_data <- readRDS(list.files(data_path, full.names = TRUE))

# switch to RNA assay for viewing expression data
DefaultAssay(seurat_data) <- "RNA"

seurat_data <- DietSeurat(seurat_data, counts = TRUE, data = TRUE, scale.data = TRUE, assays = 'RNA')


# seurat_data <- subset(seurat_data, cells = colnames(seurat_data)[1:2500])


# strip end of cell names as this is incorrectly reformated in Antler
seurat_data <- RenameCells(seurat_data, new.names = gsub('-', '_', colnames(seurat_data)))

antler_data <- data.frame(row.names = colnames(seurat_data),
                          "timepoint" = as.numeric(substring(colnames(seurat_data), 3, 3)),
                          "treatment" = rep("null", ncol(seurat_data)),
                          "replicate_id" = rep(1, ncol(seurat_data))
)

# save pheno data
write.table(antler_data, file = paste0(antler_path, "phenoData.csv"), row.names = T, sep = "\t", col.names = T)

# save count data
write.table(GetAssayData(seurat_data, assay = "RNA", slot = "counts"), file = paste0(antler_path, "assayData.csv"), row.names = T, sep = "\t", col.names = T, quote = F)

########################################################################################################
#                            Load Antler data and generate correlation matrix                          #
########################################################################################################

# Change plot path
antler_data <- Antler$new(output_folder = plot_path, num_cores = ncores)
antler_data$load_dataset(folder_path = antler_path)
antler_data$exclude_unexpressed_genes(min_cells=10, min_level=1, verbose=T, data_status='Raw')

antler_data$normalize(method = 'MR')

########################################################################################################################################################

# Calculate 200 GMs
antler_data$gene_modules$identify(
  name                  = "GMs200",
  corr_t                = 0.3,  # the Spearman correlation treshold
  corr_min              = 3,    # min. number of genes a gene must correlate with
  mod_min_cell          = 10,   # min. number of cells expressing the module
  mod_consistency_thres = 0.4,  # ratio of expressed genes among "positive" cells
  num_initial_gms       = 200,
  process_plots         = TRUE)

# plot all gene modules
ncell = ncol(seurat_data)
ngene = length(unlist(antler_data$gene_modules$lists$GMs200$content))

metadata = c("stage", "seurat_clusters", "run")

png(paste0(plot_path, 'allmodules_200.png'), height = round(ngene/10), width = 75, units = 'cm', res = 600)
GeneModulePheatmap(seurat_obj = seurat_data, metadata = metadata, gene_modules = antler_data$gene_modules$lists$GMs200$content,
                   show_rownames = FALSE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 13)
graphics.off()


# use bait genes to filter mods
bait_genes = c("PAX7", "SOX2", "SOX21", "SOX10", "EYA2", "GBX2", "PAX6", "PAX2", "SIX3", "FRZB", "MSX1", "WNT1", "DLX5", "TFAP2A", "TFAP2B", "AXUD1", "GATA2", "HOMER2", "SIX1", "EYA2", "ETS1")
temp_gms <- lapply(antler_data$gene_modules$lists$GMs200$content, function(x) if(any(bait_genes %in% x)){x})
temp_gms <- temp_gms[!sapply(temp_gms, is.null)]

ncell = ncol(seurat_data)
ngene = length(unlist(temp_gms))

png(paste0(plot_path, 'bait_allmodules_200.png'), height = round(ngene/2), width = 75, units = 'cm', res = 600)
GeneModulePheatmap(seurat_obj = seurat_data, metadata = metadata, gene_modules = temp_gms,
                   show_rownames = TRUE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
graphics.off()

# Plot gene modules with at least 50% of genes DE > 0.25 logFC & FDR < 0.001
gms <- DEGeneModules(seurat_data, antler_data$gene_modules$get("GMs200"), logfc = 0.25, pval = 0.001, selected_gene_proportion = 0.5)

ngene = length(unlist(gms))

png(paste0(plot_path, 'DE_rownames_allmodules_200.png'), height = round(ngene/2), width = 75, units = 'cm', res = 600)
GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = gms,
                   show_rownames = TRUE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
graphics.off()

png(paste0(plot_path, 'DE_allmodules_200.png'), height = round(ngene/2), width = 75, units = 'cm', res = 600)
GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = gms,
                   show_rownames = FALSE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
graphics.off()

# Filter gene modules which are deferentially expressed across batches
if(length(unique(seurat_data$run)) > 1){
  gms <- gms[!names(gms) %in% names(DEGeneModules(seurat_data, antler_data$gene_modules$get("GMs200"),
                                                  logfc = 0.25,
                                                  pval = 0.001,
                                                  selected_gene_proportion = 0.5,
                                                  active_ident = 'run'))]

  ngene = length(unlist(gms))

  png(paste0(plot_path, 'DE_rownames_allmodules_200_batchfilt.png'), height = round(ngene/2), width = 75, units = 'cm', res = 600)
  GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = gms,
                    show_rownames = TRUE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
  graphics.off()

  png(paste0(plot_path, 'DE_allmodules_200_batchfilt.png'), height = round(ngene/2), width = 75, units = 'cm', res = 600)
  GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = gms,
                    show_rownames = FALSE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
  graphics.off()
}


########################################################################################################################################################
# Calculate GMs unbiasedly
antler_data$gene_modules$identify(
  name                  = "unbiasedGMs",
  corr_t                = 0.3,  # the Spearman correlation treshold
  corr_min              = 3,    # min. number of genes a gene must correlate with
  mod_min_cell          = 10,   # min. number of cells expressing the module
  mod_consistency_thres = 0.4,  # ratio of expressed genes among "positive" cells
  process_plots         = TRUE)

# plot all gene modules
ncell = ncol(seurat_data)
ngene = length(unlist(antler_data$gene_modules$lists$unbiasedGMs$content))

metadata = c("stage", "seurat_clusters", "run")

png(paste0(plot_path, 'allmodules_unbiased.png'), height = round(ngene/10), width = 75, units = 'cm', res = 600)
GeneModulePheatmap(seurat_obj = seurat_data, metadata = metadata, gene_modules = antler_data$gene_modules$lists$unbiasedGMs$content,
                   show_rownames = FALSE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15)
graphics.off()


# use bait genes to filter mods
bait_genes = c("PAX7", "SOX2", "SOX21", "SOX10", "EYA2", "GBX2", "PAX6", "PAX2", "SIX3", "FRZB", "MSX1", "WNT1", "DLX5", "TFAP2A", "TFAP2B", "AXUD1", "GATA2", "HOMER2", "SIX1", "EYA2", "ETS1")
temp_gms <- lapply(antler_data$gene_modules$lists$unbiasedGMs$content, function(x) if(any(bait_genes %in% x)){x})
temp_gms <- temp_gms[!sapply(temp_gms, is.null)]

ncell = ncol(seurat_data)
ngene = length(unlist(temp_gms))

png(paste0(plot_path, 'bait_allmodules_unbiased.png'), height = round(ngene/2), width = 75, units = 'cm', res = 600)
GeneModulePheatmap(seurat_obj = seurat_data, metadata = metadata, gene_modules = temp_gms,
                   show_rownames = TRUE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
graphics.off()




# Plot gene modules with at least 50% of genes DE > 0.25 logFC & FDR < 0.001
gms <- DEGeneModules(seurat_data, antler_data$gene_modules$get("unbiasedGMs"), logfc = 0.25, pval = 0.001, selected_gene_proportion = 0.5)

ngene = length(unlist(gms))

png(paste0(plot_path, 'DE_rownames_allmodules_unbiased.png'), height = round(ngene/2), width = 75, units = 'cm', res = 600)
GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = gms,
                   show_rownames = TRUE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
graphics.off()

png(paste0(plot_path, 'DE_allmodules_unbiased.png'), height = round(ngene/2), width = 75, units = 'cm', res = 600)
GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = gms,
                   show_rownames = FALSE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
graphics.off()

# Filter gene modules which are deferentially expressed across batches
if(length(unique(seurat_data$run)) > 1){
  gms <- gms[!names(gms) %in% names(DEGeneModules(seurat_data, antler_data$gene_modules$get("unbiasedGMs"),
                                                  logfc = 0.25,
                                                  pval = 0.001,
                                                  selected_gene_proportion = 0.5,
                                                  active_ident = 'run'))]

  ngene = length(unlist(gms))

  png(paste0(plot_path, 'DE_rownames_allmodules_unbiased_batchfilt.png'), height = round(ngene/2), width = 75, units = 'cm', res = 600)
  GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = gms,
                    show_rownames = TRUE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
  graphics.off()

  png(paste0(plot_path, 'DE_allmodules_unbiased_batchfilt.png'), height = round(ngene/2), width = 75, units = 'cm', res = 600)
  GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = gms,
                    show_rownames = FALSE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
  graphics.off()
}
