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
    
    plot_path = "./output/NF-downstream_analysis_stacas/antler/gene_modules/plots/"
    rds_path = "./output/NF-downstream_analysis_stacas/antler/gene_modules/rds_files/"
    gm_path = "./output/NF-downstream_analysis_stacas/antler/gene_modules/gene_module_lists/"
    antler_path = "./output/NF-downstream_analysis_stacas/antler/gene_modules/antler_data/"
    data_path = "./output/NF-downstream_analysis_stacas/seurat/6_contamination_filt/rds_files/"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    gm_path = "./gene_module_lists/"
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
  dir.create(gm_path, recursive = T)
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

png(paste0(plot_path, 'unbiasedGMs.png'), height = round(ngene/8), width = 75, units = 'cm', res = 600)
GeneModulePheatmap(seurat_obj = seurat_data, metadata = metadata, gene_modules = antler_data$gene_modules$lists$unbiasedGMs$content,
                   show_rownames = FALSE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15)
graphics.off()


########## DE GMs ##############
# Plot gene modules with at least 50% of genes DE > 0.25 logFC & FDR < 0.001
gms <- DEGeneModules(seurat_data, antler_data$gene_modules$get("unbiasedGMs"), logfc = 0.5, pval = 0.001, selected_gene_proportion = 0.5)

# save unbiasedGMs_DE in antler object
antler_data$gene_modules$set(name= "unbiasedGMs_DE", content = gms)

ncell = ncol(seurat_data)
ngene = length(unlist(antler_data$gene_modules$lists$unbiasedGMs_DE$content))

png(paste0(plot_path, 'unbiasedGMs_DE_rownames.png'), height = round(ngene/3), width = 75, units = 'cm', res = 200)
GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = antler_data$gene_modules$lists$unbiasedGMs_DE$content,
                   show_rownames = TRUE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
graphics.off()

png(paste0(plot_path, 'unbiasedGMs_DE.png'), height = round(ngene/8), width = 75, units = 'cm', res = 600)
GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = antler_data$gene_modules$lists$unbiasedGMs_DE$content,
                   show_rownames = FALSE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
graphics.off()

########## DE batch filter GMs ##############
# Filter gene modules which are deferentially expressed across batches
if(length(unique(seurat_data$run)) > 1){
  batch_gms <- DEGeneModules(seurat_data, antler_data$gene_modules$get("unbiasedGMs"), logfc = 0.5, pval = 0.001, selected_gene_proportion = 0.5, active_ident = 'run')
  gms <- antler_data$gene_modules$lists$unbiasedGMs_DE$content[!names(antler_data$gene_modules$lists$unbiasedGMs_DE$content) %in% names(batch_gms)]
  
  # save unbiasedGMs_DE in antler object
  antler_data$gene_modules$set(name= "unbiasedGMs_DE_batchfilt", content = gms)
  
  ncell = ncol(seurat_data)
  ngene = length(unlist(antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content))
  
  png(paste0(plot_path, 'unbiasedGMs_DE_batchfilt_rownames.png'), height = round(ngene/3), width = 75, units = 'cm', res = 600)
  GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content,
                     show_rownames = TRUE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
  graphics.off()
  
  png(paste0(plot_path, 'unbiasedGMs_DE_batchfilt.png'), height = round(ngene/8), width = 60, units = 'cm', res = 1000)
  GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content,
                     show_rownames = FALSE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
  graphics.off()
}


########## Bait GMs ##############
# use bait genes to filter mods
bait_genes = c("PAX7", "SOX2", "SOX21", "SOX10", "EYA2", "GBX2", "PAX6", "PAX2", "SIX3", "FRZB", "MSX1", "WNT1", "DLX5", "TFAP2A", "TFAP2B", "AXUD1", "GATA2", "HOMER2", "SIX1", "EYA2", "ETS1")
gms <- lapply(antler_data$gene_modules$lists$unbiasedGMs$content, function(x) if(any(bait_genes %in% x)){x})
gms <- gms[!sapply(gms, is.null)]

# save unbiasedGMs_DE in antler object
antler_data$gene_modules$set(name= "unbiasedGMs_bait", content = gms)

ncell = ncol(seurat_data)
ngene = length(unlist(antler_data$gene_modules$lists$unbiasedGMs_bait$content))

png(paste0(plot_path, 'unbiasedGMs_bait.png'), height = round(ngene/2), width = 75, units = 'cm', res = 600)
GeneModulePheatmap(seurat_obj = seurat_data, metadata = metadata, gene_modules = antler_data$gene_modules$lists$unbiasedGMs_bait$content,
                   show_rownames = TRUE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
graphics.off()


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

png(paste0(plot_path, 'GMs200.png'), height = round(ngene/8), width = 75, units = 'cm', res = 600)
GeneModulePheatmap(seurat_obj = seurat_data, metadata = metadata, gene_modules = antler_data$gene_modules$lists$GMs200$content,
                   show_rownames = FALSE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 13)
graphics.off()


########## DE GMs ##############
# Plot gene modules with at least 50% of genes DE > 0.25 logFC & FDR < 0.001
gms <- DEGeneModules(seurat_data, antler_data$gene_modules$get("GMs200"), logfc = 0.5, pval = 0.001, selected_gene_proportion = 0.5)

# save GMs200_DE in antler object
antler_data$gene_modules$set(name= "GMs200_DE", content = gms)

ncell = ncol(seurat_data)
ngene = length(unlist(antler_data$gene_modules$lists$GMs200_DE$content))

png(paste0(plot_path, 'GMs200_DE_rownames.png'), height = round(ngene/3), width = 75, units = 'cm', res = 200)
GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = antler_data$gene_modules$lists$GMs200_DE$content,
                   show_rownames = TRUE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
graphics.off()

png(paste0(plot_path, 'GMs200_DE.png'), height = round(ngene/8), width = 75, units = 'cm', res = 600)
GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = antler_data$gene_modules$lists$GMs200_DE$content,
                   show_rownames = FALSE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
graphics.off()

########## DE batch filter GMs ##############
# Filter gene modules which are deferentially expressed across batches
if(length(unique(seurat_data$run)) > 1){
  batch_gms <- DEGeneModules(seurat_data, antler_data$gene_modules$get("GMs200"), logfc = 0.5, pval = 0.001, selected_gene_proportion = 0.5, active_ident = 'run')
  gms <- antler_data$gene_modules$lists$GMs200_DE$content[!names(antler_data$gene_modules$lists$GMs200_DE$content) %in% names(batch_gms)]
  
  # save GMs200_DE in antler object
  antler_data$gene_modules$set(name= "GMs200_DE_batchfilt", content = gms)
  
  ncell = ncol(seurat_data)
  ngene = length(unlist(antler_data$gene_modules$lists$GMs200_DE_batchfilt$content))
  
  png(paste0(plot_path, 'GMs200_DE_batchfilt_rownames.png'), height = round(ngene/2), width = 75, units = 'cm', res = 600)
  GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = antler_data$gene_modules$lists$GMs200_DE_batchfilt$content,
                     show_rownames = TRUE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
  graphics.off()
  
  png(paste0(plot_path, 'GMs200_DE_batchfilt.png'), height = round(ngene/2), width = 75, units = 'cm', res = 600)
  GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = antler_data$gene_modules$lists$GMs200_DE_batchfilt$content,
                     show_rownames = FALSE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
  graphics.off()
}


########## Bait GMs ##############
# use bait genes to filter mods
bait_genes = c("PAX7", "SOX2", "SOX21", "SOX10", "EYA2", "GBX2", "PAX6", "PAX2", "SIX3", "FRZB", "MSX1", "WNT1", "DLX5", "TFAP2A", "TFAP2B", "AXUD1", "GATA2", "HOMER2", "SIX1", "EYA2", "ETS1")
gms <- lapply(antler_data$gene_modules$lists$GMs200$content, function(x) if(any(bait_genes %in% x)){x})
gms <- gms[!sapply(gms, is.null)]

# save GMs200_DE in antler object
antler_data$gene_modules$set(name= "GMs200_bait", content = gms)

ncell = ncol(seurat_data)
ngene = length(unlist(antler_data$gene_modules$lists$GMs200_bait$content))

png(paste0(plot_path, 'GMs200_bait.png'), height = round(ngene/2), width = 75, units = 'cm', res = 600)
GeneModulePheatmap(seurat_obj = seurat_data, metadata = metadata, gene_modules = antler_data$gene_modules$lists$GMs200_bait$content,
                   show_rownames = TRUE, col_order = metadata, col_ann_order = metadata, gaps_col = "stage", fontsize = 15, fontsize_row = 10)
graphics.off()


########## Write GMs ##############
export_antler_modules <- function(antler_object, publish_dir, names_list){
  for(gm_list in names_list){
    mods = antler_data$gene_modules$lists[[gm_list]]$content
    for (i in seq(length(mods))) {
      modname = base::names(mods)[i]
      if (is.null(modname)) {
        modname = paste0("GM: ", i)
      }
      write(paste0(modname, "; ", paste0(mods[[i]], collapse = ", ")), file = paste0(publish_dir, '/', gm_list, '.txt'), append = TRUE)
    }
  }
}

export_antler_modules(antler_data, publish_dir = gm_path, names_list = c('unbiasedGMs', 'unbiasedGMs_DE', 'unbiasedGMs_DE_batchfilt', 'unbiasedGMs_bait',
                                                                         'GMs200', 'GMs200_DE', 'GMs200_DE_batchfilt', 'GMs200_bait'))

########## Save Antler object ##############
saveRDS(antler_data, paste0(rds_path, 'antler_out.RDS'))
