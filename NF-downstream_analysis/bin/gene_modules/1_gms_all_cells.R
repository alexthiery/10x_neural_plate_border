#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
reticulate::use_python('/usr/bin/python3.7')
library(Seurat)

library(pheatmap)
library(tidyverse)
library(Antler)

spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer",
  'custom_functions', 'm', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set run location
if(length(commandArgs(trailingOnly = TRUE)) == 0){
  cat('No command line arguments provided, user defaults paths are set for running interactively in Rstudio on docker\n')
  opt$runtype = "user"
} else {
  if(is.null(opt$runtype)){
    stop("--runtype must be either 'user' or 'nextflow'")
  }
  if(tolower(opt$runtype) != "user" & tolower(opt$runtype) != "nextflow"){
    stop("--runtype must be either 'user' or 'nextflow'")
  }
  if(tolower(opt$runtype) == "nextflow"){
    if(is.null(opt$custom_functions) | opt$custom_functions == "null"){
      stop("--custom_functions path must be specified in process params config")
    }
  }
}

# Set paths and load data
{
  if (opt$runtype == "user"){
    sapply(list.files('./NF-downstream_analysis/bin/custom_functions/', full.names = T), source)
    plot_path = "./output/NF-downstream_analysis/1_gms_all_cells/plots/"
    rds_path = "./output/NF-downstream_analysis/1_gms_all_cells/rds_files/"
    antler_path = "./output/NF-downstream_analysis/1_gms_all_cells/antler_data/"
    data_path = "./output/NF-downstream_analysis/6_contamination_filt/rds_files/"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    sapply(list.files(opt$custom_functions, full.names = T), source)
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    antler_path = "./antler_data/"
    data_path = "./input/rds_files/"
    ncores = opt$cores

    # Multi-core when running from command line
    plan("multiprocess", workers = ncores)
    options(future.globals.maxSize = 32* 1024^3) # 32gb
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
  dir.create(antler_path, recursive = T)
}

seurat_data <- readRDS(paste0(data_path, 'contamination_filt_data.RDS'))

# strip end of cell names as this is incorrectly reformated in Antler
seurat_data <- RenameCells(seurat_data, new.names = sub('-', '_', colnames(seurat_data)))

seurat_data <- RenameCells(seurat_data, new.names = sub('-.*', '', colnames(seurat_data)))

antler_data <- data.frame(row.names = colnames(seurat_data),
                          "timepoint" = as.numeric(substring(colnames(seurat_out), 3, 3)),
                          "treatment" = rep("null", ncol(seurat_out)),
                          "replicate_id" = rep(1, ncol(seurat_out))
)

# save pheno data
write.table(antler_data, file = paste0(antler.path, "phenoData.csv"), row.names = T, sep = "\t", col.names = T)

# save count data
write.table(GetAssayData(seurat_out, assay = "RNA", slot = "counts"), file = paste0(antler_path, "assayData.csv"), row.names = T, sep = "\t", col.names = T, quote = F)

########################################################################################################
#                            Load Antler data and generate correlation matrix                          #
########################################################################################################

# Change plot path
antler_data <- Antler$new(output_folder = plot_path, num_cores = ncores)
antler_data$load_dataset(folder_path = antler_path)
antler_data$exclude_unexpressed_genes(min_cells=10, min_level=1, verbose=T, data_status='Raw')

antler_data$normalize(method = 'MR')

antler_data$gene_modules$identify(
  name                  = "unbiasedGMs",
  corr_t                = 0.3,  # the Spearman correlation treshold
  corr_min              = 3,    # min. number of genes a gene must correlate with
  mod_min_cell          = 10,   # min. number of cells expressing the module
  mod_consistency_thres = 0.4,  # ratio of expressed genes among "positive" cells
  num_initial_gms       = 200,
  process_plots         = TRUE)

saveRDS(antler_data, paste0(rds_path, "antler_all.RDS"))
# antler <- readRDS(paste0(rds.path, "antler_all.RDS"))

# plot all gene modules
png(paste0(plot_path, 'allmodules.200.png'), height = 100, width = 80, units = 'cm', res = 400)
GM.plot(data = seurat_data, metadata = c("stage", "seurat_clusters"), gene_modules = antler_data$gene_modules$lists$unbiasedGMs$content,
        show_rownames = F, col_order = c("stage", "seurat_clusters"))
graphics.off()

# use bait genes to filter mods
bait.genes = c("PAX7", "SOX2", "SOX21", "SOX10", "EYA2", "GBX2", "PAX6", "PAX2", "SIX3", "FRZB", "MSX1", "WNT1", "DLX5", "TFAP2A", "TFAP2B", "AXUD1", "GATA2", "HOMER2", "SIX1", "EYA2", "ETS1")
temp.gms = lapply(antler$gene_modules$lists$unbiasedGMs$content, function(x) if(any(bait.genes %in% x)){x})

png(paste0(plot.path, 'DE.GM.200.png'), height = 50, width = 80, units = 'cm', res = 400)
GM.plot(data = seurat_out, metadata = c("stage", "orig.ident", "seurat_clusters"), gene_modules = temp.gms, gaps_col = "stage",
        show_rownames = T, col_order = c("stage", "seurat_clusters"))
graphics.off()








antler$gene_modules$identify(
  name                  = "unbiasedGMs",
  corr_t                = 0.3,  # the Spearman correlation treshold
  corr_min              = 3,    # min. number of genes a gene must correlate with
  mod_min_cell          = 10,   # min. number of cells expressing the module
  mod_consistency_thres = 0.4,  # ratio of expressed genes among "positive" cells
  process_plots         = TRUE)

saveRDS(antler, paste0(rds.path, "antler_all.RDS"))
# antler <- readRDS(paste0(rds.path, "antler_all.RDS"))

# plot all gene modules
png(paste0(plot.path, 'allmodules.unbiased.png'), height = 100, width = 80, units = 'cm', res = 400)
GM.plot(data = seurat_out, metadata = c("stage", "seurat_clusters"), gene_modules = antler$gene_modules$lists$unbiasedGMs$content,
        show_rownames = F, col_order = c("stage", "seurat_clusters"))
graphics.off()

# use bait genes to filter mods
bait.genes = c("PAX7", "SOX2", "SOX21", "SOX10", "EYA2", "GBX2", "PAX6", "PAX2", "SIX3", "FRZB", "MSX1", "WNT1", "DLX5", "TFAP2A", "TFAP2B", "AXUD1", "GATA2", "HOMER2", "SIX1", "EYA2", "ETS1")
temp.gms = lapply(antler$gene_modules$lists$unbiasedGMs$content, function(x) if(any(bait.genes %in% x)){x})

png(paste0(plot.path, 'DE.GM.unbiased.png'), height = 50, width = 80, units = 'cm', res = 400)
GM.plot(data = seurat_out, metadata = c("stage", "orig.ident", "seurat_clusters"), gene_modules = temp.gms, gaps_col = "stage",
        show_rownames = T, col_order = c("stage", "seurat_clusters"))
graphics.off()