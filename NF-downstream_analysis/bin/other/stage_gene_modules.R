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
  'cores', 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
data_path = "./input/rds_files/"
stage_cluster_data <- readRDS(list.files(data_path, full.names = TRUE))

# Set stage var based on input
if(length(unique(stage_cluster_data$stage)) == 1){
  stage = unique(stage_cluster_data$stage)
} else {
  stop('Input RDS object contains more than 1 developmental stage')
}

plot_path = paste0("./plots/", stage, "/")
rds_path = "./rds_files/"
antler_path = paste0("./antler_data/", stage, "/")
dir.create(plot_path, recursive = T)
dir.create(rds_path, recursive = T)
dir.create(antler_path, recursive = T)

# Multi-core when running from command line
ncores = opt$cores
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 4* 1024^3)
cat(paste0("script ran with ", ncores, " cores\n"))

seurat_data <- readRDS(list.files(data_path, full.names = TRUE))

# switch to RNA assay for viewing expression data
DefaultAssay(seurat_data) <- "RNA"

seurat_data <- DietSeurat(seurat_data, counts = TRUE, data = TRUE, scale.data = TRUE, assays = 'RNA')

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

# antler_data$gene_modules$identify(
#   name                  = "unbiasedGMs",
#   corr_t                = 0.3,  # the Spearman correlation treshold
#   corr_min              = 3,    # min. number of genes a gene must correlate with
#   mod_min_cell          = 10,   # min. number of cells expressing the module
#   mod_consistency_thres = 0.4,  # ratio of expressed genes among "positive" cells
#   num_initial_gms       = 100,
#   process_plots         = TRUE)

saveRDS(antler_data, paste0(rds_path, "antler_all.RDS"))
# antler <- readRDS(paste0(rds.path, "antler_all.RDS"))

# # plot all gene modules
# png(paste0(plot_path, 'allmodules_100.png'), height = 100, width = 80, units = 'cm', res = 400)
# GeneModulePheatmap(seurat_obj = seurat_data, metadata = c("stage", "seurat_clusters"), gene_modules = antler_data$gene_modules$lists$unbiasedGMs$content,
#         show_rownames = F, col_order = c("stage", "seurat_clusters"))
# graphics.off()









antler_data$gene_modules$identify(
  name                  = "unbiasedGMs",
  corr_t                = 0.3,  # the Spearman correlation treshold
  corr_min              = 3,    # min. number of genes a gene must correlate with
  mod_min_cell          = 10,   # min. number of cells expressing the module
  mod_consistency_thres = 0.4,  # ratio of expressed genes among "positive" cells
  process_plots         = TRUE)

saveRDS(antler_data, paste0(rds.path, stage, "_antler.RDS"))
# antler <- readRDS(paste0(rds.path, "antler_all.RDS"))

# plot all gene modules
png(paste0(plot.path, 'allmodules_unbiased.png'), height = 100, width = 80, units = 'cm', res = 400)
GeneModulePheatmap(seurat_obj = seurat_out, metadata = c("stage", "seurat_clusters"), gene_modules = antler_data$gene_modules$lists$unbiasedGMs$content,
        show_rownames = T, col_order = c("stage", "seurat_clusters"))
graphics.off()


# Plot gene modules with at least 50% of genes DE > 0.25 logFC & FDR < 0.001

# Find DEGs
DE_genes <- FindAllMarkers(seurat_data, only.pos = T, logfc.threshold = 0.25) %>% filter(p_val_adj < 0.001)

# Get automated cluster order based on percentage of cells in adjacent stages
cluster_order <- OrderCellClusters(seurat_object = seurat_data, col_to_sort = seurat_clusters, sort_by = stage)

# Filter GMs with 50% genes DE logFC > 0.25 & FDR < 0.001
gms <- SubsetGeneModules(antler_data$gene_modules$get("unbiasedGMs"), selected_genes = DE_genes$gene, keep_mod_ID = T, selected_gene_ratio = 0.5)

png(paste0(plot_path, 'DE_gms.png'), height = 160, width = 80, units = 'cm', res = 500)
GeneModulePheatmap(seurat_obj = seurat_data, metadata = c("stage", "seurat_clusters"), gene_modules = gms, gaps_col = "seurat_clusters", custom_order = cluster_order, custom_order_column = "seurat_clusters", fontsize = 25, fontsize_row = 10)
graphics.off()