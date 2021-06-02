#!/usr/bin/env Rscript

# Load packages
library(getopt)
library(Seurat)
library(future)
# library(cowplot)
# library(clustree)
# library(gridExtra)
# library(grid)
# library(pheatmap)
# library(RColorBrewer)
# library(tidyverse)
# library(scHelper)

spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores', 'c', 2, "integer"
#   'clustres', 'r', 2, "numeric",
#   'npcs', 'p', 2, "integer",
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    plot_path = "./output/NF-downstream_analysis_stacas/seurat/7_split_stage/plots/"
    rds_path = "./output/NF-downstream_analysis_stacas/seurat/7_split_stage/rds_files/"
    data_path = "./output/NF-downstream_analysis_stacas/seurat/6_contamination_filt/rds_files/"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multiprocess", workers = ncores)
    options(future.globals.maxSize = 16* 1024^3) # 32gb
    
    # if(is.null(opt.clustres) && is.null(opt.npcs)){
    #   stop("--clustres and --npcs must be specified")
    # }
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

split_stage_data <- readRDS(list.files(data_path, full.names = TRUE))

#####################################################################################
#                         Split stage, scale and re-cluster                         #
#####################################################################################
# Set RNA to default assay
# DefaultAssay(split_stage_data) <- "RNA"

split_stage_data <- SplitObject(split_stage_data, split.by = 'stage')

# save RDS object for each stage
for(stage in names(split_stage_data)){
  saveRDS(split_stage_data[[stage]], paste0(rds_path, stage, "_split_stage_data.RDS"), compress = FALSE)
}


# # Re-run findvariablefeatures and scaling
# split_stage_data <- lapply(split_stage_data, function(x){FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, assay = 'RNA')})
# split_stage_data <- lapply(split_stage_data, function(x){ScaleData(x, features = rownames(x), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))})

# # Set Integrated to default assay
# split_stage_data <- lapply(split_stage_data, function(x){
#   assign(DefaultAssay(x),  "integrated")
# })

# # set default assay to integrated data
# for(stage in names(split_stage_data)){DefaultAssay(split_stage_data[[stage]]) <- 'integrated'}

# # Rescale data on integrated assay
# split_stage_data <- lapply(split_stage_data, function(x){ScaleData(x, features = rownames(x), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))})

# # PCA
# split_stage_data <- lapply(split_stage_data, RunPCA, verbose = FALSE)

# # Plot PC info
# for(stage in names(split_stage_data)){
#   png(paste0(plot_path, '/', stage, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
#   DimHeatmap(split_stage_data[[name]], dims = 1:30, balanced = TRUE, cells = 500)
#   graphics.off()
  
#   png(paste0(plot_path, '/', stage, "ElbowCutoff.png"), width=30, height=20, units = 'cm', res = 200)
#   ElbowCutoff(split_stage_data[[name]], return = 'plot')
#   graphics.off()
  
#   png(paste0(plot_path, '/', stage, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
#   PCALevelComparison(split_stage_data[[name]], PCA_levels = c(10, 20, 30, 40), cluster_res = 0.5)
#   graphics.off()
# }

# # Re-run 
# split_stage_data <- lapply(split_stage_data, function(x){
#   pc_cutoff <- ElbowCutoff(x)
#   FindNeighbors(x, dims = 1:pc_cutoff, verbose = FALSE)
#   RunUMAP(x, dims = 1:pc_cutoff, verbose = FALSE)
# })

# # Find optimal cluster resolution
# for(stage in names(split_stage_data)){
#   png(paste0(plot_path, '/', stage, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
#   ClustRes(seurat_object = split_stage_data[[stage]], by = 0.2, prefix = "integrated_snn_res.")
#   graphics.off()
# }

# # Use clustering resolution = 1 in order to identify any remaining poor quality cells
# split_stage_data <- lapply(split_stage_data, FindClusters, resolution = 1)

# # Plot UMAP for clusters and developmental stage
# for(stage in names(split_stage_data)){
#   png(paste0(plot_path, '/', stage, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
#   ClustStagePlot(split_stage_data[[stage]], stage_col = "stage")
#   graphics.off()
  
#   # Plot QC for each cluster
#   png(paste0(plot_path, stage, "QCPlot.png"), width=40, height=28, units = 'cm', res = 200)
#   QCPlot(split_stage_data[[stage]])
#   graphics.off()
# }


# # Find differentially expressed genes and plot heatmap of top DE genes for each cluster
# for(stage in names(split_stage_data)){
#   markers <- FindAllMarkers(split_stage_data[[stage]], only.pos = T, logfc.threshold = 0.25)
#   # get automated cluster order based on percentage of cells in adjacent stages
#   cluster_order = OrderCellClusters(seurat_object = split_stage_data[[stage]], col_to_sort = seurat_clusters, sort_by = stage)
#   # Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
#   top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) %>% arrange(factor(cluster, levels = cluster_order))
  
#   png(paste0(plot_path, stage, '/HM.top15.DE.split_stage_data.png'), height = 75, width = 100, units = 'cm', res = 500)
#   TenxPheatmap(data = split_stage_data[[stage]], metadata = c("seurat_clusters", "stage"), custom_order_column = "seurat_clusters",
#                custom_order = cluster_order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters", assay = 'integrated')
#   graphics.off()
# }

# # save RDS object for each stage
# for(stage in names(split_stage_data)){
#   saveRDS(split_stage_data[[stage]], paste0(rds_path, stage, "_split_stage_data.RDS"), compress = FALSE)
# }
