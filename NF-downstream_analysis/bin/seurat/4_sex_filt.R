#!/usr/bin/env Rscript

# Load packages
library(getopt)
library(Seurat)
library(future)
library(cowplot)
library(clustree)
library(gridExtra)
library(grid)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
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

    plot_path = "./output/NF-downstream_analysis_stacas/seurat_filtering/4_sex_filt/plots/"
    rds_path = "./output/NF-downstream_analysis_stacas/seurat_filtering/4_sex_filt/rds_files/"
    data_path = "./output/NF-downstream_analysis_stacas/seurat_filtering/3_integration_qc/rds_files/"
    
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

  } else {
    stop("--runtype must be set to 'nextflow'")
  }

  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

pre_sex_filt_data <- readRDS(list.files(data_path, full.names = TRUE))

png(paste0(plot_path, 'ClusterDimplot_pre-sex_filt.png'), height = 40, width = 70, units = 'cm', res = 500)
ClusterDimplot(pre_sex_filt_data, clusters = c(1, 6), xlim = c(-10, 15), ylim = c(-10, 15))
graphics.off()

DefaultAssay(pre_sex_filt_data) <- "RNA"

# Log normalize data and find variable features
pre_sex_filt_data <- NormalizeData(pre_sex_filt_data, normalization.method = "LogNormalize", scale.factor = 10000)
pre_sex_filt_data <- FindVariableFeatures(pre_sex_filt_data, selection.method = "vst", nfeatures = 2000)

pre_sex_filt_data <- ScaleData(pre_sex_filt_data, features = rownames(pre_sex_filt_data), vars.to.regress = "percent.mt")


# There is a strong sex effect - this plot shows DE genes between top HH4 clusters. Clustering is driven by sex genes

# Select two largest HH4 clusters to look at sex effect
sex_clusters <- pre_sex_filt_data@meta.data %>% filter(grepl("HH4", orig.ident)) %>% count(seurat_clusters) %>%
  top_n(n, n = 2) %>% pull(seurat_clusters) %>% as.character()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindMarkers(pre_sex_filt_data, ident.1 = sex_clusters[1], ident.2 = sex_clusters[2])

top30 <- markers %>% rownames_to_column("gene") %>% mutate(pos = avg_log2FC>0) %>% group_by(pos) %>% top_n(n = 30, wt = avg_log2FC)

png(paste0(plot_path, 'HM.top30.DE.pre-sex_filt.png'), height = 40, width = 70, units = 'cm', res = 500)
TenxPheatmap(data = pre_sex_filt_data[,rownames(filter(pre_sex_filt_data@meta.data, seurat_clusters %in% sex_clusters))],
              metadata = c("seurat_clusters", "orig.ident"), selected_genes = unique(top30$gene), hclust_rows = T, gaps_col = "seurat_clusters", assay = "RNA")
graphics.off()

# plot dimplot for main W gene
png(paste0(plot_path, 'Sex_genes_multifeature_pre-sex_filt.png'), height = 40, width = 80, units = 'cm', res = 200)
MultiFeaturePlot(pre_sex_filt_data, plot_stage = TRUE, gene_list = c("W-Wpkci-7", "W-HNRNPKL", "W-NIPBLL", "Z-MRPS30", "Z-FXN", "Z-RCL1"), n_col = 4)
graphics.off()

#####################################################################################################
#     Heatmap clearly shows clusters segregate by sex - check this and regress out the sex effect   #
#####################################################################################################

# Use W chromosome genes to K-means cluster the cells into male (zz) and female (zw)
W_genes <- as.matrix(pre_sex_filt_data@assays$RNA[grepl("W-", rownames(pre_sex_filt_data@assays$RNA)),])
k_clusters <- kmeans(t(W_genes), 2)
k_clusters <- data.frame(k_clusters$cluster)
pre_sex_filt_data@meta.data$k_clusters <- k_clusters[match(colnames(pre_sex_filt_data@assays$RNA), rownames(k_clusters)),]

# Get rownames for kmeans clusters 1 and 2
k_clus_1 <- rownames(pre_sex_filt_data@meta.data[pre_sex_filt_data@meta.data$k_clusters == 1,])
k_clus_2 <- rownames(pre_sex_filt_data@meta.data[pre_sex_filt_data@meta.data$k_clusters == 2,])

# K clustering identities are stochastic, so I mneed to identify which cluster is male and female
# Sum of W genes is order of magnitude greater in cluster 2 - these are the female cells
sumclus1 <- sum(W_genes[,k_clus_1])
sumclus2 <- sum(W_genes[,k_clus_2])

if(sumclus1 < sumclus2){
  k_male <- k_clus_1
  k_female <- k_clus_2
} else {
  k_female <- k_clus_1
  k_male <- k_clus_2
}

# Add sex data to meta.data
pre_sex_filt_data@meta.data$sex <- unlist(lapply(rownames(pre_sex_filt_data@meta.data), function(x)
  if(x %in% k_male){"male"} else if(x %in% k_female){"female"} else{stop("cell sex is not assigned")}))


###### Next - subset autosomal genes and Z genes - then calculate the average for each gene for both kmeans clustered cells and plot in order to compare whether the two groups show significant differences in their expression.

# This is a test on autosomal genes to try and calculate and compare FC betweeen clusters
# Calculating median is tricky as there are a lot of dropouts in 10x data so you end up with either 0s (when the median  = 0) or 1 (when the median expression in both clusters is the same - probably a result of normalisation resulting in a UMI of 0 or 1 being normalised to a nominal value)

# Make dataframe for mean Z expression in male cells
mean_Z_male <- data.frame(Z.mean = apply(pre_sex_filt_data@assays$RNA[grepl("Z-", rownames(pre_sex_filt_data@assays$RNA)), k_male], 1, mean))
# add 1 before log2 as log2(1) = 0
mean_Z_male <- log2(mean_Z_male + 1)

# Make dataframe for mean Z expression in female cells
mean_Z_female <- data.frame(Z.mean = apply(pre_sex_filt_data@assays$RNA[grepl("Z-", rownames(pre_sex_filt_data@assays$RNA)), k_female], 1, mean))
mean_Z_female <- log2(mean_Z_female + 1)

# Make dataframe for mean autosomal expression in male cells
mean_auto_male <- data.frame(auto.mean = apply(pre_sex_filt_data@assays$RNA[!grepl("Z-", rownames(pre_sex_filt_data@assays$RNA)) & !grepl("W-", rownames(pre_sex_filt_data@assays$RNA)), k_male], 1, mean))
mean_auto_male <- log2(mean_auto_male + 1)

# Make dataframe for mean autosomal expression in male cells
mean_auto_female <- data.frame(auto.mean = apply(pre_sex_filt_data@assays$RNA[!grepl("Z-", rownames(pre_sex_filt_data@assays$RNA)) & !grepl("W-", rownames(pre_sex_filt_data@assays$RNA)), k_female], 1, mean))
mean_auto_female <- log2(mean_auto_female + 1)

# Calculate FC by subtracting log2 expression from each other
FC <- list()
FC$Z <- mean_Z_male - mean_Z_female
FC$auto <-  mean_auto_male - mean_auto_female

# Plot boxplot of Z gene and autosomal expression in male vs female cells
png(paste0(plot_path,"sex_kmeans_log2FC_boxplot.png"), height = 18, width = 18, units = "cm", res = 200)
boxplot(c(FC$Z, FC$auto),  ylab = "male - female log2 FC (mean normalised UMI +1)", names = c("Z chromosome genes", "autosomal genes"))
graphics.off()

# Z genes are upregulated within male cells relative to female cells whereas autosomal genes have a normal distribution of logFCs

#####################################################################################################
#                                     Regress sex effect                                            #
#####################################################################################################

# Init sexscale object 
sex_filt_data <- pre_sex_filt_data

DefaultAssay(sex_filt_data) <- "integrated"

sex_filt_data <- ScaleData(sex_filt_data, features = rownames(sex_filt_data), vars.to.regress = c("percent.mt", "sex"), verbose = FALSE)

# PCA
sex_filt_data <- RunPCA(object = sex_filt_data, verbose = FALSE)

png(paste0(plot_path, "dimHM.png"), width=30, height=65, units = 'cm', res = 200)
DimHeatmap(sex_filt_data, dims = 1:40, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(plot_path, "ElbowCutoff.png"), width=30, height=20, units = 'cm', res = 200)
ElbowCutoff(sex_filt_data, return = 'plot')
graphics.off()

pc_cutoff <- ElbowCutoff(sex_filt_data)

png(paste0(plot_path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCALevelComparison(sex_filt_data, PCA_levels = c(pc_cutoff-5, pc_cutoff, pc_cutoff+5, pc_cutoff+10), cluster_res = 0.5)
graphics.off()

sex_filt_data <- FindNeighbors(sex_filt_data, dims = 1:pc_cutoff, verbose = FALSE)
sex_filt_data <- RunUMAP(sex_filt_data, dims = 1:pc_cutoff, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
ClustRes(seurat_object = sex_filt_data, by = 0.1, prefix = "integrated_snn_res.")
graphics.off()

# Use default clustering resolution (0.5)
sex_filt_data <- FindClusters(sex_filt_data, resolution = 0.5, verbose = FALSE)

# Plot UMAP for clusters and developmental stage
png(paste0(plot_path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
ClustStagePlot(sex_filt_data)
graphics.off()

# Plot QC for each cluster
png(paste0(plot_path, "QCPlot.png"), width=28, height=28, units = 'cm', res = 200)
QCPlot(sex_filt_data, plot_quantiles = TRUE)
graphics.off()

# switch to RNA assay for viewing expression data
DefaultAssay(sex_filt_data) <- "RNA"

# Find variable features and scale data on RNA assay
sex_filt_data <- FindVariableFeatures(sex_filt_data, selection.method = "vst", nfeatures = 2000)
sex_filt_data <- ScaleData(sex_filt_data, features = rownames(sex_filt_data), vars.to.regress = c("percent.mt", "sex"), verbose = FALSE)

# Save RDS
saveRDS(sex_filt_data, paste0(rds_path, "sex_filt_data.RDS"), compress = FALSE)

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(sex_filt_data, only.pos = T, logfc.threshold = 0.25, assay = "RNA")
# get automated cluster order based on percentage of cells in adjacent stages
cluster_order <- OrderCellClusters(seurat_object = sex_filt_data, col_to_sort = 'seurat_clusters', sort_by = 'orig.ident')
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) %>% arrange(factor(cluster, levels = cluster_order))

png(paste0(plot_path, 'HM.top15.DE.post-sex_filt.png'), height = 75, width = 100, units = 'cm', res = 500)
TenxPheatmap(data = sex_filt_data, metadata = c("seurat_clusters", "orig.ident"), custom_order_column = "seurat_clusters",
              custom_order = cluster_order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters", assay = "RNA")
graphics.off()