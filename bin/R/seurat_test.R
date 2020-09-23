#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer",
  'custom_functions', 'm', 2, "character",
  'networkGenes', 'd', 2, "character"
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
    sapply(list.files('./bin/R/custom_functions/', full.names = T), source)
    plot.path = "./results/plots/"
    rds.path = "./results/RDS.files/"
    data.path = "./alignment_out"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    sapply(list.files(opt$custom_functions, full.names = T), source)
    plot.path = "plots/"
    rds.path = "RDS.files/"
    data.path = "."
    
    ncores = opt$cores
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  
  dir.create(plot.path, recursive = T)
  dir.create(rds.path, recursive = T)
  
  # Load packages - packages are stored within renv in the repository
  reticulate::use_python('/usr/bin/python3.7')
  library(Seurat)
  
  library(future)
  library(dplyr)
  library(cowplot)
  library(clustree)
  library(gridExtra)
  library(grid)
  library(pheatmap)
  library(RColorBrewer)
}


# read all files from dir
files <- list.files(data.path, recursive = T, full.names = T)
# remove file suffix
file.path <- dirname(files)[!duplicated(dirname(files))]
# make dataframe with tissue matching directory
sample = c("THI300A1" = "hh4_1", "THI300A3" = "ss4_1", "THI300A4" = "ss8_1", "THI300A6" = "hh6_1",
           "THI725A1" = "hh5_2", "THI725A2" = "hh6_2", "THI725A3" = "hh7_2", "THI725A4" = "ss4_2")
matches <- sapply(names(sample), function(x) file.path[grep(pattern = x, x = file.path)])
sample.paths <- data.frame(row.names = sample, sample = sample, tissue = names(matches), path = matches, run = gsub(".*_", "", sample))

# Make Seurat objects for each of the different samples and then merge
seurat_data <- apply(sample.paths, 1, function(x) CreateSeuratObject(counts= Read10X(data.dir = x[["path"]]), project = x[["sample"]]))
seurat_data <- merge(x = seurat_data[[1]], y=seurat_data[-1], add.cell.ids = gsub("_.*", "", names(seurat_data)), project = "chick.10x")

# Remove genes expressed in fewer than 3 cells
seurat_data <- DietSeurat(seurat_data, features = names(which(Matrix::rowSums(GetAssayData(seurat_data) > 0) >=3)))

# Store mitochondrial percentage in object meta data
seurat_data <- PercentageFeatureSet(seurat_data, pattern = "^MT-", col.name = "percent.mt")

# Remove data which do not pass filter threshold
seurat_data <- subset(seurat_data, subset = c(nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 15))

# Add metadata col for seq run
seurat_data@meta.data[["seq_run"]] <- gsub(".*_", "", as.character(seurat_data@meta.data$orig.ident))

#####################################################################################################
#                           Integrate data from different 10x runs                                  #
#####################################################################################################

# Split object by run and find integration points
seurat_data_integrated <- SplitObject(seurat_data, split.by = "seq_run")

# Log normalize data and find variable features
seurat_data_integrated <- lapply(seurat_data_integrated, function(x) NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000))
seurat_data_integrated <- lapply(seurat_data_integrated, function(x) FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000))

plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 3000 * 1024^2)
seurat_data_integrated <- FindIntegrationAnchors(object.list = seurat_data_integrated, dims = 1:30)
seurat_data_integrated <- IntegrateData(anchorset = seurat_data_integrated, dims = 1:30)

# set inegrated count data as default
DefaultAssay(seurat_data_integrated) <- "integrated"

# Scale data and regress out MT content
# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)
seurat_data_integrated <- ScaleData(seurat_data_integrated, features = rownames(seurat_data_integrated), vars.to.regress = "percent.mt")

# Save RDS after scaling as this step takes time
saveRDS(seurat_data_integrated, paste0(rds.path, "seurat_data_integrated.RDS"))


# Change plot path
curr.plot.path <- paste0(plot.path, '0_seurat_data_integrated/')
dir.create(curr.plot.path)

# Run PCA analysis on the each set of data
seurat_data_integrated <- RunPCA(object = seurat_data_integrated, verbose = FALSE)

#####################################################################################################
#   Seurat's clustering algorithm is based on principle components, so we need to ensure that only the informative PCs are kept!                   #
#####################################################################################################

# Plot heatmap of top variable genes across top principle components
png(paste0(curr.plot.path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(seurat_data_integrated, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

# another heuristic method is ElbowPlot which ranks PCs based on the % variance explained by each PC
png(paste0(curr.plot.path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(seurat_data_integrated, ndims = 40))
graphics.off()

# Run clustering and UMAP at different PCA cutoffs - save this output to compare the optimal number of PCs to be used
png(paste0(curr.plot.path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(seurat_data_integrated, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
# Use clustering resolution = 0.5 for filtering
seurat_data_integrated <- FindNeighbors(seurat_data_integrated, dims = 1:15, verbose = FALSE)
seurat_data_integrated <- RunUMAP(seurat_data_integrated, dims = 1:15, verbose = FALSE)
seurat_data_integrated <- FindClusters(seurat_data_integrated, resolution = 0.5, verbose = FALSE)

# Plot UMAP for clusters and developmental stage
png(paste0(curr.plot.path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(seurat_data_integrated)
graphics.off()

# Plot QC for each cluster
png(paste0(curr.plot.path, "cluster.QC.png"), width=40, height=14, units = 'cm', res = 200)
QC.plot(seurat_data_integrated)
graphics.off()

























#####################################################################################################
#                               Run scaling on non integrated object                                #
#####################################################################################################

# Log normalize data and find variable features
norm.data <- NormalizeData(seurat_data, normalization.method = "LogNormalize", scale.factor = 10000)
norm.data <- FindVariableFeatures(norm.data, selection.method = "vst", nfeatures = 2000)

# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)

# Scale data and regress out MT content
norm.data <- ScaleData(norm.data, features = rownames(norm.data), vars.to.regress = "percent.mt")

# Save RDS after scaling as this step takes time
saveRDS(norm.data, paste0(rds.path, "norm.data.RDS"))
#####################################################################################################
#                    Perform dimensionality reduction by PCA and UMAP embedding                    #
#####################################################################################################

# Read in RDS data if needed
# norm.data <- readRDS(paste0(rds.path, "norm.data.RDS"))

# Change plot path
curr.plot.path <- paste0(plot.path, '0_filt_data/')
dir.create(curr.plot.path)

# Run PCA analysis on the each set of data
norm.data <- RunPCA(object = norm.data, verbose = FALSE)

#####################################################################################################
#   Seurat's clustering algorithm is based on principle components, so we need to ensure that only the informative PCs are kept!                   #
#####################################################################################################

# Plot heatmap of top variable genes across top principle components
png(paste0(curr.plot.path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(norm.data, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

# another heuristic method is ElbowPlot which ranks PCs based on the % variance explained by each PC
png(paste0(curr.plot.path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(norm.data, ndims = 40))
graphics.off()

# Run clustering and UMAP at different PCA cutoffs - save this output to compare the optimal number of PCs to be used
png(paste0(curr.plot.path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(norm.data, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
# Use clustering resolution = 0.5 for filtering
norm.data <- FindNeighbors(norm.data, dims = 1:15, verbose = FALSE)
norm.data <- RunUMAP(norm.data, dims = 1:15, verbose = FALSE)
norm.data <- FindClusters(norm.data, resolution = 0.5, verbose = FALSE)

# Plot UMAP for clusters and developmental stage
png(paste0(curr.plot.path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(norm.data)
graphics.off()

# Plot QC for each cluster
png(paste0(curr.plot.path, "cluster.QC.png"), width=40, height=14, units = 'cm', res = 200)
QC.plot(norm.data)
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(norm.data, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster.order = order.cell.stage.clust(seurat_object = norm.data, col.to.sort = seurat_clusters, sort.by = orig.ident)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) %>% arrange(factor(cluster, levels = cluster.order))

png(paste0(curr.plot.path, 'HM.top15.DE.png'), height = 50, width = 75, units = 'cm', res = 700)
tenx.pheatmap(data = norm.data, metadata = c("seurat_clusters", "orig.ident"), custom_order_column = "seurat_clusters",
              custom_order = cluster.order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters")
graphics.off()

#####################################################################################################
#     Heatmap clearly shows clusters segregate by sex - check this and regress out the sex effect   #
#####################################################################################################

# Change plot path
curr.plot.path <- paste0(plot.path, '1_sex_filt/')
dir.create(curr.plot.path)

# There is a strong sex effect - this plot shows DE genes between clusters 1 and 2 which are preodominantly hh4 clusters. Clustering is driven by sex genes
png(paste0(curr.plot.path, 'HM.top15.DE.pre-sexfilt.png'), height = 40, width = 70, units = 'cm', res = 500)
tenx.pheatmap(data = norm.data[,rownames(norm.data@meta.data[norm.data$seurat_clusters == 1 | norm.data$seurat_clusters == 2,])],
              metadata = c("seurat_clusters", "orig.ident"), selected_genes = rownames(FindMarkers(norm.data, ident.1 = 1, ident.2 = 2)),
              hclust_rows = T, gaps_col = "seurat_clusters")
graphics.off()

# Use W chromosome genes to K-means cluster the cells into male (zz) and female (zw)
W_genes <- as.matrix(norm.data@assays$RNA[grepl("W-", rownames(norm.data@assays$RNA)),])
k_clusters <- kmeans(t(W_genes), 2)
k_clusters <- data.frame(k_clusters$cluster)
norm.data@meta.data$k_clusters <- k_clusters[match(colnames(norm.data@assays$RNA), rownames(k_clusters)),]

# Get rownames for kmeans clusters 1 and 2
k_clus_1 <- rownames(norm.data@meta.data[norm.data@meta.data$k_clusters == 1,])
k_clus_2 <- rownames(norm.data@meta.data[norm.data@meta.data$k_clusters == 2,])

# K clustering identities are stochastic, so I mneed to identify which cluster is male and female
# Sum of W genes is order of magnitude greater in cluster 2 - these are the female cells
sumclus1 <- sum(W_genes[,k_clus_1])
sumclus2 <- sum(W_genes[,k_clus_2])

if(sumclus1 < sumclus2){
  k.male <- k_clus_1
  k.female <- k_clus_2
} else {
  k.female <- k_clus_1
  k.male <- k_clus_2
}

cell.sex.ID <- list("male.cells" = k.male, "female.cells" = k.female)
saveRDS(cell.sex.ID, paste0(rds.path, "sex_kmeans.RDS"))

# Add sex data to meta.data
norm.data@meta.data$sex <- unlist(lapply(rownames(norm.data@meta.data), function(x)
  if(x %in% k.male){"male"} else if(x %in% k.female){"female"} else{stop("cell sex is not assigned")}))


###### Next - subset autosomal genes and Z genes - then calculate the average for each gene for both kmeans clustered cells and plot in order to compare whether the two groups show significant differences in their expression.

# This is a test on autosomal genes to try and calculate and compare FC betweeen clusters
# Calculating median is tricky as there are a lot of dropouts in 10x data so you end up with either 0s (when the median  = 0) or 1 (when the median expression in both clusters is the same - probably a result of normalisation resulting in a UMI of 0 or 1 being normalised to a nominal value)

# Make dataframe for mean Z expression in male cells
mean.Z.male <- data.frame(Z.mean = apply(norm.data@assays$RNA[grepl("Z-", rownames(norm.data@assays$RNA)), k.male], 1, mean))
# add 1 before log2 as log2(1) = 0
mean.Z.male <- log2(mean.Z.male + 1)

# Make dataframe for mean Z expression in female cells
mean.Z.female <- data.frame(Z.mean = apply(norm.data@assays$RNA[grepl("Z-", rownames(norm.data@assays$RNA)), k.female], 1, mean))
mean.Z.female <- log2(mean.Z.female + 1)

# Make dataframe for mean autosomal expression in male cells
mean.auto.male <- data.frame(auto.mean = apply(norm.data@assays$RNA[!grepl("Z-", rownames(norm.data@assays$RNA)) & !grepl("W-", rownames(norm.data@assays$RNA)), k.male], 1, mean))
mean.auto.male <- log2(mean.auto.male + 1)

# Make dataframe for mean autosomal expression in male cells
mean.auto.female <- data.frame(auto.mean = apply(norm.data@assays$RNA[!grepl("Z-", rownames(norm.data@assays$RNA)) & !grepl("W-", rownames(norm.data@assays$RNA)), k.female], 1, mean))
mean.auto.female <- log2(mean.auto.female + 1)

# Calculate FC by subtracting log2 expression from each other
FC <- list()
FC$Z <- mean.Z.male - mean.Z.female
FC$auto <-  mean.auto.male - mean.auto.female

# Plot boxplot of Z gene and autosomal expression in male vs female cells
png(paste0(curr.plot.path,"sex_kmeans_log2FC_boxplot.png"), height = 18, width = 18, units = "cm", res = 200)
boxplot(c(FC$Z, FC$auto),  ylab = "male - female log2 FC (mean normalised UMI +1)", names = c("Z chromosome genes", "autosomal genes"))
graphics.off()

# Z genes are upregulated within male genes relative to female genes whereas autosomal genes have a normal distribution of logFCs

#####################################################################################################
#                                     Regress sex effect                                            #
#####################################################################################################

# Init sexscale object 
norm.data.sexscale <- norm.data

# Re-run findvariablefeatures and scaling
norm.data.sexscale <- FindVariableFeatures(norm.data.sexscale, selection.method = "vst", nfeatures = 2000)
# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)

norm.data.sexscale <- ScaleData(norm.data.sexscale, features = rownames(norm.data.sexscale), vars.to.regress = c("percent.mt", "sex"))

# Save RDS
saveRDS(norm.data.sexscale, paste0(rds.path, "norm.data.sexscale.RDS"))

# Read in RDS data if needed
# norm.data.sexscale <- readRDS(paste0(rds.path, "norm.data.sexscale.RDS"))

# Set plot path
curr.plot.path <- paste0(plot.path, '1_sex_filt/')

# PCA
norm.data.sexscale <- RunPCA(object = norm.data.sexscale, verbose = FALSE)

png(paste0(curr.plot.path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(norm.data.sexscale, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(curr.plot.path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(norm.data.sexscale, ndims = 40))
graphics.off()

png(paste0(curr.plot.path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(norm.data.sexscale, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
norm.data.sexscale <- FindNeighbors(norm.data.sexscale, dims = 1:15, verbose = FALSE)
norm.data.sexscale <- RunUMAP(norm.data.sexscale, dims = 1:15, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(curr.plot.path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
clust.res(seurat.obj = norm.data.sexscale, by = 0.1)
graphics.off()

# Use clustering resolution = 0.5 to look for contamination clusters
norm.data.sexscale <- FindClusters(norm.data.sexscale, resolution = 0.5, verbose = FALSE)

# Plot UMAP for clusters and developmental stage
png(paste0(curr.plot.path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(norm.data.sexscale)
graphics.off()

# Plot QC for each cluster
png(paste0(curr.plot.path, "cluster.QC.png"), width=40, height=14, units = 'cm', res = 200)
QC.plot(norm.data.sexscale)
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(norm.data.sexscale, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster.order = order.cell.stage.clust(seurat_object = norm.data.sexscale, col.to.sort = seurat_clusters, sort.by = orig.ident)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) %>% arrange(factor(cluster, levels = cluster.order))

png(paste0(curr.plot.path, 'HM.top15.DE.post-sexfilt.png'), height = 75, width = 100, units = 'cm', res = 500)
tenx.pheatmap(data = norm.data.sexscale, metadata = c("seurat_clusters", "orig.ident"), custom_order_column = "seurat_clusters",
              custom_order = cluster.order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters")
graphics.off()

#####################################################################################################
#                           Identify and remove contamination (mesoderm and PGCs)                   #
#####################################################################################################

# Change plot path
curr.plot.path <- paste0(plot.path, "2_contamination_filt/")
dir.create(curr.plot.path)

# Identify mesoderm and PGCs
# UMAP plots GOI
genes <- c("EYA2", "SIX1", "TWIST1", "PITX2", "SOX17", "DAZL", "DND1", "CXCR4")

ncol = 4
png(paste0(curr.plot.path, "UMAP_GOI.png"), width = ncol*10, height = 10*ceiling((length(genes)+1)/ncol), units = "cm", res = 200)
multi.feature.plot(seurat.obj = norm.data.sexscale, gene.list = genes, plot.clusters = T,
                   plot.stage = T, label = "", cluster.col = "RNA_snn_res.0.5", n.col = ncol)
graphics.off()

# Dotplot for identifying PGCs, Early mesoderm and Late mesoderm
png(paste0(curr.plot.path, "dotplot.GOI.png"), width = 20, height = 12, units = "cm", res = 200)
DotPlot(norm.data.sexscale, features = c( "SOX17", "CXCR4","EYA2", "TWIST1", "SIX1",  "PITX2", "DAZL"))
graphics.off()
