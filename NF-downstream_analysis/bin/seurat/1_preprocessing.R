#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
library(Seurat)
library(STACAS)
library(future)
library(tidyverse)
library(grid)
library(gridExtra)
library(clustree)
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
    
    plot_path = "./output/NF-downstream_analysis_stacas/seurat/1_preprocessing/plots/"
    rds_path = "./output/NF-downstream_analysis_stacas/seurat/1_preprocessing/rds_files/"
    data_path = "./output/NF-scRNAseq_alignment/cellranger/count/filtered_feature_bc_matrix"
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/filtered_feature_bc_matrix"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multiprocess", workers = ncores)
    options(future.globals.maxSize = 16* 1024^3) # 16gb
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

# Make dataframe with stage and replicate info extracted from path
input <- list.dirs(data_path, recursive = FALSE, full.names = TRUE)
input <- data.frame(sample = sub('.*/', '', input), run = str_split(sub('.*/', '', input), pattern = "-", simplify = T)[,2], path = input)

# Init list of seurat objects then merge
seurat_list <- apply(input, 1, function(x) CreateSeuratObject(counts= Read10X(data.dir = x[["path"]]), project = x[["sample"]]))
names(seurat_list) <- input$sample
seurat_all <- merge(x = seurat_list[[1]], y=seurat_list[-1], add.cell.ids = names(seurat_list), project = "chick.10x")

# Add metadata col for seq run
seurat_all@meta.data[["run"]] <- gsub(".*-", "", as.character(seurat_all@meta.data$orig.ident))
seurat_all@meta.data[["stage"]] <- gsub("-.*", "", as.character(seurat_all@meta.data$orig.ident))

# Convert metadata character cols to factors
seurat_all@meta.data[sapply(seurat_all@meta.data, is.character)] <- lapply(seurat_all@meta.data[sapply(seurat_all@meta.data, is.character)], as.factor)

# Make seurat gene annotation dataframe and save
annotations <- read.table(paste0(input[1,'path'], '/features.tsv.gz'), col.names = c('Accession', 'Gene', 'V3', 'V4'))[,1:2]
# make gene names unique in annotations dataframe in order to match seurat annotations
annotations$Gene <- make.unique(annotations$Gene)
# Save annnotation dataframe
write.table(annotations, 'seurat_annotations.csv', row.names=FALSE, quote=FALSE, sep=',')

# Remove genes expressed in fewer than 5 cells
seurat_all <- DietSeurat(seurat_all, features = names(which(Matrix::rowSums(GetAssayData(seurat_all) > 0) >=5)))

# Store mitochondrial percentage in object meta data
seurat_all <- PercentageFeatureSet(seurat_all, pattern = "^MT-", col.name = "percent.mt")

# make dataframe with different filtering parameters which can be put into a loop for carrying out downstream analysis
filter_thresholds <- data.frame(gene_min = c(0, 1000, 1500, 2000), gene_max = c(Inf, 7000, 6500, 6000), MT_max = c(Inf, 15, 15, 15), row.names = c("unfilt", "low", "med", "high"))

# Plot remaining cells following different filter thresholds
filter_qc <- lapply(rownames(filter_thresholds), function(condition){
  seurat_all@meta.data %>%
    filter(nFeature_RNA > filter_thresholds[condition,'gene_min']) %>%
    filter(nFeature_RNA < filter_thresholds[condition,'gene_max']) %>%
    filter(percent.mt < filter_thresholds[condition,'MT_max']) %>%
    group_by(orig.ident) %>%
    tally() %>%
    rename(!!condition := n)
})

filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc)

png(paste0(plot_path, 'remaining_cell_table.png'), height = 10, width = 18, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(filter_qc, rows=NULL, theme = ttheme_minimal()))
graphics.off()

png(paste0(plot_path, 'remaining_cell_bar.png'), height = 15, width = 21, units = 'cm', res = 400)
ggplot(filter_qc %>% reshape2::melt(), aes(x=variable, y=value, fill=orig.ident)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  xlab("Filter Condition") +
  ylab("Cell Count") +
  ggtitle("Cell count after filtering") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
graphics.off()


# Plot median gene count per cell following different filter thresholds
filter_qc <- lapply(rownames(filter_thresholds), function(condition){
  seurat_all@meta.data %>%
    filter(nFeature_RNA > filter_thresholds[condition,'gene_min']) %>%
    filter(nFeature_RNA < filter_thresholds[condition,'gene_max']) %>%
    filter(percent.mt < filter_thresholds[condition,'MT_max']) %>%
    group_by(orig.ident) %>%
    summarise(median = median(nFeature_RNA, na.rm = TRUE)) %>%
    mutate(median = as.integer(median)) %>%
    rename(!!condition := median)
})

filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc)

png(paste0(plot_path, 'median_gene_count_table.png'), height = 10, width = 18, units = 'cm', res = 400)
grid.arrange(top=textGrob("Median Gene Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(filter_qc, rows=NULL, theme = ttheme_minimal()))
graphics.off()

png(paste0(plot_path, 'median_gene_count_bar.png'), height = 15, width = 21, units = 'cm', res = 400)
ggplot(filter_qc %>% reshape2::melt(), aes(x=variable, y=value, group=orig.ident)) +
  geom_line(aes(colour = orig.ident)) +
  xlab("Filter Condition") +
  ylab("Median Gene Count") +
  ggtitle("Median gene count per cell after filtering") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
graphics.off()


# Plot median gene count on simulated minimum filter thresholds
tests <- data.frame(gene_cutoff = seq(from = 0, to = 3000, by = 10))

filter_qc <- lapply(seq(from = 0, to = 3000, by = 10), function(cutoff){
  seurat_all@meta.data %>%
    filter(nFeature_RNA > cutoff) %>%
    group_by(orig.ident) %>%
    summarise(median = median(nFeature_RNA, na.rm = TRUE)) %>%
    mutate(median = as.integer(median)) %>%
    rename(!! paste(cutoff) := median)
})

filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc) %>% reshape2::melt() %>% mutate(variable = as.integer(variable)*10)

png(paste0(plot_path, 'median_gene_count_bar.png'), height = 15, width = 21, units = 'cm', res = 400)
ggplot(filter_qc, aes(x=variable, y=value, group=orig.ident)) +
  geom_line(aes(colour = orig.ident)) +
  xlab("Lower Gene Threshold") +
  ylab("Median Gene Count") +
  ggtitle("Median gene counts at simulated minimum filter thresholds") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
graphics.off()


# Plot violins for nCount, nFeature and percent.mt at different filtering thresholds
filter_qc <- lapply(rownames(filter_thresholds), function(condition){
  seurat_all@meta.data %>%
    filter(nFeature_RNA > filter_thresholds[condition,'gene_min']) %>%
    filter(nFeature_RNA < filter_thresholds[condition,'gene_max']) %>%
    filter(percent.mt < filter_thresholds[condition,'MT_max']) %>%
    dplyr::select(orig.ident, nCount_RNA, nFeature_RNA, percent.mt) %>%
    mutate(filter_condition = !!condition)
})

filter_qc <- do.call(rbind, filter_qc) %>%
  mutate(filter_condition = factor(filter_condition, rownames(filter_thresholds))) %>%
  mutate(nCount_RNA = ifelse(nCount_RNA >= 100000, 100000, nCount_RNA)) %>% # limit max RNA to 100k
  reshape2::melt()

png(paste0(plot_path, 'violins_filter_thresholds.png'), height = 18, width = 30, units = 'cm', res = 400)
ggplot(filter_qc, aes(x = filter_condition, y = value, fill = orig.ident)) +
  geom_violin() +
  facet_wrap(~ variable, nrow = 3, strip.position = "left", scales = "free_y") +
  theme_minimal() +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(strip.placement = "outside")
graphics.off()

##########################################################################
############# Remove data which do not pass filter threshold #############

seurat_all <- subset(seurat_all, subset = c(nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 15))

# Split object by run and find integration points
seurat_split <- SplitObject(seurat_all, split.by = "run")

# Log normalize data, find variable features and scale data
seurat_split <- lapply(seurat_split, NormalizeData, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_split <- lapply(seurat_split, FindVariableFeatures, selection.method = "vst", nfeatures = 2000)
seurat_split <- lapply(seurat_split, function(x) ScaleData(x, features = rownames(x), vars.to.regress = "percent.mt"))

# Run PCA analysis
seurat_split <- lapply(seurat_split, RunPCA)


# Plot heatmap of top variable genes across top principle components
for(run in names(seurat_split)){
  png(paste0(plot_path, "dimHM_run_", run, ".png"), width=30, height=65, units = 'cm', res = 200)
  DimHeatmap(seurat_split[[run]], dims = 1:40, balanced = TRUE, cells = 500)
  graphics.off()
}


# In order to calculate the PC cutoff for downstream analysis, we unbiasedly calculate the point at which PC elbow starts.
# First we take the larger value of the point where the principal components only contribute 5% of standard deviation and the point where the principal components cumulatively contribute 90% of the standard deviation.
# Next we take the point where the percent change in variation between the consecutive PCs is less than 0.1%.
# The smaller out of these two values is determined at the elbow cutoff
for(run in names(seurat_split)){
  png(paste0(plot_path, "ElbowCutoff_run", run, ".png"), width=30, height=20, units = 'cm', res = 200)
  print(ElbowCutoff(seurat_split[[run]], return = 'plot'))
  graphics.off()
}

pc_cutoff <- lapply(seurat_split, ElbowCutoff)

# Run clustering and UMAP at different PCA cutoffs - save this output to compare the optimal number of PCs to be used
for(run in names(seurat_split)){
  png(paste0(plot_path, "UMAP_PCA_comparison_run_", run, ".png"), width=40, height=30, units = 'cm', res = 200)
  print(PCALevelComparison(seurat_split[[run]], PCA_levels = c(pc_cutoff[[run]]-5, pc_cutoff[[run]], pc_cutoff[[run]]+5, pc_cutoff[[run]]+10), cluster_res = 1))
  graphics.off()
}

# Use clustering resolution = 0.5 for filtering
seurat_split <- lapply(names(pc_cutoff), function(x) RunUMAP(seurat_split[[x]], dims = 1:pc_cutoff[[x]], verbose = FALSE))
names(seurat_split) <- names(pc_cutoff)

seurat_split <- lapply(names(pc_cutoff), function(x) FindNeighbors(seurat_split[[x]], dims = 1:pc_cutoff[[x]], verbose = FALSE))
names(seurat_split) <- names(pc_cutoff)

# Find optimal cluster resolution
for(run in names(seurat_split)){
  png(paste0(plot_path, "clustree_run_", run, ".png"), width=70, height=35, units = 'cm', res = 200)
  ClustRes(seurat_object = seurat_split[[run]], by = 0.2)
  graphics.off()
}

############################## Identify poor quality clusters #######################################

# Use higher cluster resolution for filtering poor clusters
seurat_split <- lapply(seurat_split, FindClusters, resolution = 1, verbose = FALSE)

# Plot UMAP for clusters and developmental stage
for(run in names(seurat_split)){
  png(paste0(plot_path, "UMAP_run_", run, ".png"), width=40, height=20, units = 'cm', res = 200)
  ClustStagePlot(seurat_split[[run]])
  graphics.off()
}

# Plot QC for each cluster
for(run in names(seurat_split)){
  png(paste0(plot_path, "QCPlot_run_", run, ".png"), width=32, height=28, units = 'cm', res = 200)
  QCPlot(seurat_split[[run]], quantiles = c(0.25, 0.75))
  graphics.off()
}

# Automatically find poor quality clusters
poor_clusters <- lapply(seurat_split, IdentifyOutliers, metrics = c('nCount_RNA', 'nFeature_RNA'), quantiles = c(0.25, 0.75))

# Plot UMAP for poor quality clusters
for(run in names(seurat_split)){
  png(paste0(plot_path, "PoorClusters_run_", run, ".png"), width=60, height=20, units = 'cm', res = 200)
  ClusterDimplot(seurat_split[[run]], clusters = poor_clusters[[run]], plot_title = 'poor quality clusters')
  graphics.off()
}

# Filter poor quality clusters
preprocessing_data <- lapply(names(seurat_split), function(x) {
  subset(seurat_split[[x]], cells = rownames(filter(seurat_split[[x]]@meta.data, seurat_clusters %in% poor_clusters[[x]])), invert = T)
})
names(preprocessing_data) <- names(seurat_split)

saveRDS(preprocessing_data, paste0(rds_path, "preprocessing_data.RDS"), compress = FALSE)