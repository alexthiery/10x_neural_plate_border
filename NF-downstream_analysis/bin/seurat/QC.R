#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
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
}

# Set paths and load data
{
  if (opt$runtype == "user"){
    
    # work dir is base of repo
      
    input_path = "./output/NF-scRNAseq_alignment/" #through nf this is read in using the samplesheet.csv
    output_path = "./results/scRNAseq/seurat_QC/"
    dir.create(output_path, recursive = T)

    sapply(list.files('./NF-downstream_analysis_test/bin/custom_functions/', full.names = T), source)


  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    input_path = "./input/"
    output_path = "./"

    sapply(list.files(opt$custom_functions, full.names = T), source)

  }
  
}

library(Seurat)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
#library(rlist) #dont have this package
library(dplyr)
library(pheatmap)
library(gridExtra)
library(grid)
library(reshape2)
library(viridis)
library(tidyr)

#setwd("/home/rstudio/NF-downstream_analysis")

files <- list.files(input_path, recursive = T, full.names = T)
# remove file suffix
file.path <- dirname(files)[!duplicated(dirname(files))]
# make dataframe with stage matching directory 
### NEED TO MAKE THIS BIT GENERIC

### Non-generic version:
sample = c("THI300A1" = "hh4-1", "THI300A3" = "ss4-1")
matches <- sapply(names(sample), function(x) file.path[grep(pattern = x, x = file.path)])
sample.paths <- data.frame(row.names = sample, sample = sample, stage = names(matches), path = matches, run = gsub(".*-", "", sample))
seurat <- apply(sample.paths, 1, function(x) CreateSeuratObject(counts= Read10X(data.dir = x[["path"]]), project = x[["sample"]]))
seurat <- merge(x = seurat[[1]], y=seurat[-1], add.cell.ids = names(seurat), project = "chick.10x")

# store mitochondrial percentage in object meta data
seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")



######################################### Filtering  ##################################################################

seurat_filtered_low <- subset(seurat, subset = c(nFeature_RNA > 600 & nFeature_RNA < 8000 & percent.mt < 30))
seurat_filtered_med <- subset(seurat, subset = c(nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 15))
seurat_filtered_high <- subset(seurat, subset = c(nFeature_RNA > 1500 & nFeature_RNA < 10000 & percent.mt < 8))

seurat_filtered_low@meta.data["filtering"] <- "low"
seurat_filtered_med@meta.data["filtering"] <- "med"
seurat_filtered_high@meta.data["filtering"] <- "high"
seurat@meta.data["filtering"] <- "unfiltered"
seurat_list <- list(seurat, seurat_filtered_low, seurat_filtered_med, seurat_filtered_high)

######################################### QC plots  ##################################################################

############### Violin plots of QC stats ###############

extract_md <- function(x){
  seurat_meta <- x@meta.data
  md <- seurat_meta[, c("nCount_RNA", "nFeature_RNA", "percent.mt", "filtering", "orig.ident")]
  rownames(md) <- c()
  return(md)
}
meta_data_list <- lapply(seurat_list, extract_md)
meta_data_df <- do.call("rbind", meta_data_list)
meta_data_long <- gather(meta_data_df, QC_metric, value, nCount_RNA:percent.mt, factor_key = TRUE)
meta_data_long$filtering <- factor(meta_data_long$filtering)
meta_data_long$filtering <- factor(meta_data_long$filtering, levels = c("unfiltered", "low", "med", "high"))

nCount_RNA <-  filter(meta_data_long, QC_metric == "nCount_RNA")
png(paste0(output_path, "nCount_RNA_violin.png"), width=30, height=20, units = 'cm', res = 200)
ggplot(nCount_RNA, aes(x = filtering, y = value, fill = filtering)) + 
  geom_violin(trim = TRUE) + 
  scale_fill_viridis(discrete = T) +
  ggtitle("Number of RNA transcripts after filtering") +
  xlab("Filtering threshold") +
  ylab("Number of RNA transcripts") +
  theme_classic() +
  theme(legend.position = "none")
graphics.off()

nFeature_RNA <-  filter(meta_data_long, QC_metric == "nFeature_RNA")
png(paste0(output_path, "nFeature_RNA_violin.png"), width=30, height=20, units = 'cm', res = 200)
ggplot(nFeature_RNA, aes(x = filtering, y = value, fill = filtering)) + 
  geom_violin() + 
  scale_fill_viridis(discrete = T) +
  ggtitle("Number of genes after filtering") +
  xlab("Filtering threshold") +
  ylab("Number of genes") +
  theme_classic() +
  theme(legend.position = "none")
graphics.off()

percent.mt <-  filter(meta_data_long, QC_metric == "percent.mt")
png(paste0(output_path, "percent.mt_violin.png"), width=30, height=20, units = 'cm', res = 200)
ggplot(percent.mt, aes(x = filtering, y = value, fill = filtering)) + 
  geom_violin() + 
  scale_fill_viridis(discrete = T) +
  ggtitle("Percentage of genes which are mitochondrial after filtering") +
  xlab("Filtering threshold") +
  ylab("% mitochondrial genes") +
  theme_classic() +
  theme(legend.position = "none")
graphics.off()

## can facet wrap all the violin plots but because of different scales doesnt look good
# ggplot(data = meta_data_long, aes(filtering, value)) +
#   geom_violin() + 
#   facet_wrap(~ QC_metric)

## also have stage information so can put that in if of interest

############### Bar charts of cell counts ###############

extract_cell_count <- function(x){
  m <- as.data.frame(table(x$orig.ident))
  m$filtering <- x@meta.data$filtering[1]
  return(m)
}
cell_counts_list <- lapply(seurat_list, extract_cell_count)
cell_counts_df <- do.call("rbind", cell_counts_list)
cell_counts_df$filtering <- factor(cell_counts_df$filtering)
cell_counts_df$filtering <- factor(cell_counts_df$filtering, levels = c("unfiltered", "low", "med", "high"))

png(paste0(output_path,"cell_counts_barchart.png"), width=30, height=20, units = 'cm', res = 200)
ggplot(cell_counts_df, aes(x = filtering, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis(discrete = T) +
  ggtitle("Cell counts after filtering") +
  xlab("Filtering threshold") +
  ylab("Cell Count") +
  theme_classic()
graphics.off()


############### Mean gene counts per cell ###############
means <- aggregate(nFeature_RNA, by = list(nFeature_RNA$filtering, nFeature_RNA$orig.ident), FUN = mean)

png(paste0(output_path, "mean_gene_counts_linechart.png"), width=30, height=20, units = 'cm', res = 200)
means %>%
  tail(10) %>%
  ggplot( aes(x=Group.1, y=value, color=Group.2)) +
  geom_line(aes(group = Group.2)) +
  geom_point() +
  scale_color_viridis(discrete = T) +
  ggtitle("Mean gene counts per cell after filtering") +
  xlab("Filtering threshold") +
  ylab("Mean gene count") +
  theme_classic() +
  theme(legend.position = "none")
graphics.off()


## cant get this to work for medians too??
#medians <- aggregate(nFeature_RNA, by = list(nFeature_RNA$filtering, nFeature_RNA$orig.ident), FUN = median)



##############    Scaling, PCA and clustering     #####################

# PCA
seurat_process_PCA <- function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x)
  x <- RunPCA(x, feature = VariableFeatures(object = x))
  return(x)
}

seurat_list_PCAs <- lapply(seurat_list, seurat_process_PCA)

# PCA plots
PCA1 <- DimPlot(seurat_list_PCAs[[1]], reduction = "pca") +
  theme(legend.position = "none") +
  ggtitle("Unfiltered")
PCA2 <- DimPlot(seurat_list_PCAs[[2]], reduction = "pca") +
  theme(legend.position = "none") +
  ggtitle("Low")
PCA3 <- DimPlot(seurat_list_PCAs[[3]], reduction = "pca") +
  theme(legend.position = "none") +
  ggtitle("Med")
PCA4 <- DimPlot(seurat_list_PCAs[[4]], reduction = "pca") +
  theme(legend.position = "none") +
  ggtitle("High")

gPCA <- DimPlot(seurat_list_PCAs[[4]], reduction = "pca")
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(gPCA)

png(paste0(output_path, "PCA_plots.png"), width=30, height=30, units = 'cm', res = 200)
grid.arrange(arrangeGrob(PCA1, PCA2, PCA3, PCA4, nrow=2), mylegend, nrow = 2,
             heights=c(13, 1))
graphics.off()

# PCA loadings
Load1 <- VizDimLoadings(seurat_list_PCAs[[1]], dims = 1, reduction = "pca") +
  ggtitle("Unfiltered")
Load2 <- VizDimLoadings(seurat_list_PCAs[[2]], dims = 1, reduction = "pca") +
  ggtitle("Low")
Load3 <- VizDimLoadings(seurat_list_PCAs[[3]], dims = 1, reduction = "pca") +
  ggtitle("Med")
Load4 <- VizDimLoadings(seurat_list_PCAs[[4]], dims = 1, reduction = "pca") +
  ggtitle("High")

png(paste0(output_path,"PCA_loadings.png"), width=30, height=30, units = 'cm', res = 200)
grid.arrange(Load1, Load2, Load3, Load4, nrow = 2)
graphics.off()

# Other plots that might be of interest
#DimHeatmap(seurat, dims = 1:10, cells = 500, balanced = TRUE)
#ElbowPlot(seurat)

# Clustering and differential expression
seurat_process_cluster <- function(x){
  x <- FindNeighbors(x, dims = 1:10)
  x <- FindClusters(x, resolution = 0.5)
  x <- RunUMAP(x, dims = 1:10)
  return(x)
}

seurat_list_clusters <- lapply(seurat_list_PCAs, seurat_process_cluster)

#might want some UMAPs before get onto the heatmaps
UMAP1 <- DimPlot(seurat_list_clusters[[1]], reduction = "umap") +
  ggtitle("Unfiltered")
UMAP2 <- DimPlot(seurat_list_clusters[[2]], reduction = "umap") +
  ggtitle("Low")
UMAP3 <- DimPlot(seurat_list_clusters[[3]], reduction = "umap") +
  ggtitle("Med")
UMAP4 <- DimPlot(seurat_list_clusters[[4]], reduction = "umap") +
  ggtitle("High")

png(paste0(output_path, "UMAP_plots.png"), width=30, height=30, units = 'cm', res = 200)
grid.arrange(UMAP1, UMAP2, UMAP3, UMAP4, nrow = 2)
graphics.off()

# Differential expression
seurat_process_top_10 <- function(x){
  seurat.markers <- FindAllMarkers(x, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  top10 <- seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  return(top10)
  }
 
top10_list <- lapply(seurat_list_clusters, seurat_process_top_10)
 
#DoHeatmap(seurat_list_clusters[[1]], features = top10_list[[1]]$gene) + NoLegend()

#GM.plot(seurat_list_clusters[[1]], features = top10_list[[1]]$gene)
 
 
#tenx.pheatmap(seurat_list_clusters[[1]], metadata = "orig.ident", selected_genes = top10_list[[1]]$gene)



