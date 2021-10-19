
coexpression_highlight_cells = function(seurat_object, gm_1, gm_2, bin_number = 10){
  # check gms are not null and are arrays
  if(sum(is.null(gm_1), is.null(gm_2)) != 0){
    stop("one or more gene modules lists are empty")
  }
  if(class(gm_1 ) != 'character' | class(gm_1 ) != 'character'){
    stop("one or more gene modules lists are not characters")
  }
  # calculate expression aggregates and products per cell and use that to order them
  x = seurat_object@meta.data[,'scHelper_cell_type']
  gm_1_sum <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data'))[,gm_1] %>% rowSums(.)
  gm_2_sum <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data'))[,gm_2] %>% rowSums(.)
  plot_data <- data.frame(gm_1_sum = gm_1_sum, gm_2_sum = gm_2_sum, x = x)
  plot_data <- plot_data %>% mutate(ratio = gm_1_sum/(gm_1_sum + gm_2_sum)) %>% arrange(ratio)
  ordered_cells <- rownames(plot_data)
  # create dimplots of cells in each of the bins
  plots <- c()
  bin_size <- length(ordered_cells)/bin_number
  for (i in 0:(bin_number-1)){
    start <- bin_size*i + 1
    end <- (bin_size * (i+1))-1
    bin <- ordered_cells[start:end]
    plots[[i+1]] <- DimPlot(seurat_object, cells.highlight = bin, combine = TRUE) + NoLegend() + ggtitle(paste("Bin ", i+1))
  }
  return(grid.arrange(grobs = plots))
}


coexpression = function(seurat_object, gm_1, gm_2, meta_data='scHelper_cell_type', show_bins = FALSE, bin_number = 10, binarise = FALSE){
  
  x = seurat_object@meta.data[,meta_data,drop=FALSE]
  
  if(binarise == TRUE){
    gm_1_sum <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'counts'))[,gm_1]
    gm_1_sum[gm_1_sum>0] <- 1
    gm_1_sum <- rowSums(gm_1_sum)
    
    gm_2_sum <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'counts'))[,gm_2]
    gm_2_sum[gm_2_sum>0] <- 1
    gm_2_sum <- rowSums(gm_2_sum)
  } else {
    gm_1_sum <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data'))[,gm_1] %>% rowSums(.)
    gm_2_sum <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data'))[,gm_2] %>% rowSums(.)
  }

  
  gm_product <- gm_1_sum * gm_2_sum
  
  plot_data <- data.frame(gm_1_sum = gm_1_sum, gm_2_sum = gm_2_sum, gm_product = gm_product, x)
  
  plot_data <- plot_data[!(plot_data$gm_1_sum == 0 & plot_data$gm_2_sum == 0),]
  
  plot_data <- plot_data %>% mutate(ratio = gm_1_sum/(gm_1_sum + gm_2_sum)) %>% arrange(ratio)
  
  plot_data$x <- 1:nrow(plot_data)
  
  p1 = ggplot(plot_data, aes(x = x, y = gm_1_sum)) +
    geom_smooth(method = "gam", se = FALSE)
  
  p2 = ggplot(plot_data, aes(x = x, y = gm_2_sum)) +
    geom_smooth(method = "gam", se = FALSE)
  
  p3 = ggplot(plot_data, aes(x = x, y = gm_product)) +
    geom_smooth(method = "gam", se = FALSE)
  
  p4 = ggplot(plot_data, aes(x = x, y = gm_product, colour = !!sym(meta_data))) +
    geom_point()
  
  if (show_bins == TRUE){
    bin_size <- length(plot_data$x)/bin_number
    start = c()
    end = c()
    for (i in 0:(bin_number-1)){
      start <- c(start, bin_size*i + 1)
      end <- c(end, bin_size * (i+1))
    }
    bins <- data.frame(start = start, end = end)
    bins <- bins %>% mutate(x = ifelse(row_number() %% 2 == 0, 'even', 'odd'))
    p1 <- p1 +
      geom_rect(data = bins, inherit.aes = FALSE, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = x)) +
      theme(legend.position = "none") +
      scale_fill_manual(values = alpha(c("green", "skyblue"), alpha = 0.07))
    p2 <- p2 +
      geom_rect(data = bins, inherit.aes = FALSE, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = x)) +
      theme(legend.position = "none") +
      scale_fill_manual(values = alpha(c("green", "skyblue"), alpha = 0.07))
    p3 <- p3 +
      geom_rect(data = bins, inherit.aes = FALSE, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = x)) +
      theme(legend.position = "none") +
      scale_fill_manual(values = alpha(c("green", "skyblue"), alpha = 0.07))
  }
  
  return(p1/p2/p3/p4)
}

subset_seurat = function(seurat_data, population, split_by = "scHelper_cell_type", rerun_UMAP = FALSE){
  subset <- subset(seurat_data, cells = rownames(seurat_data@meta.data)[seurat_data@meta.data[[split_by]] %in% population])
  Idents(subset) <- split_by
  if (rerun_UMAP == TRUE){
    subset <- RunPCA(object = subset, verbose = FALSE)
    pc_cutoff <- ElbowCutoff(subset)
    subset <- FindNeighbors(subset, dim = 1:pc_cutoff)
    subset <- RunUMAP(subset, dims = 1:pc_cutoff)
  }
  return(subset)
}

filter_genes = function(seurat_object, gene_list, ident_1, ident_2 = NULL, group_by, logfc = 0.5){
  markers <- FindMarkers(seurat_object, group.by = group_by, ident.1 = ident_1, ident.2 = ident_2, only.pos = TRUE, logfc.threshold = logfc)
  # markers <- markers[markers$avg_log2FC > 0, ]
  filt_genes <- gene_list[gene_list %in% rownames(markers)]
  print(setdiff(gene_list, filt_genes))
  return(filt_genes)
}


extract_bin <- function(seurat_object, gm_1, gm_2, meta_data='scHelper_cell_type', bin_number = 10, bin_extract){
  if(sum(is.null(gm_1), is.null(gm_2)) != 0){
    stop("one or more gene modules lists are empty")
  }
  if(class(gm_1 ) != 'character' | class(gm_1 ) != 'character'){
    stop("one or more gene modules lists are not characters")
  }
  # calculate expression aggregates and products per cell and use that to order them
  x = seurat_object@meta.data[,'scHelper_cell_type']
  gm_1_sum <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data'))[,gm_1] %>% rowSums(.)
  gm_2_sum <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data'))[,gm_2] %>% rowSums(.)
  plot_data <- data.frame(gm_1_sum = gm_1_sum, gm_2_sum = gm_2_sum, x = x)
  plot_data <- plot_data %>% mutate(ratio = gm_1_sum/(gm_1_sum + gm_2_sum)) %>% arrange(ratio)
  ordered_cells <- rownames(plot_data)
  
  bin_size <- length(ordered_cells)/bin_number
  start = bin_size*(bin_extract-1)
  end = start + (bin_size-1)
  return(ordered_cells[start:end])
}


library(optparse)
library(future)
library(Seurat)
library(pheatmap)
library(tidyverse)
library(Antler)
library(RColorBrewer)
library(scHelper)
install.packages('patchwork')
library(patchwork)

# Set paths and load data
plot_path = "./plots/"
rds_path = "./rds_files/"
data_path = "./input/"
ncores = opt$cores

dir.create(plot_path, recursive = T)
dir.create(rds_path, recursive = T)

# metadata <- read.csv(list.files(data_path, pattern = "*.csv", full.names = TRUE))
# # metadata <- read.csv('./output/NF-downstream_analysis_stacas/filtered_seurat/cellrank/NF-scRNAseq_alignment_out_metadata.csv')

# seurat_data <- readRDS(list.files(data_path, pattern = "*.RDS", full.names = TRUE)[!list.files(data_path, pattern = "*.RDS") %>% grepl('antler', .)])
seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/transfer_labels/seurat/rds_files/seurat_label_transfer.RDS')
hh7 <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/hh7_splitstage_data/seurat/stage_state_classification/rds_files/hh7_cell_state_classification.RDS')
ss8 <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/ss8_splitstage_data/seurat/stage_state_classification/rds_files/ss8_cell_state_classification.RDS')
ss4 <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/ss4_splitstage_data/seurat/stage_state_classification/rds_files/ss4_cell_state_classification.RDS')

# load antler data
# antler_data <- readRDS(list.files(data_path, pattern = "antler_out.RDS", full.names = TRUE))



# Rename cell states based on group var
# new_cell_states = c('aPPR' = 'PPR', 'PPR' = 'PPR', 'pPPR' = 'PPR', 'aNPB' = 'NPB', 'pNPB' = 'NPB', 'NC' = 'NC', 'delaminating_NC' = 'NC')
# seurat_data$scHelper_cell_type = new_cell_states[as.character(seurat_data$scHelper_cell_type)]

# Set RNA to default assay for plotting expression data
DefaultAssay(seurat_data) <- "integrated"

# Subset cell states and stages from full dataset
subset <- subset_seurat(seurat_data, population = c('hh7', 'ss4', 'ss8'), split_by = 'stage', rerun_UMAP = FALSE)
DimPlot(subset, group.by = 'scHelper_cell_type')
subset <- subset_seurat(subset, population = c('pPPR', 'aPPR', 'NC', 'delaminating_NC', 'PPR', 'aNPB', 'pNPB'), split_by = 'scHelper_cell_type', rerun_UMAP = TRUE)
DimPlot(subset, group.by = 'scHelper_cell_type')


# Set RNA to default assay for plotting expression data
DefaultAssay(seurat_data) <- "RNA"

####################################################################################################
# Run coexpression using PPR and NC gene modules from ss4
antler_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/ss4_splitstage_data/antler/stage_gene_modules/rds_files/antler_out.RDS')

ppr_gm <- antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content$GM23
nc_gm <- antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content$GM29

coexpression_highlight_cells(subset, gm_1 = ppr_gm, gm_2 = nc_gm, bin_number = 20)
coexpression(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), show_bins = TRUE, bin_number = 20)

bin <- extract_bin(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = 20, bin_extract = 6)
DimPlot(ss4, cells.highlight = bin) + NoLegend()
DimPlot(ss8, cells.highlight = bin) + NoLegend()
DimPlot(hh7, cells.highlight = bin) + NoLegend()

ppr_gm <- filter_genes(subset, gene_list = ppr_gm, ident_1 = c('pPPR', 'aPPR', 'PPR'), group_by = 'scHelper_cell_type', logfc = 0.5)
nc_gm <- filter_genes(subset, gene_list = nc_gm, ident_1 = c('NC', 'delaminating_NC'), group_by = 'scHelper_cell_type', logfc = 0.5)

coexpression_highlight_cells(subset, gm_1 = ppr_gm, gm_2 = nc_gm, bin_number = 20)
coexpression(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), show_bins = TRUE, bin_number = 20)

bin <- extract_bin(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = 20, bin_extract = 6)
DimPlot(hh7, cells.highlight = bin) + NoLegend()
DimPlot(ss4, cells.highlight = bin) + NoLegend()
DimPlot(ss8, cells.highlight = bin) + NoLegend()

# Binarise expression before running coexpression as above
coexpression_highlight_cells(subset, gm_1 = ppr_gm, gm_2 = nc_gm, bin_number = 20)
coexpression(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), show_bins = TRUE, bin_number = 20)
bin <- extract_bin(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = 20, bin_extract = 7)
DimPlot(ss4, cells.highlight = bin) + NoLegend()
DimPlot(ss8, cells.highlight = bin) + NoLegend()
DimPlot(hh7, cells.highlight = bin) + NoLegend()

ppr_gm <- filter_genes(subset, gene_list = ppr_gm, ident_1 = c('pPPR', 'aPPR', 'PPR'), group_by = 'scHelper_cell_type', logfc = 0.5)
nc_gm <- filter_genes(subset, gene_list = nc_gm, ident_1 = c('NC', 'delaminating_NC'), group_by = 'scHelper_cell_type', logfc = 0.5)

coexpression_highlight_cells(subset, gm_1 = ppr_gm, gm_2 = nc_gm, bin_number = 20)
coexpression(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), show_bins = TRUE, bin_number = 20)

bin <- extract_bin(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = 10, bin_extract = 4)
DimPlot(hh7, cells.highlight = bin) + NoLegend()
DimPlot(ss4, cells.highlight = bin) + NoLegend()
DimPlot(ss8, cells.highlight = bin) + NoLegend()





start = 0
end = 10
width = end - start




res = (arr - min(arr))/(max(arr) - min(arr)) * width + start











library(tidyverse)
gm_scores <- lapply(list(ppr_gm1 = ppr_gm, nc_gm1 = nc_gm), function(x) t(GetAssayData(object = subset, assay = 'RNA', slot = 'scale.data'))[,x] %>% rowMeans(.)) %>%
  do.call('cbind', .)

subset@meta.data <- merge(seurat_data@meta.data, gm_scores, by = 0) %>% column_to_rownames('Row.names')

FeaturePlot(subset, features = c('ppr_gm1', 'nc_gm1'), blend = TRUE, blend.threshold = 0)


debug(FeaturePlot)
####################################################################################################
# Run coexpression using candidate PPR and NC genes from literature/GMs


PPR = c("SIX1", "DLX3", "HOMER2", "PITX1", "PITX2", 'ENSGALG00000044849', 'PRDM1', 'GATA2', 'GATA3', 'PAX6', 'HIF1A', 'NET1')
NC = c("WNT1", "DRAXIN", "SOX9", "ENSGALG00000030902", "FOXD3", "SOX10", "ETS1", 'WNT6', 'MYC', 'SOX5', 'LMO4', 'SOX8')

DotPlot(subset, features = PPR, group.by = 'scHelper_cell_type')
DotPlot(subset, features = NC, group.by = 'scHelper_cell_type')

coexpression_highlight_cells(subset, gm_1 = PPR, gm_2 = NC, bin_number = 20)
coexpression(subset, gm_1 = PPR, gm_2 = NC, meta_data = c('scHelper_cell_type'), show_bins = TRUE, bin_number = 20)

DimPlot(ss4, cells.highlight = bin) + NoLegend()
DimPlot(ss8, cells.highlight = bin) + NoLegend()
DimPlot(hh7, cells.highlight = bin) + NoLegend()








# PPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "PRDM1")
# NC = c("PAX7", "MSX1", "MSX2", "ETS1", "ENSGALG00000030902", "FOXD3", "TFAP2B")
# 
# PPR = c("SIX1", "DLX3", "ZNF385C", "HOMER2", "DLX6", "PITX1", "PITX2")
# NC = c("PAX7", "WNT1", "DRAXIN", "SOX9", "ENSGALG00000030902", "FOXD3", "TFAP2B")

# ENSGALG00000044849 (FOXI2)
# ENSGALG00000030902 (SNAIL2)













