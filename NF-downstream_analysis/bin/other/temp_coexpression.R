
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
  
  meta_data = seurat_object@meta.data[,meta_data,drop=FALSE]
  
  gm_1_sum <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data'))[,gm_1] %>% rowSums(.)
  gm_2_sum <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data'))[,gm_2] %>% rowSums(.)
  
  plot_data <- data.frame(gm_1_sum = gm_1_sum, gm_2_sum = gm_2_sum, meta_data)
  
  # plot_data <- plot_data[!(plot_data$gm_1_sum == 0 & plot_data$gm_2_sum == 0),]
  
  plot_data <- plot_data %>%
    mutate(ratio = gm_1_sum/(gm_1_sum + gm_2_sum)) %>%
    arrange(ratio) %>%
    pivot_longer(cols = c(gm_1_sum, gm_2_sum)) %>%
    rename(expression_magnitude = value)
  
  plot_data$cell_order <- 1:nrow(plot_data)
  
  p1 = ggplot(plot_data, aes(x = cell_order, y = expression_magnitude)) +
    theme_classic() +
    geom_point(aes(colour = scHelper_cell_type), alpha = 0.5, size = 0.5) +
    geom_smooth(method = "gam", se = FALSE, colour="purple", size=1.5) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size=10),
          strip.text.x = element_text(size = 10)) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    facet_wrap(~name, dir = "v", scales = "free") +
    xlab("Cells") + ylab("Expression level")
  
  if (show_bins == TRUE){
    bin_size = length(plot_data$cell_order) / bin_number
    
    start <- seq(min(plot_data$cell_order), max(plot_data$cell_order), by = bin_size)
    end <- seq(min(plot_data$cell_order + bin_size-1), max(plot_data$cell_order + bin_size-1), by = bin_size)
    bins <- data.frame(start, end) %>% mutate(x = ifelse(row_number() %% 2 == 0, 'even', 'odd'))
    
    p1 <- p1 +
      geom_rect(data = bins, inherit.aes = FALSE, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = x)) +
      scale_fill_manual(values = alpha(c("gray20", "gray80"), alpha = 0.07)) +
      geom_rect(data = data.frame(start = 20, end = 2000), inherit.aes = FALSE, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill='red')
  }
  
  return(p1)
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
# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/transfer_labels/seurat/rds_files/seurat_label_transfer.RDS')
# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/transfer_subset/transfer_npb_subset/seurat/transfer_cluster/rds_files/transfer_clustered_data.RDS')
hh6 <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/hh6_splitstage_data/seurat/stage_state_classification/rds_files/hh6_cell_state_classification.RDS')
hh7 <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/hh7_splitstage_data/seurat/stage_state_classification/rds_files/hh7_cell_state_classification.RDS')
ss8 <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/ss8_splitstage_data/seurat/stage_state_classification/rds_files/ss8_cell_state_classification.RDS')
ss4 <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/ss4_splitstage_data/seurat/stage_state_classification/rds_files/ss4_cell_state_classification.RDS')

# load antler data
# antler_data <- readRDS(list.files(data_path, pattern = "antler_out.RDS", full.names = TRUE))



# Rename cell states based on group var
# new_cell_states = c('aPPR' = 'PPR', 'PPR' = 'PPR', 'pPPR' = 'PPR', 'aNPB' = 'NPB', 'pNPB' = 'NPB', 'NC' = 'NC', 'dNC' = 'NC')
# seurat_data$scHelper_cell_type = new_cell_states[as.character(seurat_data$scHelper_cell_type)]

# # Set RNA to default assay for plotting expression data
# DefaultAssay(seurat_data) <- "integrated"
# # 
# # Subset cell states and stages from full dataset
# subset <- subset_seurat(seurat_data, population = c('hh7', 'ss4', 'ss8'), split_by = 'stage', rerun_UMAP = FALSE)
# DimPlot(subset, group.by = 'scHelper_cell_type')
# subset <- subset_seurat(subset, population = c('pPPR', 'aPPR', 'NC', 'dNC', 'PPR', 'aNPB', 'pNPB'), split_by = 'scHelper_cell_type', rerun_UMAP = TRUE)
# DimPlot(subset, group.by = 'scHelper_cell_type')
# 
# 
# subset <- seurat_data
# DimPlot(subset, group.by = 'seurat_clusters')
# 
# # Set RNA to default assay for plotting expression data
# DefaultAssay(seurat_data) <- "RNA"

####################################################################################################
# # Run coexpression using PPR and NC gene modules from ss4
# antler_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/ss4_splitstage_data/antler/stage_gene_modules/rds_files/antler_out.RDS')
# 
# ppr_gm <- antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content$GM23
# nc_gm <- antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content$GM29
# 
# coexpression(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), show_bins = TRUE, bin_number = 10)
# coexpression_highlight_cells(subset, gm_1 = ppr_gm, gm_2 = nc_gm, bin_number = 10)
# 
# bin <- extract_bin(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = 20, bin_extract = 6)
# DimPlot(ss4, cells.highlight = bin) + NoLegend()
# DimPlot(ss8, cells.highlight = bin) + NoLegend()
# DimPlot(hh7, cells.highlight = bin) + NoLegend()
# 
# ppr_gm <- filter_genes(subset, gene_list = ppr_gm, ident_1 = c('pPPR', 'aPPR', 'PPR'), group_by = 'scHelper_cell_type', logfc = 0.5)
# nc_gm <- filter_genes(subset, gene_list = nc_gm, ident_1 = c('NC', 'dNC'), group_by = 'scHelper_cell_type', logfc = 0.5)
# 
# coexpression_highlight_cells(subset, gm_1 = ppr_gm, gm_2 = nc_gm, bin_number = 20)
# coexpression(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), show_bins = TRUE, bin_number = 20)
# 
# bin <- extract_bin(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = 20, bin_extract = 6)
# DimPlot(hh7, cells.highlight = bin) + NoLegend()
# DimPlot(ss4, cells.highlight = bin) + NoLegend()
# DimPlot(ss8, cells.highlight = bin) + NoLegend()
# 
# # Binarise expression before running coexpression as above
# coexpression_highlight_cells(subset, gm_1 = ppr_gm, gm_2 = nc_gm, bin_number = 20)
# coexpression(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), show_bins = TRUE, bin_number = 20)
# bin <- extract_bin(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = 20, bin_extract = 7)
# DimPlot(ss4, cells.highlight = bin) + NoLegend()
# DimPlot(ss8, cells.highlight = bin) + NoLegend()
# DimPlot(hh7, cells.highlight = bin) + NoLegend()
# 
# ppr_gm <- filter_genes(subset, gene_list = ppr_gm, ident_1 = c('pPPR', 'aPPR', 'PPR'), group_by = 'scHelper_cell_type', logfc = 0.5)
# nc_gm <- filter_genes(subset, gene_list = nc_gm, ident_1 = c('NC', 'dNC'), group_by = 'scHelper_cell_type', logfc = 0.5)
# 
# coexpression_highlight_cells(subset, gm_1 = ppr_gm, gm_2 = nc_gm, bin_number = 20)
# coexpression(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), show_bins = TRUE, bin_number = 20)
# 
# bin <- extract_bin(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = 10, bin_extract = 4)
# DimPlot(hh7, cells.highlight = bin) + NoLegend()
# DimPlot(ss4, cells.highlight = bin) + NoLegend()
# DimPlot(ss8, cells.highlight = bin) + NoLegend()
# 
# 
# 
# 







####################################################################################################
# Run coexpression using PPR and NC gene modules from ss8
subset <- readRDS('./output/NF-downstream_analysis_stacas/transfer_subset/transfer_ppr_nc_subset/seurat/transfer_cluster/rds_files/transfer_clustered_data.RDS')

subset <- subset_seurat(subset, population = c('hh7', 'ss4', 'ss8'), split_by = 'stage', rerun_UMAP = FALSE)
subset <- subset_seurat(subset, population = c('pPPR', 'aPPR', 'NC', 'dNC', 'PPR', 'aNPB', 'pNPB'), split_by = 'scHelper_cell_type', rerun_UMAP = TRUE)
DimPlot(subset, group.by = 'scHelper_cell_type')
DimPlot(subset, group.by = 'stage')


antler_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/ss8_splitstage_data/antler/stage_gene_modules/rds_files/antler_out.RDS')
ppr_gm <- unlist(antler_data$gene_modules$lists$unbiasedGMs_DE$content[c('GM5')])
nc_gm <- unlist(antler_data$gene_modules$lists$unbiasedGMs_DE$content[c('GM2')])

coexpression(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), show_bins = TRUE, bin_number = 10)
coexpression_highlight_cells(subset, gm_1 = ppr_gm, gm_2 = nc_gm, bin_number = 10)

bin<- extract_bin(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = 10, bin_extract = 3)
# bin2 <- extract_bin(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = 10, bin_extract = 8)
# bin = c(bin1, bin2)
DimPlot(hh7, cells.highlight = bin) + NoLegend()
DimPlot(ss4, cells.highlight = bin) + NoLegend()
DimPlot(ss8, cells.highlight = bin) + NoLegend()

####################################################################################################
# Run coexpression using MB and NC gene modules from ss8
subset <- readRDS('./output/NF-downstream_analysis_stacas/transfer_subset/transfer_npb_subset/seurat/transfer_cluster/rds_files/transfer_clustered_data.RDS')
subset <- subset_seurat(subset, population = c('NC', 'dNC', 'iNP', 'MB'), split_by = 'scHelper_cell_type', rerun_UMAP = TRUE)
DimPlot(subset, group.by = 'scHelper_cell_type')
DimPlot(subset, group.by = 'stage')


antler_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/ss8_splitstage_data/antler/stage_gene_modules/rds_files/antler_out.RDS')
mbhb_gm <- unlist(antler_data$gene_modules$lists$unbiasedGMs_DE$content[c('GM16', 'GM15')])
nc_gm <- unlist(antler_data$gene_modules$lists$unbiasedGMs_DE$content[c('GM2', 'GM1')])

coexpression(subset, gm_1 = mbhb_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), show_bins = TRUE, bin_number = 10)
coexpression_highlight_cells(subset, gm_1 = mbhb_gm, gm_2 = nc_gm, bin_number = 10)

bin1 <- extract_bin(subset, gm_1 = mbhb_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = 10, bin_extract = 5)
bin2 <- extract_bin(subset, gm_1 = mbhb_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = 10, bin_extract = 5)
bin = c(bin1, bin2)
DimPlot(hh7, cells.highlight = bin) + NoLegend()
DimPlot(ss4, cells.highlight = bin) + NoLegend()
DimPlot(ss8, cells.highlight = bin) + NoLegend()


####################################################################################################
# Run coexpression using FB and PPR gene modules from ss8
subset <- readRDS('./output/NF-downstream_analysis_stacas/transfer_subset/transfer_ppr_fb_subset/seurat/transfer_cluster/rds_files/transfer_clustered_data.RDS')
# subset <- subset_seurat(subset, population = c('NC', 'dNC', 'iNP', 'MB'), split_by = 'scHelper_cell_type', rerun_UMAP = TRUE)
subset <- subset_seurat(subset, population = c('hh6', 'hh7', 'ss4', 'ss8'), split_by = 'stage', rerun_UMAP = TRUE)

DimPlot(subset, group.by = 'scHelper_cell_type')
DimPlot(subset, group.by = 'stage')

antler_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/ss8_splitstage_data/antler/stage_gene_modules/rds_files/antler_out.RDS')
aPPR_gm <- unlist(antler_data$gene_modules$lists$unbiasedGMs_DE$content[c('GM5', 'GM6')])
fb_gm <- unlist(antler_data$gene_modules$lists$unbiasedGMs_DE$content[c('GM11', 'GM10')])


coexpression_highlight_cells(subset, gm_1 = aPPR_gm, gm_2 = fb_gm, bin_number = 10)
coexpression(subset, gm_1 = aPPR_gm, gm_2 = fb_gm, meta_data = c('scHelper_cell_type'), show_bins = TRUE, bin_number = 10)

bin1 <- extract_bin(subset, gm_1 = aPPR_gm, gm_2 = fb_gm, meta_data = c('scHelper_cell_type'), bin_number = 10, bin_extract = 1)
bin2 <- extract_bin(subset, gm_1 = aPPR_gm, gm_2 = fb_gm, meta_data = c('scHelper_cell_type'), bin_number = 10, bin_extract = 1)
bin = c(bin1, bin2)
DimPlot(ss4, cells.highlight = bin) + NoLegend()
DimPlot(ss8, cells.highlight = bin) + NoLegend()
DimPlot(hh7, cells.highlight = bin) + NoLegend()
DimPlot(hh6, cells.highlight = bin) + NoLegend()


####################################################################################################
# Run coexpression using FB and HB gene modules from ss8
subset <- readRDS('./output/NF-downstream_analysis_stacas/transfer_subset/transfer_fb_mbhb_subset/seurat/transfer_cluster/rds_files/transfer_clustered_data.RDS')
# subset <- subset_seurat(subset, population = c('NC', 'dNC', 'iNP', 'MB'), split_by = 'scHelper_cell_type', rerun_UMAP = TRUE)
subset <- subset_seurat(subset, population = c('hh6', 'hh7', 'ss4', 'ss8'), split_by = 'stage', rerun_UMAP = TRUE)
subset <- subset_seurat(subset, population = c('aNP', 'iNP', 'MB', 'FB'), split_by = 'scHelper_cell_type', rerun_UMAP = TRUE)

DimPlot(subset, group.by = 'scHelper_cell_type')

antler_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/ss8_splitstage_data/antler/stage_gene_modules/rds_files/antler_out.RDS')
mb_gm <- unlist(antler_data$gene_modules$lists$unbiasedGMs_DE$content[c('GM15')])
fb_gm <- unlist(antler_data$gene_modules$lists$unbiasedGMs_DE$content[c('GM10')])


coexpression_highlight_cells(subset, gm_1 = mb_gm, gm_2 = fb_gm, bin_number = 10)
coexpression(subset, gm_1 = mb_gm, gm_2 = fb_gm, meta_data = c('scHelper_cell_type'), show_bins = TRUE, bin_number = 10)

bin1 <- extract_bin(subset, gm_1 = mb_gm, gm_2 = fb_gm, meta_data = c('scHelper_cell_type'), bin_number = 10, bin_extract = 6)
bin2 <- extract_bin(subset, gm_1 = mb_gm, gm_2 = fb_gm, meta_data = c('scHelper_cell_type'), bin_number = 10, bin_extract = 7)
bin = c(bin1, bin2)
DimPlot(ss4, cells.highlight = bin) + NoLegend()
DimPlot(ss8, cells.highlight = bin) + NoLegend()
DimPlot(hh7, cells.highlight = bin) + NoLegend()
DimPlot(hh6, cells.highlight = bin) + NoLegend()













