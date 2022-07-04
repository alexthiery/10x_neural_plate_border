#!/usr/bin/env Rscript


############################################################################################
plot_umap_gm_coexpression <- function(seurat_object, gm_1, gm_2, col.threshold = 0, two.colors = c('#FF0000', '#00ff00'), negative.color = 'gray80', limit = 0, 
                                      highlight_cell_size = 1, module_names = c('Gene module 1', 'Gene module 2'), show_legend = TRUE,
                                      axes_label_size = 12, axes_title_size = 10, axes_tick_size = 0.15){
  start = 1
  end = 100
  width = end - start
  gm_1_means <- t(as.matrix(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data')))[,gm_1] %>% rowMeans(.)
  gm_1_scaled <- (gm_1_means - min(gm_1_means))/(max(gm_1_means) - min(gm_1_means)) * width + start
  gm_2_means <- t(as.matrix(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data')))[,gm_2] %>% rowMeans(.)
  gm_2_scaled <- (gm_2_means - min(gm_2_means))/(max(gm_2_means) - min(gm_2_means)) * width + start
  dat <- data.frame(gm_1_scaled, gm_2_scaled, row.names = names(gm_1_scaled))
  dat <-  round(dat, 0)
  # col.mat <- expand.grid(a=seq(0,100,by=1), b=seq(0,100,by=1))
  # col.mat <- within(col.mat, mix <- rgb(green = a, red = a, blue = 0, maxColorValue = 100))
  col_mat = Seurat:::BlendMatrix(n = 100, col.threshold = col.threshold, two.colors =  two.colors, negative.color = negative.color)
  col_mat <- as.data.frame.table(col_mat, responseName = "value") %>% mutate_if(is.factor, as.integer)
  col_mat[!(col_mat$Var1 > limit*100 & col_mat$Var2 > limit*100), 'value'] <- negative.color
  colnames(col_mat) <- c('a', 'b', 'mix')
  
  cell_cols <- unlist(apply(dat, 1, function(x){filter(col_mat, a == x[[1]] & b == x[[2]])[[3]]}))
  col_mat[,1:2] <- col_mat[,1:2]/100
  key_plot <- ggplot(col_mat %>% filter(mix != !!negative.color), aes(x = a, y = b)) +
    xlab(module_names[1]) +
    ylab(module_names[2]) +
    geom_tile(aes(fill = mix)) +
    scale_fill_identity() +
    scale_x_continuous(breaks = c(limit, 1), expand = c(0.01, 0.01)) +
    scale_y_continuous(breaks = c(limit, 1), expand = c(0.01, 0.01)) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_text(size = axes_title_size),
          axis.text.x = element_text(size = axes_label_size),
          axis.title.y = element_text(size = axes_title_size),
          axis.text.y = element_text(size = axes_label_size),
          axis.ticks.length=unit(axes_tick_size,"cm"))
  
  plot_data <- as.data.frame(seurat_object[["umap"]]@cell.embeddings)
  
  # add cell colours to plot_data
  plot_data <- merge(plot_data, as.data.frame(cell_cols), by=0) %>% column_to_rownames('Row.names')
  
  positive_cells = filter(plot_data, cell_cols != negative.color)
  
  umap_plot <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, colour = rownames(positive_cells))) +
    geom_point(colour = negative.color, size = 2) +
    geom_point(data = positive_cells, size = highlight_cell_size) +
    scale_colour_manual(breaks = rownames(positive_cells), values=positive_cells$cell_cols) +
    theme_void() +
    NoLegend()
  
  if (show_legend == FALSE){
    print(umap_plot)
  } else {
    layout <- '
    BA
    B#
    '
    wrap_plots(A = key_plot, B = umap_plot, design = layout, widths = c(4,1), heights = c(1,3))
  }
}
#######################################################################################################


# Define arguments for Rscript
library(getopt)
library(future)
library(Seurat)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library(scHelper)
library(gridExtra)
library(grid)
set.seed(100)
library(mgcv)
library(viridis)
library(ggnewscale)
library(patchwork)

spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
plot_path = "./plots/"
rds_path = "./rds_files/"
data_path = "./input/"
ncores = opt$cores

dir.create(plot_path, recursive = T)
dir.create(rds_path, recursive = T)

antler_data <- readRDS(list.files(data_path, pattern = "antler_out.RDS", full.names = TRUE))
subset <- readRDS(list.files(data_path, pattern = "^transfer_clustered_data", full.names = TRUE))

######################################################
# Calculate co-expression of placodal and NC modules calculated from the NPB subset

ppr_gm <- unlist(antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content[c('GM13')])
for(gm in c('GM40', 'GM42', 'GM44', 'GM43')){
  nc_gm <- unlist(antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content[c(gm)])
  
  plot = plot_umap_gm_coexpression(subset, gm_1 = ppr_gm, gm_2 = nc_gm, col.threshold = 0, two.colors = c("red", "blue"),
                                   negative.color = 'gray90', limit = 0.3, module_names = c('PPR module', 'NC module'), highlight_cell_size = 2,
                                   show_legend = TRUE, axes_label_size = 16, axes_title_size = 14, axes_tick_size = 0.25)
  
  png(paste0(plot_path, 'NPB_UMAP_coexpression_', gm, '.png'), width = 20, height = 14, res = 600, units = 'cm')
  print(plot)
  graphics.off()
  
  # Plot co-expression no legend
  plot = plot_umap_gm_coexpression(subset, gm_1 = ppr_gm, gm_2 = nc_gm, col.threshold = 0, two.colors = c("red", "blue"),
                                   negative.color = 'gray90', limit = 0.3, module_names = c('PPR module', 'NC module'), highlight_cell_size = 2,
                                   show_legend = FALSE, axes_label_size = 16, axes_title_size = 14, axes_tick_size = 0.25)
  
  png(paste0(plot_path, 'NPB_UMAP_coexpression_no_legend_', gm, '.png'), width = 14, height = 14, res = 600, units = 'cm')
  print(plot)
  graphics.off()
}

