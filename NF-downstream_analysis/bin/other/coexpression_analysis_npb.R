#!/usr/bin/env Rscript

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

########################       CELL STATE COLOURS    ########################################
scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR',
                              'eNPB', 'NPB', 'aNPB', 'pNPB','NC', 'dNC',
                              'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                              'aNP', 'FB', 'vFB', 'node', 'streak')

scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", "#53A651", "#6D8470",
                                "#87638F", "#A5548D", "#C96555", "#ED761C", "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30",
                                "#CC9F2C", "#AD6428", "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3")

names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 'MB', 
                                       'vFB', 'aNP', 'node', 'FB', 'pEpi')
########################       STAGE COLOURS     ###########################################
stage_order <- c("HH4", "HH5", "HH6", "HH7", "ss4", "ss8")

stage_colours = c("#E78AC3", "#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order
############################################################################################

#######################################################################################################
#########################################   FUNCTIONS   ###############################################
named_list_to_vector <- function(list){
  if(is.null(names(list))){stop('list must be named')}
  unlist(lapply(lapply(names(list), function(x) rep(x, length(list[[x]]))), function(y) names(y) <- x))
}

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
  gm_1_sum <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data'))[,gm_1] %>% rowMeans(.)
  gm_2_sum <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data'))[,gm_2] %>% rowMeans(.)
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

coexpression = function(seurat_object, gms, order_1 = 1, order_2 = length(gms), meta_col='scHelper_cell_type', meta_col_colours = NULL,
                        show_bins = FALSE, bin_number = 10, extract_bins = NULL, bin_colour = 'Set2', save_bins = NULL){
  
  meta_data = seurat_object@meta.data[,meta_col,drop=FALSE]
  
  gm_means = lapply(gms, function(x) t(as.matrix(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data')))[,x] %>% rowMeans(.))
  
  if(is.null(names(gm_means))){names(gm_means) <- paste0('GM:', 1:length(gm_means), ' average')}
  
  plot_data <- cbind(do.call(cbind, gm_means), meta_data)
  
  plot_data <- plot_data %>%
    rownames_to_column(var = 'cell_name') %>%
    dplyr::select(-cell_name, cell_name) %>%
    mutate(order_1 = rowMeans(.[, order_1, drop=FALSE])) %>%
    mutate(order_2 = rowMeans(.[, order_2, drop=FALSE])) %>%
    mutate(ratio = order_1/(order_1+order_2)) %>%
    arrange(ratio) %>%
    mutate(cell_order = row_number())
  
  if(show_bins == TRUE | !is.null(save_bins)){
    bin_size = nrow(meta_data) / bin_number
    
    start <- floor(seq(min(plot_data$cell_order), max(plot_data$cell_order), by = bin_size))
    end <- floor(seq(min(plot_data$cell_order + bin_size-1), max(plot_data$cell_order + bin_size-1), by = bin_size))
    bins <- data.frame(start, end) %>% mutate(x = ifelse(row_number() %% 2 == 0, 'even', 'odd'))
    
    if(!is.null(save_bins)){
      plots = lapply(1:nrow(bins), function(x) DimPlot(seurat_object, cells.highlight = plot_data$cell_name[c(bins[x, 'start'] : bins[x, 'end'])], combine = TRUE) + NoLegend() + ggtitle(paste("Bin ", x)))
      dir.create(dirname(save_bins), recursive = TRUE)
      png(save_bins, width = floor(sqrt(bin_number)) * 7, height = ceiling(sqrt(bin_number)) * 7, units = 'cm', res = 200)
      grid.arrange(grobs = plots)
      graphics.off()
    }
    
    if(!is.null(extract_bins)){
      if(class(extract_bins) == 'list'){
        bin_plot_data <- data.frame()
        for(bin in extract_bins){
          sub_bin = bins[bin,]
          bin_plot_data <- rbind(bin_plot_data, c(sub_bin[1,][,1], sub_bin[nrow(sub_bin),][,2]))
        }
        bin_plot_data[['bin']] <- 1:nrow(bin_plot_data)
        colnames(bin_plot_data) <- c('start', 'end', 'bin')
      } else if(class(extract_bins) == 'numeric')
        bin_plot_data <- bins[extract_bins,1:2]
      bin_plot_data[['bin']] <- 1:nrow(bin_plot_data)
    } else {
      stop("'extract_bins' must be a numeric array or a list containing numeric arrays")
    }
    
    plot_data = plot_data %>%
      dplyr::select(-c(order_1, order_2, ratio)) %>%
      pivot_longer(cols = -c(meta_col, cell_order, cell_name), names_to = 'gms', values_to = 'expression_magnitude')
    
    if (is.null(meta_col_colours)){
      cols = ggPlotColours(length(unique((plot_data[[meta_col]]))))
    } else {cols = meta_col_colours}
    
    p1 = ggplot(plot_data, aes(x = cell_order, y = expression_magnitude)) +
      theme_classic() +
      geom_point(aes(colour = as.factor(scHelper_cell_type)), alpha = 0.5, size = 0.5) +
      scale_colour_manual(values = cols, guide = guide_legend(override.aes = list(size=6), title = '')) +
      geom_smooth(method = "gam", se = FALSE, colour="purple", size=1.5) +
      facet_wrap(~gms, dir = "v", scales = "free", ncol=1) +
      xlab("Cells") + ylab("Expression level") +
      theme(axis.text.x=element_blank(),
            axis.ticks.x = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(size=10),
            strip.text.x = element_text(size = 10))
    
    
    p1 <- p1 +
      new_scale_color() +
      geom_rect(data = bins, inherit.aes = FALSE, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = x), show.legend = F) +
      scale_fill_manual(values = alpha(c("gray20", "gray80"), alpha = 0.07)) +
      new_scale_color() +
      geom_rect(data = bin_plot_data, inherit.aes = FALSE, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, colour = as.factor(bin)), fill=NA, show.legend = F) +
      scale_colour_manual(values = c("#FDE725FF", "#440154FF", "#21908CFF"))
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
  gm_1_sum <- t(as.matrix(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data')))[,gm_1] %>% rowSums(.)
  gm_2_sum <- t(as.matrix(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data')))[,gm_2] %>% rowSums(.)
  plot_data <- data.frame(gm_1_sum = gm_1_sum, gm_2_sum = gm_2_sum, x = x)
  plot_data <- plot_data %>% mutate(ratio = gm_1_sum/(gm_1_sum + gm_2_sum)) %>% arrange(ratio)
  ordered_cells <- rownames(plot_data)
  
  bin_size <- length(ordered_cells)/bin_number
  start = bin_size*(bin_extract[1]-1)
  end = bin_size*(bin_extract[length(bin_extract)])-1
  # end = start + (bin_size-1)
  return(ordered_cells[start:end])
}

plot_umap_gm_coexpression <- function(seurat_object, gm_1, gm_2, col.threshold = 0.25, two.colors = c('#FF0000', '#00ff00'), negative.color = 'gray80', limit = 0, 
                                      highlight_cell_size = 1, module_names = c('Gene module 1', 'Gene module 2'), show_legend = TRUE){
  start = 1
  end = 100
  width = end - start
  gm_1_means <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'scale.data'))[,gm_1] %>% rowMeans(.)
  gm_1_scaled <- (gm_1_means - min(gm_1_means))/(max(gm_1_means) - min(gm_1_means)) * width + start
  gm_2_means <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'scale.data'))[,gm_2] %>% rowMeans(.)
  gm_2_scaled <- (gm_2_means - min(gm_2_means))/(max(gm_2_means) - min(gm_2_means)) * width + start
  dat <- data.frame(gm_1_scaled, gm_2_scaled, row.names = names(gm_1_scaled))
  dat <-  round(dat, 0)
  # col.mat <- expand.grid(a=seq(0,100,by=1), b=seq(0,100,by=1))
  # col.mat <- within(col.mat, mix <- rgb(green = a, red = a, blue = 0, maxColorValue = 100))
  col_mat = Seurat:::BlendMatrix(n = 100, col.threshold = 0, two.colors =  two.colors, negative.color = negative.color)
  col_mat <- as.data.frame.table(col_mat, responseName = "value") %>% mutate_if(is.factor, as.integer)
  col_mat[!(col_mat$Var1 > limit & col_mat$Var2 > limit), 'value'] <- negative.color
  colnames(col_mat) <- c('a', 'b', 'mix')
  
  cell_cols <- unlist(apply(dat, 1, function(x){filter(col_mat, a == x[[1]] & b == x[[2]])[[3]]}))
  col_mat[,1:2] <- col_mat[,1:2]/100
  key_plot <- ggplot(col_mat, aes(x = col_mat[,1], y = col_mat[,2])) +
    xlab(module_names[1]) +
    ylab(module_names[2]) +
    geom_tile(aes(fill = mix)) +
    scale_fill_identity() +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  plot_data <- as.data.frame(seurat_object[["umap"]]@cell.embeddings)
  
  # add cell colours to plot_data
  plot_data <- merge(plot_data, as.data.frame(cell_cols), by=0) %>% column_to_rownames('Row.names')
  
  umap_plot <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, colour = cell_cols)) +
    geom_point(colour = negative.color, size = 2) +
    geom_point(data = plot_data %>% filter(cell_cols != negative.color), size = highlight_cell_size) +
    scale_colour_manual(values=plot_data %>% filter(cell_cols != negative.color) %>% dplyr::pull(cell_cols))+
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
HH5 <- readRDS(list.files(data_path, pattern = "^HH5", full.names = TRUE))
HH6 <- readRDS(list.files(data_path, pattern = "^HH6", full.names = TRUE))
HH7 <- readRDS(list.files(data_path, pattern = "^HH7", full.names = TRUE))
ss4 <- readRDS(list.files(data_path, pattern = "^ss4", full.names = TRUE))
ss8 <- readRDS(list.files(data_path, pattern = "^ss8", full.names = TRUE))
subset <- readRDS(list.files(data_path, pattern = "^transfer_clustered_data", full.names = TRUE))

# seurat_data <- readRDS(list.files(data_path, pattern = "*.RDS", full.names = TRUE)[!list.files(data_path, pattern = "*.RDS") %>% grepl('antler', .)])
# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/transfer_labels/seurat/rds_files/seurat_label_transfer.RDS')
# subset <- readRDS('./output/NF-downstream_analysis_stacas/transfer_subset/transfer_ppr_nc_subset/seurat/transfer_cluster/rds_files/transfer_clustered_data.RDS')
# HH5 <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/HH5_splitstage_data/seurat/stage_state_classification/rds_files/HH5_cell_state_classification.RDS')
# HH6 <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/HH6_splitstage_data/seurat/stage_state_classification/rds_files/HH6_cell_state_classification.RDS')
# HH7 <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/HH7_splitstage_data/seurat/stage_state_classification/rds_files/HH7_cell_state_classification.RDS')
# ss4 <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/ss4_splitstage_data/seurat/stage_state_classification/rds_files/ss4_cell_state_classification.RDS')
# ss8 <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/ss8_splitstage_data/seurat/stage_state_classification/rds_files/ss8_cell_state_classification.RDS')
# antler_data <- readRDS("~/output/NF-downstream_analysis_stacas/stage_split/ss8_splitstage_data/antler/stage_gene_modules/rds_files/antler_out.RDS")

####################################################################################################
# Run coexpression using PPR and NC gene modules from ss8
ppr_gm <- unlist(antler_data$gene_modules$lists$unbiasedGMs_DE$content[c('GM5')])
nc_gm <- unlist(antler_data$gene_modules$lists$unbiasedGMs_DE$content[c('GM2')])
bin_number = 10

# Plot module coexpression dynamics
scHelper_cols_subset <- scHelper_cell_type_colours[levels(droplevels(subset@meta.data$scHelper_cell_type))]
extract_bins = c(1, 3, 10)

coexpression_plot = coexpression(subset, gms = list('PPR GM' = ppr_gm, 'NC GM' = nc_gm), meta_col = c('scHelper_cell_type'), meta_col_colours = scHelper_cols_subset,
                    order_1 = 1, order_2 = c(2), show_bins = TRUE, bin_number = bin_number, extract_bins = extract_bins,
                    save_bins = paste0(plot_path, 'bins.png'))

png(paste0(plot_path, 'coexpression_dynamics.png'), width = 20, height = 20, res = 200, units = 'cm')
coexpression_plot
graphics.off()

# Plot bins on UMAP
for(i in 1:length(extract_bins)){
  png(paste0(plot_path, 'bin_', extract_bins[[i]], '.png'), width = 10, height = 10, res = 200, units = 'cm')
  print(DimPlot(object = subset, cells.highlight = extract_bin(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = bin_number, bin_extract = extract_bins[[i]]), 
                cols.highlight = c("#FDE725FF", "#440154FF", "#21908CFF")[i]) +
          theme_void() +
          NoLegend() +
          ggtitle(paste0('Bin ', extract_bins[[i]])) +
          theme(plot.title = element_text(hjust = 0.5, colour = c("#FDE725FF", "#440154FF", "#21908CFF")[i], face = 'bold', size = 20)))
  graphics.off()
}


################ DimPlot of scHelper_cell_types and stage in npb subset
scHelper_cols_subset <- scHelper_cell_type_colours[levels(droplevels(subset@meta.data$scHelper_cell_type))]
stage_cols <- stage_colours[levels(droplevels(subset@meta.data$stage))]

png(paste0(plot_path, "npb_subset_scHelper_celltype_umap.png"), width=29, height=20, units = 'cm', res = 200)
DimPlot(subset, group.by = 'scHelper_cell_type', label = FALSE, 
        label.size = 10,
        label.box = TRUE, repel = TRUE,
        pt.size = 1.2, 
        cols = scHelper_cols_subset, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(plot.title = element_blank(),
                 legend.key.size = unit(10, 'cm'), #change legend key size
                 legend.key.height = unit(1, 'cm'), #change legend key height
                 legend.key.width = unit(1, 'cm'), #change legend key width
                 legend.text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size = 8)))
graphics.off()

png(paste0(plot_path, "npb_subset_scHelper_celltype_umap_no_legend.png"), width=26, height=20, units = 'cm', res = 200)
DimPlot(subset, group.by = 'scHelper_cell_type', label = FALSE, 
        label.size = 10,
        label.box = TRUE, repel = TRUE,
        pt.size = 1.2, 
        cols = scHelper_cols_subset, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(plot.title = element_blank(),
                 legend.position = "NULL")
graphics.off()

png(paste0(plot_path, "npb_subset_stage_umap.png"), width=26, height=20, units = 'cm', res = 200)
DimPlot(subset, group.by = 'stage', label = TRUE, 
        label.size = 15,
        label.box = TRUE, repel = TRUE,
        pt.size = 1.2, 
        cols = stage_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

#############################################################################
############################# Stage subsets #################################
bin_extract = 3
bin_sub <- lapply(extract_bins, function(x) extract_bin(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = bin_number, bin_extract = bin_extract))


# HH5
############################################
plot_data = merge(as.data.frame(HH5[["umap"]]@cell.embeddings), HH5@meta.data, by = 0) %>% column_to_rownames('Row.names')
plot_data <- plot_data %>% mutate(scHelper_cols = ifelse(scHelper_cell_type %in% c('dNC', 'NC', 'pNPB', 'aNPB', 'aPPR', 'PPR', 'pPPR'), scHelper_cell_type_colours[as.vector(scHelper_cell_type)], "gray90"))

# subset all one colour on full UMAP
png(paste0(plot_path, "HH5_subset_full_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = 'gray90', size = 2) +
  geom_point(data = subset(plot_data, scHelper_cell_type %in% c('dNC', 'NC', 'pNPB', 'aNPB', 'aPPR', 'PPR', 'pPPR')), size = 2, aes(colour = "red")) +
  theme_void() +
  theme(legend.position = "none")
graphics.off()

# subset coloured by scHelper_cell_state on full UMAP
png(paste0(plot_path, "HH5_subset_cell_state_full_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = plot_data$scHelper_cols, size = 2) +
  scale_fill_manual(values=as.character(plot_data$scHelper_cols)) +
  theme_void() +
  theme(legend.position = "none")
graphics.off()

# subset coloured by scHelper_cell_state on subset UMAP
plot_data = subset(plot_data, scHelper_cell_type %in% c('dNC', 'NC', 'pNPB', 'aNPB', 'aPPR', 'PPR', 'pPPR'))

png(paste0(plot_path, "HH5_subset_cell_state_subset_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = plot_data$scHelper_cols, size = 2) +
  scale_fill_manual(values=as.character(plot_data$scHelper_cols)) +
  theme_void() +
  theme(legend.position = "none")
graphics.off()

# Highlight bin 3
plot_data$bin_class <- apply(plot_data %>% rownames_to_column, 1, function(x) if(x["rowname"] %in% bin_sub[[1]]){"a"} else if(x["rowname"] %in% bin_sub[[2]]){"b"} else if(x["rowname"] %in% bin_sub[[3]]){"c"} else {NA})

png(paste0(plot_path, 'HH5_highlight_bin.png'), width = 12, height = 12, res = 200, units = 'cm')
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = 'gray90', size = 2) +
  geom_point(data = plot_data %>% filter(!is.na(bin_class)), inherit.aes = FALSE, aes(x = UMAP_1, y = UMAP_2), colour = "#440154FF") +
  theme_void() +
  NoLegend()
graphics.off()

# Plot co-expression with legend
coexpression <- subset(HH5, cells = rownames(plot_data))

png(paste0(plot_path, 'HH5_UMAP_coexpression.png'), width = 14, height = 10, res = 400, units = 'cm')
plot_umap_gm_coexpression(coexpression, gm_1 = ppr_gm, gm_2 = nc_gm, col.threshold = 0, two.colors = c("red", "blue"),
                          negative.color = 'gray90', limit = 30, module_names = c('PPR module', 'NC module'), highlight_cell_size = 2,
                          show_legend = TRUE)
graphics.off()

# Plot co-expression no legend
png(paste0(plot_path, 'HH5_UMAP_coexpression_no_legend.png'), width = 12, height = 12, res = 200, units = 'cm')
plot_umap_gm_coexpression(coexpression, gm_1 = ppr_gm, gm_2 = nc_gm, col.threshold = 0, two.colors = c("red", "blue"),
                          negative.color = 'gray90', limit = 30, module_names = c('PPR module', 'NC module'), highlight_cell_size = 2,
                          show_legend = FALSE)
graphics.off()



# HH6
############################################
plot_data = merge(as.data.frame(HH6[["umap"]]@cell.embeddings), HH6@meta.data, by = 0) %>% column_to_rownames('Row.names')
plot_data <- plot_data %>% mutate(scHelper_cols = ifelse(scHelper_cell_type %in% c('dNC', 'NC', 'pNPB', 'aNPB', 'aPPR', 'PPR', 'pPPR'), scHelper_cell_type_colours[as.vector(scHelper_cell_type)], "gray90"))

# subset all one colour on full UMAP
png(paste0(plot_path, "HH6_subset_full_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = 'gray90', size = 2) +
  geom_point(data = subset(plot_data, scHelper_cell_type %in% c('dNC', 'NC', 'pNPB', 'aNPB', 'aPPR', 'PPR', 'pPPR')), size = 2, aes(colour = "red")) +
  theme_void() +
  theme(legend.position = "none")
graphics.off()

# subset coloured by scHelper_cell_state on full UMAP
png(paste0(plot_path, "HH6_subset_cell_state_full_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = plot_data$scHelper_cols, size = 2) +
  scale_fill_manual(values=as.character(plot_data$scHelper_cols)) +
  theme_void() +
  theme(legend.position = "none")
graphics.off()

# subset coloured by scHelper_cell_state on subset UMAP
plot_data = subset(plot_data, scHelper_cell_type %in% c('dNC', 'NC', 'pNPB', 'aNPB', 'aPPR', 'PPR', 'pPPR'))

png(paste0(plot_path, "HH6_subset_cell_state_subset_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = plot_data$scHelper_cols, size = 2) +
  scale_fill_manual(values=as.character(plot_data$scHelper_cols)) +
  theme_void() +
  theme(legend.position = "none")
graphics.off()

# Highlight bin 3
plot_data$bin_class <- apply(plot_data %>% rownames_to_column, 1, function(x) if(x["rowname"] %in% bin_sub[[1]]){"a"} else if(x["rowname"] %in% bin_sub[[2]]){"b"} else if(x["rowname"] %in% bin_sub[[3]]){"c"} else {NA})

png(paste0(plot_path, 'HH6_highlight_bin.png'), width = 12, height = 12, res = 200, units = 'cm')
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = 'gray90', size = 2) +
  geom_point(data = plot_data %>% filter(!is.na(bin_class)), inherit.aes = FALSE, aes(x = UMAP_1, y = UMAP_2), colour = "#440154FF") +
  theme_void() +
  NoLegend()
graphics.off()

# Plot co-expression with legend
coexpression <- subset(HH6, cells = rownames(plot_data))

png(paste0(plot_path, 'HH6_UMAP_coexpression.png'), width = 14, height = 10, res = 400, units = 'cm')
plot_umap_gm_coexpression(coexpression, gm_1 = ppr_gm, gm_2 = nc_gm, col.threshold = 0, two.colors = c("red", "blue"),
                          negative.color = 'gray90', limit = 30, module_names = c('PPR module', 'NC module'), highlight_cell_size = 2,
                          show_legend = TRUE)
graphics.off()

# Plot co-expression no legend
png(paste0(plot_path, 'HH6_UMAP_coexpression_no_legend.png'), width = 12, height = 12, res = 200, units = 'cm')
plot_umap_gm_coexpression(coexpression, gm_1 = ppr_gm, gm_2 = nc_gm, col.threshold = 0, two.colors = c("red", "blue"),
                          negative.color = 'gray90', limit = 30, module_names = c('PPR module', 'NC module'), highlight_cell_size = 2,
                          show_legend = FALSE)
graphics.off()



# HH7
############################################
plot_data = merge(as.data.frame(HH7[["umap"]]@cell.embeddings), HH7@meta.data, by = 0) %>% column_to_rownames('Row.names')
plot_data <- plot_data %>% mutate(scHelper_cols = ifelse(scHelper_cell_type %in% c('dNC', 'NC', 'pNPB', 'aNPB', 'aPPR', 'PPR', 'pPPR'), scHelper_cell_type_colours[as.vector(scHelper_cell_type)], "gray90"))

# subset all one colour on full UMAP
png(paste0(plot_path, "HH7_subset_full_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = 'gray90', size = 2) +
  geom_point(data = subset(plot_data, scHelper_cell_type %in% c('dNC', 'NC', 'pNPB', 'aNPB', 'aPPR', 'PPR', 'pPPR')), size = 2, aes(colour = "red")) +
  theme_void() +
  theme(legend.position = "none")
graphics.off()

# subset coloured by scHelper_cell_state on full UMAP
png(paste0(plot_path, "HH7_subset_cell_state_full_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = plot_data$scHelper_cols, size = 2) +
  scale_fill_manual(values=as.character(plot_data$scHelper_cols)) +
  theme_void() +
  theme(legend.position = "none")
graphics.off()

# subset coloured by scHelper_cell_state on subset UMAP
plot_data = subset(plot_data, scHelper_cell_type %in% c('dNC', 'NC', 'pNPB', 'aNPB', 'aPPR', 'PPR', 'pPPR'))

png(paste0(plot_path, "HH7_subset_cell_state_subset_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = plot_data$scHelper_cols, size = 2) +
  scale_fill_manual(values=as.character(plot_data$scHelper_cols)) +
  theme_void() +
  theme(legend.position = "none")
graphics.off()


# Highlight bin 3
plot_data$bin_class <- apply(plot_data %>% rownames_to_column, 1, function(x) if(x["rowname"] %in% bin_sub[[1]]){"a"} else if(x["rowname"] %in% bin_sub[[2]]){"b"} else if(x["rowname"] %in% bin_sub[[3]]){"c"} else {NA})

png(paste0(plot_path, 'HH7_highlight_bin.png'), width = 12, height = 12, res = 200, units = 'cm')
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = 'gray90', size = 2) +
  geom_point(data = plot_data %>% filter(!is.na(bin_class)), inherit.aes = FALSE, aes(x = UMAP_1, y = UMAP_2), colour = "#440154FF") +
  theme_void() +
  NoLegend()
graphics.off()

# Plot co-expression with legend
coexpression <- subset(HH7, cells = rownames(plot_data))

png(paste0(plot_path, 'HH7_UMAP_coexpression.png'), width = 14, height = 10, res = 400, units = 'cm')
plot_umap_gm_coexpression(coexpression, gm_1 = ppr_gm, gm_2 = nc_gm, col.threshold = 0, two.colors = c("red", "blue"),
                          negative.color = 'gray90', limit = 30, module_names = c('PPR module', 'NC module'), highlight_cell_size = 2,
                          show_legend = TRUE)
graphics.off()

# Plot co-expression no legend
png(paste0(plot_path, 'HH7_UMAP_coexpression_no_legend.png'), width = 12, height = 12, res = 200, units = 'cm')
plot_umap_gm_coexpression(coexpression, gm_1 = ppr_gm, gm_2 = nc_gm, col.threshold = 0, two.colors = c("red", "blue"),
                          negative.color = 'gray90', limit = 30, module_names = c('PPR module', 'NC module'), highlight_cell_size = 2,
                          show_legend = FALSE)
graphics.off()


# SS4
############################################
plot_data = merge(as.data.frame(ss4[["umap"]]@cell.embeddings), ss4@meta.data, by = 0) %>% column_to_rownames('Row.names')
plot_data <- plot_data %>% mutate(scHelper_cols = ifelse(scHelper_cell_type %in% c('dNC', 'NC', 'pNPB', 'aNPB', 'aPPR', 'PPR', 'pPPR'), scHelper_cell_type_colours[as.vector(scHelper_cell_type)], "gray90"))

# subset all one colour on full UMAP
png(paste0(plot_path, "ss4_subset_full_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = 'gray90', size = 2) +
  geom_point(data = subset(plot_data, scHelper_cell_type %in% c('dNC', 'NC', 'pNPB', 'aNPB', 'aPPR', 'PPR', 'pPPR')), size = 2, aes(colour = "red")) +
  theme_void() +
  theme(legend.position = "none")
graphics.off()

# subset coloured by scHelper_cell_state on full UMAP
png(paste0(plot_path, "ss4_subset_cell_state_full_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = plot_data$scHelper_cols, size = 2) +
  scale_fill_manual(values=as.character(plot_data$scHelper_cols)) +
  theme_void() +
  theme(legend.position = "none")
graphics.off()

# subset coloured by scHelper_cell_state on subset UMAP
plot_data = subset(plot_data, scHelper_cell_type %in% c('dNC', 'NC', 'pNPB', 'aNPB', 'aPPR', 'PPR', 'pPPR'))

png(paste0(plot_path, "ss4_subset_cell_state_subset_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = plot_data$scHelper_cols, size = 2) +
  scale_fill_manual(values=as.character(plot_data$scHelper_cols)) +
  theme_void() +
  theme(legend.position = "none")
graphics.off()


# Highlight bin 3
plot_data$bin_class <- apply(plot_data %>% rownames_to_column, 1, function(x) if(x["rowname"] %in% bin_sub[[1]]){"a"} else if(x["rowname"] %in% bin_sub[[2]]){"b"} else if(x["rowname"] %in% bin_sub[[3]]){"c"} else {NA})

png(paste0(plot_path, 'ss4_highlight_bin.png'), width = 12, height = 12, res = 200, units = 'cm')
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = 'gray90', size = 2) +
  geom_point(data = plot_data %>% filter(!is.na(bin_class)), inherit.aes = FALSE, aes(x = UMAP_1, y = UMAP_2), colour = "#440154FF") +
  theme_void() +
  NoLegend()
graphics.off()

# Plot co-expression with legend
coexpression <- subset(ss4, cells = rownames(plot_data))

png(paste0(plot_path, 'ss4_UMAP_coexpression.png'), width = 14, height = 10, res = 400, units = 'cm')
plot_umap_gm_coexpression(coexpression, gm_1 = ppr_gm, gm_2 = nc_gm, col.threshold = 0, two.colors = c("red", "blue"),
                          negative.color = 'gray90', limit = 30, module_names = c('PPR module', 'NC module'), highlight_cell_size = 2,
                          show_legend = TRUE)
graphics.off()

# Plot co-expression no legend
png(paste0(plot_path, 'ss4_UMAP_coexpression_no_legend.png'), width = 12, height = 12, res = 200, units = 'cm')
plot_umap_gm_coexpression(coexpression, gm_1 = ppr_gm, gm_2 = nc_gm, col.threshold = 0, two.colors = c("red", "blue"),
                          negative.color = 'gray90', limit = 30, module_names = c('PPR module', 'NC module'), highlight_cell_size = 2,
                          show_legend = FALSE)
graphics.off()


# SS8
############################################
plot_data = merge(as.data.frame(ss8[["umap"]]@cell.embeddings), ss8@meta.data, by = 0) %>% column_to_rownames('Row.names')
plot_data <- plot_data %>% mutate(scHelper_cols = ifelse(scHelper_cell_type %in% c('dNC', 'NC', 'pNPB', 'aNPB', 'aPPR', 'PPR', 'pPPR'), scHelper_cell_type_colours[as.vector(scHelper_cell_type)], "gray90"))

# subset all one colour on full UMAP
png(paste0(plot_path, "ss8_subset_full_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = 'gray90', size = 2) +
  geom_point(data = subset(plot_data, scHelper_cell_type %in% c('dNC', 'NC', 'pNPB', 'aNPB', 'aPPR', 'PPR', 'pPPR')), size = 2, aes(colour = "red")) +
  theme_void() +
  theme(legend.position = "none")
graphics.off()

# subset coloured by scHelper_cell_state on full UMAP
png(paste0(plot_path, "ss8_subset_cell_state_full_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = plot_data$scHelper_cols, size = 2) +
  scale_fill_manual(values=as.character(plot_data$scHelper_cols)) +
  theme_void() +
  theme(legend.position = "none")
graphics.off()

# subset coloured by scHelper_cell_state on subset UMAP
plot_data = subset(plot_data, scHelper_cell_type %in% c('dNC', 'NC', 'pNPB', 'aNPB', 'aPPR', 'PPR', 'pPPR'))

png(paste0(plot_path, "ss8_subset_cell_state_subset_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = plot_data$scHelper_cols, size = 2) +
  scale_fill_manual(values=as.character(plot_data$scHelper_cols)) +
  theme_void() +
  theme(legend.position = "none")
graphics.off()


# Highlight bin 3
plot_data$bin_class <- apply(plot_data %>% rownames_to_column, 1, function(x) if(x["rowname"] %in% bin_sub[[1]]){"a"} else if(x["rowname"] %in% bin_sub[[2]]){"b"} else if(x["rowname"] %in% bin_sub[[3]]){"c"} else {NA})

png(paste0(plot_path, 'ss8_highlight_bin.png'), width = 12, height = 12, res = 200, units = 'cm')
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = 'gray90', size = 2) +
  geom_point(data = plot_data %>% filter(!is.na(bin_class)), inherit.aes = FALSE, aes(x = UMAP_1, y = UMAP_2), colour = "#440154FF") +
  theme_void() +
  NoLegend()
graphics.off()

# Plot co-expression with legend
coexpression <- subset(ss8, cells = rownames(plot_data))

png(paste0(plot_path, 'ss8_UMAP_coexpression.png'), width = 14, height = 10, res = 400, units = 'cm')
plot_umap_gm_coexpression(coexpression, gm_1 = ppr_gm, gm_2 = nc_gm, col.threshold = 0, two.colors = c("red", "blue"),
                          negative.color = 'gray90', limit = 30, module_names = c('PPR module', 'NC module'), highlight_cell_size = 2,
                          show_legend = TRUE)
graphics.off()

# Plot co-expression no legend
png(paste0(plot_path, 'ss8_UMAP_coexpression_no_legend.png'), width = 12, height = 12, res = 200, units = 'cm')
plot_umap_gm_coexpression(coexpression, gm_1 = ppr_gm, gm_2 = nc_gm, col.threshold = 0, two.colors = c("red", "blue"),
                          negative.color = 'gray90', limit = 30, module_names = c('PPR module', 'NC module'), highlight_cell_size = 2,
                          show_legend = FALSE)
graphics.off()

