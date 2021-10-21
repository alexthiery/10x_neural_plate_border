
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



coexpression = function(seurat_object, gm_1, gm_2, facet_names = NULL, meta_col='scHelper_cell_type', show_bins = FALSE, bin_number = 10, extract_bins = NULL, bin_colour = 'Set2'){
  
  meta_data = seurat_object@meta.data[,meta_col,drop=FALSE]
  gm_1_sum <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data'))[,gm_1] %>% rowSums(.)
  gm_2_sum <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data'))[,gm_2] %>% rowSums(.)
  
  plot_data <- data.frame(gm_1_sum = gm_1_sum, gm_2_sum = gm_2_sum, meta_data)
  
  plot_data <- plot_data %>%
    mutate(ratio = gm_1_sum/(gm_1_sum + gm_2_sum)) %>%
    arrange(ratio) %>%
    pivot_longer(cols = c(gm_1_sum, gm_2_sum)) %>%
    rename(expression_magnitude = value)
  
  # Replace GM names for facet labels
  if(!is.null(facet_names)){
    if(length(facet_names) != 2){
      stop('facet_names must be a character vector of length 2')
    }
    plot_data$name <- c('gm_1_sum' = facet_names[1], 'gm_2_sum' = facet_names[2])[plot_data$name]
  }
  
  plot_data$cell_order <- 1:nrow(plot_data)
  
  p1 = ggplot(plot_data, aes(x = cell_order, y = expression_magnitude)) +
    theme_classic() +
    geom_point(aes(colour = scHelper_cell_type), alpha = 0.5, size = 0.5) +
    scale_colour_manual(values = ggPlotColours(length(unique((plot_data[[meta_col]])))), guide = guide_legend(override.aes = list(size=3),
                                                                                                              title = '')) +
    geom_smooth(method = "gam", se = FALSE, colour="purple", size=1.5) +
    facet_wrap(~name, dir = "v", scales = "free") +
    xlab("Cells") + ylab("Expression level") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size=10),
          strip.text.x = element_text(size = 10))
  
  if (show_bins == TRUE){
    bin_size = length(plot_data$cell_order) / bin_number
    
    start <- seq(min(plot_data$cell_order), max(plot_data$cell_order), by = bin_size)
    end <- seq(min(plot_data$cell_order + bin_size-1), max(plot_data$cell_order + bin_size-1), by = bin_size)
    bins <- data.frame(start, end) %>% mutate(x = ifelse(row_number() %% 2 == 0, 'even', 'odd'))
    
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
    
    p1 <- p1 +
      new_scale_color() +
      geom_rect(data = bins, inherit.aes = FALSE, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = x), show.legend = F) +
      scale_fill_manual(values = alpha(c("gray20", "gray80"), alpha = 0.07)) +
      new_scale_color() +
      geom_rect(data = bin_plot_data, inherit.aes = FALSE, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, colour = as.factor(bin)), fill=NA, show.legend = F) +
      viridis::scale_color_viridis(discrete=TRUE)
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
  start = bin_size*(bin_extract[1]-1)
  end = bin_size*(bin_extract[length(bin_extract)])-1
  # end = start + (bin_size-1)
  return(ordered_cells[start:end])
}


plot_umap_gm_coexpression <- function(seurat_object, gm_1, gm_2){
  start = 0
  end = 100
  width = end - start
  gm_1_means <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'scale.data'))[,gm_1] %>% rowMeans(.)
  gm_1_scaled <- (gm_1_means - min(gm_1_means))/(max(gm_1_means) - min(gm_1_means)) * width + start
  gm_2_means <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'scale.data'))[,gm_2] %>% rowMeans(.)
  gm_2_scaled <- (gm_2_means - min(gm_2_means))/(max(gm_2_means) - min(gm_2_means)) * width + start
  dat <- data.frame(gm_1_scaled, gm_2_scaled, row.names = names(gm_1_scaled))
  dat <-  round(dat, 0)
  col.mat <- expand.grid(a=seq(0,100,by=1), b=seq(0,100,by=1))
  col.mat <- within(col.mat, mix <- rgb(green = a, red = b, blue = b, maxColorValue = 100))
  cell.cols <- unlist(apply(dat, 1, function(x){filter(col.mat, a == x[[1]] & b == x[[2]])[[3]]}))
  col.mat[,1:2] <- col.mat[,1:2]/100
  key.plot <- ggplot(col.mat, aes(x = col.mat[,1], y = col.mat[,2])) +
    xlab("Gene module 1") +
    ylab("Gene module 2") +
    geom_tile(aes(fill = mix)) +
    scale_fill_identity() +
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  plot <- as.data.frame(seurat_object[["umap"]]@cell.embeddings)
  umap.plot <- ggplot(plot, aes(x = plot[,1], y = plot[,2], color = rownames(plot))) +
    geom_point() +
    xlab("UMAP 1")+
    ylab("UMAP 2")+
    scale_color_manual(values=cell.cols)+
    scale_fill_manual(values=cell.cols)+
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  gridExtra::grid.arrange(umap.plot, key.plot, layout_matrix = rbind(c(1,1,2),
                                                                     c(1,1,NA)))
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
install.packages("ggnewscale")
library(ggnewscale)

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


####################################################################################################
# Run coexpression using PPR and NC gene modules from ss8
subset <- readRDS('./output/NF-downstream_analysis_stacas/transfer_subset/transfer_ppr_nc_subset/seurat/transfer_cluster/rds_files/transfer_clustered_data.RDS')

subset <- subset_seurat(subset, population = c('hh7', 'ss4', 'ss8'), split_by = 'stage', rerun_UMAP = FALSE)
subset <- subset_seurat(subset, population = c('pPPR', 'aPPR', 'NC', 'delaminating_NC', 'PPR', 'aNPB', 'pNPB'), split_by = 'scHelper_cell_type', rerun_UMAP = TRUE)
DimPlot(subset, group.by = 'scHelper_cell_type')
DimPlot(subset, group.by = 'stage')


antler_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/ss8_splitstage_data/antler/stage_gene_modules/rds_files/antler_out.RDS')
ppr_gm <- unlist(antler_data$gene_modules$lists$unbiasedGMs_DE$content[c('GM5')])
nc_gm <- unlist(antler_data$gene_modules$lists$unbiasedGMs_DE$content[c('GM2')])


coexpression_highlight_cells(subset, gm_1 = ppr_gm, gm_2 = nc_gm, bin_number = 10)


# Plot module expression together with UMAPS
extract_bins = list(c(1), c(4), c(10))
plots = list()
plots[[1]] <- coexpression(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_col = c('scHelper_cell_type'), show_bins = TRUE, bin_number = 10, extract_bins = extract_bins, facet_names = c('PPR module', 'NC module'))
for(i in 1:length(extract_bins)){
  plots[[i+1]] <- DimPlot(object = subset, cells.highlight = extract_bin(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = 10, bin_extract = extract_bins[[i]]), cols.highlight = viridis::viridis(3)[i]) +
    theme_void() +
    NoLegend()
}

layout <- '
AAAAAAAA
BB#CC#DD
'
png(paste0(plot_path, 'binned_mod_dyn.png'), width = 20, height = 20, res = 200, units = 'cm')
wrap_plots(A = plots[[1]], B = plots[[2]], C = plots[[3]], D = plots[[4]], design = layout, heights = c(2,1))
graphics.off()


# Highlight bins on different stage subsets
bin_sub <- lapply(extract_bins, function(x) extract_bin(subset, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = 10, bin_extract = x))

# hh7
hh7@meta.data$bin_class <- apply(hh7@meta.data %>% rownames_to_column, 1, function(x) if(x["rowname"] %in% bin_sub[[1]]){"a"} else if(x["rowname"] %in% bin_sub[[2]]){"b"} else if(x["rowname"] %in% bin_sub[[3]]){"c"} else {NA})
plot_data <- merge(as.data.frame(hh7[["umap"]]@cell.embeddings), hh7@meta.data[,'bin_class',drop=FALSE], by = 0) %>% column_to_rownames('Row.names')

png(paste0(plot_path, 'highlight_bins_hh7.png'), width = 8, height = 8, res = 200, units = 'cm')
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = 'gray80', size = 1) +
  geom_point(data = plot_data %>% filter(!is.na(bin_class)), inherit.aes = FALSE, aes(x = UMAP_1, y = UMAP_2, colour = factor(bin_class))) +
  viridis::scale_color_viridis(discrete=TRUE) +
  theme_void() +
  NoLegend() +
  ggtitle('hh7') +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 15))
graphics.off()

# 4ss
ss4@meta.data$bin_class <- apply(ss4@meta.data %>% rownames_to_column, 1, function(x) if(x["rowname"] %in% bin_sub[[1]]){"a"} else if(x["rowname"] %in% bin_sub[[2]]){"b"} else if(x["rowname"] %in% bin_sub[[3]]){"c"} else {NA})
plot_data <- merge(as.data.frame(ss4[["umap"]]@cell.embeddings), ss4@meta.data[,'bin_class',drop=FALSE], by = 0) %>% column_to_rownames('Row.names')

png(paste0(plot_path, 'highlight_bins_ss4.png'), width = 8, height = 8, res = 200, units = 'cm')
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = 'gray80', size = 1) +
  geom_point(data = plot_data %>% filter(!is.na(bin_class)), inherit.aes = FALSE, aes(x = UMAP_1, y = UMAP_2, colour = factor(bin_class))) +
  viridis::scale_color_viridis(discrete=TRUE) +
  theme_void() +
  NoLegend() +
  ggtitle('4ss') +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 15))
graphics.off()

# 8ss
ss8@meta.data$bin_class <- apply(ss8@meta.data %>% rownames_to_column, 1, function(x) if(x["rowname"] %in% bin_sub[[1]]){"a"} else if(x["rowname"] %in% bin_sub[[2]]){"b"} else if(x["rowname"] %in% bin_sub[[3]]){"c"} else {NA})
plot_data <- merge(as.data.frame(ss8[["umap"]]@cell.embeddings), ss8@meta.data[,'bin_class',drop=FALSE], by = 0) %>% column_to_rownames('Row.names')

png(paste0(plot_path, 'highlight_bins_ss8.png'), width = 8, height = 8, res = 200, units = 'cm')
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(colour = 'gray80', size = 1) +
  geom_point(data = plot_data %>% filter(!is.na(bin_class)), inherit.aes = FALSE, aes(x = UMAP_1, y = UMAP_2, colour = factor(bin_class))) +
  viridis::scale_color_viridis(discrete=TRUE) +
  theme_void() +
  NoLegend() +
  ggtitle('8ss') +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 15))
graphics.off()



plot_umap_gm_coexpression <- function(seurat_object, gm_1, gm_2, col.threshold = 0.25, two.colors = c('#FF0000', '#00ff00'), negative.color = 'gray80', limit = 0){
  start = 0
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
  col.mat = Seurat:::BlendMatrix(n = 100, col.threshold = 0, two.colors =  two.colors, negative.color = 'gray90')
  col.mat <- as.data.frame.table(col.mat, responseName = "value") %>% mutate_if(is.factor, as.integer)
  col.mat[!(col.mat$Var1 > limit & col.mat$Var2 > limit), 'value'] <- 'gray90'
  colnames(col.mat) <- c('a', 'b', 'mix')
  
  cell.cols <- unlist(apply(dat, 1, function(x){filter(col.mat, a == x[[1]] & b == x[[2]])[[3]]}))
  col.mat[,1:2] <- col.mat[,1:2]/100
  key.plot <- ggplot(col.mat, aes(x = col.mat[,1], y = col.mat[,2])) +
    xlab("Gene module 1") +
    ylab("Gene module 2") +
    geom_tile(aes(fill = mix)) +
    scale_fill_identity() +
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  plot <- as.data.frame(seurat_object[["umap"]]@cell.embeddings)
  umap.plot <- ggplot(plot, aes(x = plot[,1], y = plot[,2], color = rownames(plot))) +
    geom_point() +
    xlab("UMAP 1")+
    ylab("UMAP 2")+
    scale_color_manual(values=cell.cols)+
    scale_fill_manual(values=cell.cols)+
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  gridExtra::grid.arrange(umap.plot, key.plot, layout_matrix = rbind(c(1,1,2),
                                                                     c(1,1,NA)))
}


plot_umap_gm_coexpression(ss4, gm_1 = ppr_gm, gm_2 = nc_gm, col.threshold = 0, two.colors = c("green", 
                                                                                              "red"), negative.color = 'gray90', limit = 25)

