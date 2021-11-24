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

########################       CELL STATE COLOURS    ########################################
scHelper_all_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR',
                        'eNPB', 'NPB', 'aNPB', 'pNPB','NC', 'dNC',
                        'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                        'aNP', 'FB', 'vFB', 'node', 'streak')

scHelper_all_colours <- c("#676060", "#AD2828", "#551616", "#FF0000", "#DE4D00", "#FF8300",
                          "#C8E81E", "#A5E702", "#6EE702", "#16973F", "#19B4A1", "#10E0E8",
                          "#BA3CA5", "#8A4FC5", "#0A0075", "#3B0075", "#8000FF", "#D800FF",
                          "#FF00D4", "#F16DDB", "#FFBAF3", "#B672AA", "#BBBEBE", "#787878")
names(scHelper_all_colours) <- scHelper_all_order
########################       STAGE COLOURS     ###########################################
stage_all_order <- c("hh4", "hh5", "hh6", "hh7", "ss4", "ss8")

stage_all_colours = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F")
names(stage_all_colours) <- stage_all_order
############################################################################################

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    plot_path = "./test/state_classification/plots/"
    rds_path = "./test/state_classification/rds_files/"
    data_path = "./output/NF-downstream_analysis_stacas/seurat_filtering/6_contamination_filt/rds_files/"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
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

# Retrieve seurat object label
label <- sub('_.*', '', list.files(data_path, pattern = '*.RDS'))

seurat_data <- readRDS(list.files(data_path, full.names = TRUE, pattern = '*.RDS'))
# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/HH4_splitstage_data/seurat/stage_cluster/rds_files/HH4_clustered_data.RDS')
# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/HH5_splitstage_data/seurat/stage_cluster/rds_files/HH5_clustered_data.RDS')
# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/HH6_splitstage_data/seurat/stage_cluster/rds_files/HH6_clustered_data.RDS')
# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/HH7_splitstage_data/seurat/stage_cluster/rds_files/HH7_clustered_data.RDS')
# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/ss4_splitstage_data/seurat/stage_cluster/rds_files/ss4_clustered_data.RDS')
# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/ss8_splitstage_data/seurat/stage_cluster/rds_files/ss8_clustered_data.RDS')

########################################################################################################
#                                      Cell state classification                                    #
#######################################################################################
# Convert knowledge matrix to gene list
# cell_state_markers <- read.csv('./NF-downstream_analysis/binary_knowledge_matrix.csv', row.names = 1) %>% select(!c(evidence, PPR, NP, NPB, iNP))
cell_state_markers <- read.csv(list.files(data_path, full.names = TRUE, pattern = '*.csv'), row.names = 1) %>% select(!c(evidence))

cell_state_markers <- apply(cell_state_markers, 2, function(x) rownames(cell_state_markers)[x > 0])

cell_states = list(
  HH4 = c('NNE', 'node', 'streak', 'EE', 'eNPB', 'eN', 'eCN'),

  HH5 = c('NNE', 'node', 'streak', 'EE', 'eNPB', 'eN', 'eCN',
          'NPB', 'aNPB', 'pNPB', 'NP', 'pNP', 'iNP', 'aNP', 'PPR', 'aPPR', 'pPPR'),
  
  HH6 = c('NNE', 'node', 'streak', 'eN', 'eCN',
          'NPB', 'aNPB', 'pNPB', 'NP', 'pNP', 'iNP', 'aNP', 'PPR', 'aPPR', 'pPPR'),

  HH7 = c('pEpi', 'NPB', 'aNPB', 'pNPB', 'NC', 'dNC', 'NP', 'pNP', 'iNP',
          'aNP', 'HB', 'MB', 'FB', 'vFB', 'PPR', 'aPPR', 'pPPR'),

  ss4 = c('pEpi', 'NPB', 'aNPB', 'pNPB', 'NC', 'dNC', 'NP', 'pNP', 'iNP',
          'aNP', 'HB', 'MB', 'FB', 'vFB', 'PPR', 'aPPR', 'pPPR'),

  ss8 = c('pEpi', 'NPB', 'aNPB', 'pNPB', 'NC', 'dNC', 'NP', 'pNP', 'iNP',
          'aNP', 'HB', 'MB', 'FB', 'vFB', 'PPR', 'aPPR', 'pPPR')
)

cell_state_markers <- lapply(cell_states, function(x) cell_state_markers[names(cell_state_markers) %in% x])


# Run classification using different resolutions for different stages
stage = unique(seurat_data@meta.data$stage)

if(length(stage) == 1){
  cell_state_markers = cell_state_markers[[stage]]
  cluster_res = list(HH4 = 1.2, HH5 = 1.2, HH6 = 1.2, HH7 = 1.2, ss4 = 1.2, ss8 = 1.2)[[stage]]
  metadata = c('scHelper_cell_type')
} else {
  cell_state_markers = flatten(cell_state_markers)
  cell_state_markers = cell_state_markers[!duplicated(cell_state_markers)]
  metadata = c('scHelper_cell_type', 'stage')
  
  # Set cluster res to 2 for run split subsets - 3 for integrated data
  cluster_res = ifelse(length(unique(seurat_data$run)) == 1, 2, 3)
}

DefaultAssay(seurat_data) <- "integrated"

seurat_data <- FindClusters(seurat_data, resolution = cluster_res)

# Set RNA to default assay for plotting expression data
DefaultAssay(seurat_data) <- "RNA"

cell_type_df <- lapply(cell_state_markers, function(x) t(GetAssayData(object = seurat_data, assay = 'RNA', slot = 'scale.data'))[,x] %>% rowSums(.)) %>%
  do.call('cbind', .) %>%
  merge(., seurat_data@meta.data[,'seurat_clusters', drop=FALSE], by=0, all=TRUE)

cell_type_df <- cell_type_df %>% column_to_rownames('Row.names') %>%
  pivot_longer(cols = !seurat_clusters) %>%
  group_by(seurat_clusters, name)

# Plot cell classification per cluster
dir.create(paste0(plot_path, "scHelper_log/"))
png(paste0(plot_path, "scHelper_log/classification_boxplots.png"), width = 30, height = length(unique(cell_type_df$seurat_clusters)) * 3, units = 'cm', res = 200)
ggplot(cell_type_df, aes(x = name, y = value, fill = name)) +
  geom_boxplot() +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(cell_type_df$name)))) +
  facet_wrap(~seurat_clusters, ncol = 2) +
  xlab(NULL) + 
  ylab('Average scaled expression') +
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text = element_text(size = 10), 
        axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 10)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
graphics.off()

# Select top cell type per cluster
cell_type_df <- cell_type_df %>%
  summarise(value = mean(value)) %>%
  group_by(seurat_clusters) %>%
  filter(value == max(value))

# add cell_type to seurat metadata
seurat_data@meta.data[['scHelper_cell_type']] <- unlist(apply(seurat_data@meta.data, 
                                                              1, function(x) ifelse(x[['seurat_clusters']] %in% cell_type_df[['seurat_clusters']], 
                                                                                    cell_type_df[cell_type_df[['seurat_clusters']] == x[['seurat_clusters']], "name"], NA)))

# Old method
# seurat_data <- ClusterClassification(seurat_obj = seurat_data, cell_state_markers = cell_state_markers, force_assign = FALSE, quantile = 0.5, plot_path = paste0(plot_path, "scHelper_log/"))

#####################   Set levels and colours for scHelper_cell_type  ###########################################
scHelper_order <- scHelper_all_order[scHelper_all_order %in% unique(seurat_data@meta.data[["scHelper_cell_type"]])]
seurat_data@meta.data$scHelper_cell_type <- factor(seurat_data@meta.data$scHelper_cell_type, levels = scHelper_all_order)

scHelper_cols = scHelper_all_colours[unique(seurat_data@meta.data$scHelper_cell_type)]
names(scHelper_cols) <- NULL

#####################   Set levels and colours for stage   ###########################################
stage_order <- stage_all_order[stage_all_order %in% unique(seurat_data@meta.data[["stage"]])]
seurat_data@meta.data$stage <- factor(seurat_data@meta.data$stage, levels = stage_all_order)

stage_cols = stage_all_colours[levels(seurat_data@meta.data$stage)]
names(stage_cols) <- NULL

# Plot UMAP for clusters and developmental stage
png(paste0(plot_path, "ClustStagePlot_scHelper_celltype_umap.png"), width=24, height=12, units = 'cm', res = 200)
ClustStagePlot(seurat_data, stage_col = "stage", cluster_col = "scHelper_cell_type", label_clusters = TRUE)
graphics.off()

# UMAP for cell state
png(paste0(plot_path, "scHelper_celltype_umap.png"), width=12, height=12, units = 'cm', res = 200)
DimPlot(seurat_data, group.by = 'scHelper_cell_type', label = TRUE, 
        label.size = ifelse(length(unique(seurat_data$stage)) == 1, 5, 3),
        label.box = TRUE, repel = TRUE,
        pt.size = ifelse(length(unique(seurat_data$stage)) == 1, 2, 1), 
        cols = scHelper_cols) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

# UMAP for stage
png(paste0(plot_path, "stage_umap.png"), width=12, height=12, units = 'cm', res = 200)
DimPlot(seurat_data, group.by = 'stage', label = TRUE, 
        label.size = ifelse(length(unique(seurat_data$stage)) == 1, 5, 3),
        label.box = TRUE, repel = TRUE,
        pt.size = ifelse(length(unique(seurat_data$stage)) == 1, 2, 1), 
        cols = stage_cols) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

saveRDS(seurat_data, paste0(rds_path, label, "_cell_state_classification.RDS"), compress = FALSE)


# Plot stacked violins for each of the cell type classes to check genes used are good markers
curr_plot_path = paste0(plot_path, "cell_type_dotplots/")
dir.create(curr_plot_path)

for(i in names(cell_state_markers)){
  png(paste0(curr_plot_path, i, ".png"), width = (length(cell_state_markers[[i]])+2)*3, height = 15, units = 'cm', res = 200)
  print(DotPlot(seurat_data, features = cell_state_markers[[i]], group.by = "seurat_clusters"))
  graphics.off()
}


# Plot multi-feature plots for each of the cell type classes
curr_plot_path = paste0(plot_path, "cell_type_feature_plots/")
dir.create(curr_plot_path)

for(i in names(cell_state_markers)){
  ncol = ceiling((length(cell_state_markers[[i]])+1)/8)+1
  nrow = ceiling((length(cell_state_markers[[i]])+1)/ncol)

  png(paste0(curr_plot_path, i, '.png'), width = ncol*10, height = nrow*10, units = "cm", res = 200)
  MultiFeaturePlot(seurat_data, plot_stage = TRUE, stage_col = "stage", plot_celltype = TRUE, celltype_col = "seurat_clusters",
                   gene_list = cell_state_markers[[i]], n_col = ncol, label = '')
  graphics.off()
}