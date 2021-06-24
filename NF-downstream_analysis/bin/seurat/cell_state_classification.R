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
    
    plot_path = "./output/NF-downstream_analysis_stacas/seurat/cell_state_classification/plots/"
    rds_path = "./output/NF-downstream_analysis_stacas/seurat/cell_state_classification/rds_files/"
    data_path = "./output/NF-downstream_analysis_stacas/seurat/6_contamination_filt/rds_files/"
    
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

seurat_data <- readRDS(list.files(data_path, full.names = TRUE))


########################################################################################################
#                                      Cell state classification                                    #
########################################################################################################

DimPlot(seurat_data)



###################### Nerual crest and placode cell type identification ###################### 
placNC_genes <- c(
  # delaminating NC
  "ETS1", "SOX10", "SOX8", "LMO4",  "TFAP2B", "SOX9",
  # NPB
  "TFAP2A", "DRAXIN", "MSX1", "CSRNP1", "PAX7", "BMP5", "MSX2",
  # NC
  "WNT6",
  # Placodes
  "PITX1", "PITX2", "ZNF385C",  "SIX1", "EYA2", "DLX6", "HOMER2"
)

ncol = ceiling((length(placNC_genes)+1)/8)+1
nrow = ceiling((length(placNC_genes)+1)/ncol)

# plot expression of NC and placodal genes
png(paste0(plot_path, 'multi_feature_plac_nc.png'), width = ncol*10, height = nrow*10, units = "cm", res = 200)
MultiFeaturePlot(seurat_data, gene_list = placNC_genes, n_col = ncol, label = '')
graphics.off()


# Generate pie charts for cluster dev stage composition
venn_data <- seurat_data@meta.data %>%
  rownames_to_column('cell_name') %>%
  dplyr::select(c(cell_name, stage, seurat_clusters)) %>%
  group_by(seurat_clusters) %>%
  count(stage, .drop = FALSE)

venn_data <- venn_data %>%
  mutate(total_cells = sum(n)) %>%
  mutate(n = n/total_cells)

# Generate DotPlot
dotplot <- DotPlot(seurat_data, group.by = "seurat_clusters", features = placNC_genes) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position="bottom", legend.box = "horizontal",
        axis.title.x=element_blank(), legend.text=element_text(size=10),
        legend.title=element_text(size=12))

# Generate Pie charts
pies <- ggplot(venn_data, aes(x=as.numeric(total_cells)/2, y=as.numeric(n), fill=stage, width = total_cells)) +
  geom_bar(position = 'fill', stat = "identity") +
  facet_wrap( ~ seurat_clusters, nrow = nlevels(seurat_data@meta.data$seurat_clusters)) +
  coord_polar("y") +
  theme_void() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        legend.position=c(0,0), legend.direction = "horizontal",
        plot.margin = margin(t = 14, r = 0, b = 119, l = -110, unit = "pt"),
        legend.margin=margin(t = 70, r = 0, b = -100, l = -200, unit = "pt"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12)) +
  labs(fill = "Stage")

# Plot dotplot and pies
png(paste0(plot_path, "dotplot.png"), width = 32, height = 18, units = "cm", res = 200)
print(plot_grid(dotplot, pies, rel_widths = c(5,1)))
graphics.off()



# add neural crest and placodal cell types to metadata
norm.data.clustfilt.cc@meta.data$placNC_clust <- apply(norm.data.clustfilt.cc@meta.data, 1, function(x)
  if(x["seurat_clusters"] == 10){"Delaminating NC"} else if(x["seurat_clusters"] == 12){"NC Progenitors"}
  else if(x["seurat_clusters"] == 9){"Placodes"}else {NA})

# plot annotated neural crest and placodal clusters
png(paste0(plot.path, "plac.clust.png"), width = 13, height = 10, res = 200, units = "cm")
DimPlot(norm.data.clustfilt.cc, group.by = "placNC_clust")  + ggtitle("Clusters") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

# plot dotplot for neural crest and placodal genes and clusters
png(paste0(plot.path, "plac.dotplot.png"), width = 25, height = 10, res = 200, units = "cm")
DotPlot(norm.data.clustfilt.cc[, !is.na(norm.data.clustfilt.cc@meta.data$placNC_clust)], group.by = "placNC_clust", features = placNC_genes) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()







# Genes of interest used for cell state classification
GOI = rev(c( "EOMES", "ADMP","YEATS4", "MAFA", "ING5", "LIN28B", "AATF","SETD2",
             "GATA2", "DLX5", "TFAP2A", "BMP4", "SIX1", "EYA2",
             "MSX1", "PAX7", "CSRNP1", "SOX10",
             "SOX2", "SOX21", "SIX3", "OLIG2", "PAX6", "FOXA2", "SHH", "PAX2", "WNT4", "HOXB2", "HOXA2", "GBX2"))

# Change order or clusters for plotting dotplots
cluster_order <- factor(norm.data.clustfilt.cc$seurat_clusters, levels = c(12,5,2,1,
                                                                           3,8,0,
                                                                           11, 10,
                                                                           7,4,14,
                                                                           9,6,13))

# Set factor levels for plotting
norm.data.clustfilt.cc$seurat_clusters <- cluster_order

# Generate pie charts for cluster dev stage composition
venn_data <- norm.data.clustfilt.cc@meta.data %>%
  rownames_to_column('cell_name') %>%
  dplyr::select(c(cell_name, orig.ident, seurat_clusters)) %>%
  group_by(seurat_clusters) %>%
  count(orig.ident, .drop = FALSE)

venn_data <- venn_data %>%
  mutate(total_cells = sum(n)) %>%
  mutate(n = n/total_cells)


# Reverse levels to deal with seurat dotplot reversing y axis
norm.data.clustfilt.cc$seurat_clusters <- fct_rev(cluster_order)

# Generate DotPlot
dotplot <- DotPlot(norm.data.clustfilt.cc, group.by = "seurat_clusters", features = GOI) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position="bottom", legend.box = "horizontal",
        axis.title.x=element_blank(), legend.text=element_text(size=10),
        legend.title=element_text(size=12))

# Generate Pie charts
pies <- ggplot(venn_data, aes(x=as.numeric(total_cells)/2, y=as.numeric(n), fill=orig.ident, width = total_cells)) +
  geom_bar(position = 'fill', stat = "identity") +
  facet_wrap( ~ seurat_clusters, nrow = nlevels(norm.data.clustfilt.cc@meta.data$seurat_clusters)) +
  coord_polar("y") +
  theme_void() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        legend.position=c(0,0), legend.direction = "horizontal",
        plot.margin = margin(t = 14, r = 0, b = 119, l = -110, unit = "pt"),
        legend.margin=margin(t = 70, r = 0, b = -100, l = -200, unit = "pt"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12)) +
  labs(fill = "Stage")

# Plot dotplot and pies
png(paste0(curr.plot_path, "dotplot.png"), width = 32, height = 18, units = "cm", res = 200)
print(plot_grid(dotplot, pies, rel_widths = c(5,1)))
graphics.off()


# Plot dotplot with cell classification labels
dotplot <- DotPlot(norm.data.clustfilt.cc, group.by = "seurat_clusters", features = GOI) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position="bottom", legend.box = "horizontal",
        axis.title.x=element_blank(), legend.text=element_text(size=10),
        axis.title.y = element_blank(), 
        legend.title=element_text(size=12)) +
  scale_y_discrete(labels = rev(c('Node', 'HH4-1', 'HH4-2', 'HH4-3', 'Non-Neural Plate Progenitors',
                                  'Placodes', 'Neural Plate Progenitors', 'NC Progenitors', 'Delaminating NC',
                                  'Early Forebrain', 'Late Forebrain', 'Ventral Forebrain', 'Early Midbrain', 'Late Midbrain', 'Hindbrain')))

png(paste0(curr.plot_path, "dotplot_classified.png"), width = 36, height = 18, units = "cm", res = 200)
print(plot_grid(dotplot, pies, rel_widths = c(5,1)))
graphics.off()
