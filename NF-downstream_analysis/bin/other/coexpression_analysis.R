library(optparse)
library(future)
library(Seurat)
library(pheatmap)
library(tidyverse)
library(Antler)
library(RColorBrewer)
library(scHelper)


## Set paths ( this will need to be generalised )
plot_path = "./Eva_working/HH6/plots/"
rds_path = "./Eva_working/HH6/rds_files/"
antler_path = "./Eva_working/HH6/antler_data/"

## read in seurat object and antler data (will need to set correct names/paths)
seurat_data <- readRDS(paste0(rds_path, "hh6_cell_state_classification.RDS"))
DefaultAssay(seurat_data) <- "RNA"
antler_data <- readRDS(paste0(rds_path, "/antler_out.RDS"))

# set metadata for ordering (make generalisable?)
metadata <- c("scHelper_cell_type", "run")


### create cell_order object based on ordering from metadata
class_order = c('extra_embryonic', 'early_non_neural', 'non_neural', 'early_NNE', 'early_PPR', 'early_aPPR', 'aPPR', 'iPPR',
                'early_pPPR', 'pPPR', 'early_border', 'early_NPB', 'NPB', 'early_aNPB', 'aNPB', 'early_pNPB', 'pNPB', 'NC',
                'delaminating_NC', 'early_neural', 'early_neural_plate', 'early_caudal_neural', 'neural_progenitors', 'p_neural_progenitors',
                'early_hindbrain', 'hindbrain', 'early_midbrain', 'midbrain', 'a_neural_progenitors', 'early_forebrain', 'forebrain', 'node')
seurat_data@meta.data$scHelper_cell_type <- factor(seurat_data@meta.data$scHelper_cell_type, levels = class_order)
col_ann <- seurat_data@meta.data[,metadata, drop=FALSE] %>% mutate_if(is.character, as.factor)
col_ann <- col_ann[do.call("order", c(col_ann[metadata], list(decreasing = FALSE))), , drop = FALSE]
cell_order = rownames(col_ann)

################################################################################################################################
# function to extract gm counts, calculate aggregate gene expression and calculate coexpression as a product of those aggregates
GM_expression <- function(seurat_data = seurat_data, GM_1_genes, GM_2_genes, 
                          cell_order = cell_order, assay = "RNA", slot = "data"){
  
  # extract counts and order cells based on cell_order
  counts <- t(GetAssayData(object = seurat_data, assay = assay, slot = slot))
  counts <- counts[cell_order,]
  
  # extract gene counts and calculate aggregate expression per module
  GM_1 <- counts[, c(GM_1_genes)]
  GM_1 <- rowSums(GM_1)
  GM_1 <- data.frame(order = 1:length(GM_1), expression = GM_1)
  GM_2 <- counts[, c(GM_2_genes)]
  GM_2 <- rowSums(GM_2)
  GM_2 <- data.frame(order = 1:length(GM_2), expression = GM_2)
  
  # merge and calculate coexpression as product
  merge <- merge(GM_1, GM_2, by = "order")
  colnames(merge) <- c("Cells", "Gene module A", "Gene module B")
  merge$Coexpression <- merge$`Gene module A` * merge$`Gene module B`
  
  return(merge)
}
################################################################################################################################

## choose which gms want to plot from antler_data (make generalisable)
GM_1_genes <- antler_data$gene_modules$lists$unbiasedGMs_DE$content$GM17
GM_2_genes <- antler_data$gene_modules$lists$unbiasedGMs_DE$content$GM4

## run GM_expression function
merge <- GM_expression(seurat_data, cell_order = cell_order,
                       GM_1_genes = GM_1_genes, 
                       GM_2_genes = GM_2_genes)


## plot output of GM_expression
plot_merge <- pivot_longer(merge,
                     c("Gene module A", "Gene module B", "Coexpression"), 
                     names_to = "group_by", values_to = "Expression")

gm1 <- ggplot(data=plot_merge[which(plot_merge$group_by == "Gene module A"), ], aes(x=Cells, y=Expression)) +
  geom_smooth(method = "gam", se = FALSE) +
  facet_grid(group_by ~.) +
  geom_point(alpha = 0.1) +
  theme_light()

gm2 <- ggplot(data=plot_merge[which(plot_merge$group_by == "Gene module B"), ], aes(x=Cells, y=Expression)) +
  geom_smooth(method = "gam", se = FALSE) +
  facet_grid(group_by ~.) +
  geom_point(alpha = 0.1) +
  theme_light()

coexpression <- ggplot(data=plot_merge[which(plot_merge$group_by == "Coexpression"), ], aes(x=Cells, y=Expression)) +
  geom_smooth(method = "gam", se = FALSE) +
  facet_grid(group_by ~.) +
  geom_point(alpha = 0.1) +
  theme_light()

grid.arrange(gm1, gm2, coexpression)