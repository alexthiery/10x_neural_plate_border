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

metadata <- read.csv(list.files(data_path, pattern = "*.csv", full.names = TRUE))
# metadata <- read.csv('./output/NF-downstream_analysis_stacas/filtered_seurat/cellrank/NF-scRNAseq_alignment_out_metadata.csv')

seurat_data <- readRDS(list.files(data_path, pattern = "*.RDS", full.names = TRUE)[!list.files(data_path, pattern = "*.RDS") %>% grepl('antler', .)])
seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/clusters_subset/late_stage_subset/seurat/clusters_cluster/rds_files/late_clustered_data.RDS')

# load antler data
antler_data <- readRDS(list.files(data_path, pattern = "antler_out.RDS", full.names = TRUE))
antler_data <- readRDS('./output/NF-downstream_analysis_stacas/clusters_subset/late_stage_subset/antler/clusters_gene_modules/rds_files/antler_out.RDS')



metadata$CellID <- paste0(metadata$CellID, "-1")

# Rename cellIDs to match seurat data based on string match
new_names <- unlist(lapply(metadata$CellID, function(x) rownames(seurat_data@meta.data)[grep(x, rownames(seurat_data@meta.data))]))

if(length(new_names) != nrow(metadata) | length(new_names) != nrow(seurat_data@meta.data)){stop('cell IDs differ between scvelo metadata and seurat object')}

rownames(metadata) <- new_names
metadata$CellID <- NULL

# re-order metadata rows based on seurat metadata order
metadata <- metadata[ order(match(rownames(metadata), rownames(seurat_data@meta.data))), ]

# replace seurat metadata with scvelo metadata
seurat_data@meta.data <- metadata

# set boolean for whether dataset contains multiple batches
multi_run <- ifelse(length(unique(seurat_data$run)) > 1, TRUE, FALSE)


#####################################################################################################
#                           calculate module scores                   #
#####################################################################################################
# access DE gene modules (batchfilt if there are multiple batches in the dataset)
if(is.null(antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt)){
  gms <- antler_data$gene_modules$lists$unbiasedGMs_DE$content
}else{
  gms <- antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content
}

if(is.null(names(gms))){names(gms) = paste0("GM: ", 1:length(gms))}


# Set RNA to default assay for plotting expression data
DefaultAssay(seurat_data) <- "RNA"


coexpression = function(seurat_object, gm_1, gm_2){

  x = seurat_object@meta.data[,'scHelper_cell_type']
  
  gm_1_sum <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data'))[,gm_1] %>% rowSums(.)
  
  # gm_1_sum = gm_1_sum/max(gm_1_sum)
  
  gm_2_sum <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data'))[,gm_2] %>% rowSums(.)
  
  # gm_2_sum = gm_2_sum/max(gm_2_sum)
  

  # gm_1 <- unlist(gms[gm_1])
  # gm_1_sum <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data'))[,gm_1] %>% rowSums(.)
  # 
  # gm_2 <- unlist(gms[gm_2])
  # gm_2_sum <- t(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data'))[,gm_2] %>% rowSums(.)
  
  # gm_1_sum[gm_1_sum < 0] <- 0
  # gm_2_sum[gm_2_sum < 0] <- 0
  
  gm_product <- gm_1_sum * gm_2_sum
  plot_data <- data.frame(gm_1_sum = gm_1_sum, gm_2_sum = gm_2_sum, gm_product = gm_product, x = x)
  
  plot_data <- plot_data %>% mutate(ratio = gm_1_sum/(gm_1_sum + gm_2_sum)) %>% arrange(ratio)
  
  
  # plot_data <- plot_data %>% arrange(x)
  
  if(class(plot_data$x) != 'numeric'){
    plot_data$x <- 1:nrow(plot_data)
  }
  
  p1 = ggplot(plot_data, aes(x = x, y = gm_1_sum)) +
    # geom_point() +
    geom_smooth(method = "gam", se = FALSE)
  
  p2 = ggplot(plot_data, aes(x = x, y = gm_2_sum)) +
    # geom_point() +
    geom_smooth(method = "gam", se = FALSE)
  
  p3 = ggplot(plot_data, aes(x = x, y = gm_product)) +
    # geom_point() +
    geom_smooth(method = "gam", se = FALSE)
  
  return(p1/p2/p3)
}

scHelper_cell_type_order <- c('extra_embryonic', 'early_NNE', 'NNE', 'prospective_epidermis', 'PPR', 'aPPR',
                              'pPPR', 'early_border', 'early_NPB', 'NPB', 'early_aNPB', 'aNPB', 'early_pNPB', 'pNPB', 'NC', 'delaminating_NC', 'early_neural',
                              'early_neural_plate', 'early_caudal_neural', 'neural_progenitors', 'p_neural_progenitors', 'early_hindbrain', 'hindbrain',
                              'early_midbrain', 'midbrain', 'a_neural_progenitors', 'early_forebrain', 'forebrain', 'a_ventral_floorplate', 'streak', 'node')


# seurat_data <- subset(seurat_data, cells = rownames(seurat_data@meta.data)[seurat_data@meta.data[['scHelper_cell_type']] %in% c('aPPR', 'pPPR', 'aNPB', 'pNPB', 'NC')])

seurat_data@meta.data$scHelper_cell_type <- factor(seurat_data@meta.data$scHelper_cell_type, levels = scHelper_cell_type_order)

coexpression(seurat_data, gm_1 = unlist(gms['GM8']), gm_2 = unlist(gms['GM40']))


PPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "PRDM1")
NC = c("PAX7", "MSX1", "MSX2", "ETS1", "ENSGALG00000030902", "FOXD3", "TFAP2B")


MultiFeaturePlot(seurat_data, NC)


coexpression(ss4_npb, gm_1 = PPR, gm_2 = NC)
coexpression(ss8_npb, gm_1 = PPR, gm_2 = NC)

test <- merge(ss4_npb, ss8_npb)


coexpression(test, gm_1 = PPR, gm_2 = NC)


ss4_npb <- subset(seurat_data, cells = rownames(seurat_data@meta.data)[seurat_data@meta.data[['stage']] == 'ss4'])
ss4_npb <- subset(ss4_npb, cells = rownames(seurat_data@meta.data)[seurat_data@meta.data[['scHelper_cell_type']] %in% c('aPPR', 'pPPR', 'NC', 'delaminating_NC')])
ss4_npb@meta.data$scHelper_cell_type <- factor(ss4_npb@meta.data$scHelper_cell_type, levels = scHelper_cell_type_order)

coexpression(ss4_npb, gm_1 = unlist(gms['GM10']), gm_2 = unlist(gms['GM40']))

ss8_npb <- subset(seurat_data, cells = rownames(seurat_data@meta.data)[seurat_data@meta.data[['stage']] == 'ss8'])
ss8_npb <- subset(ss8_npb, cells = rownames(seurat_data@meta.data)[seurat_data@meta.data[['scHelper_cell_type']] %in% c('PPR', 'pPPR', 'NC', 'delaminating_NC')])
ss8_npb@meta.data$scHelper_cell_type <- factor(ss8_npb@meta.data$scHelper_cell_type, levels = scHelper_cell_type_order)

coexpression(ss8_npb, gm_1 = unlist(gms['GM8']), gm_2 = unlist(gms['GM1']))



# # Set gene module order
# stage_order <- c("hh4", "hh5", "hh6", "hh7", "ss4", "ss8")
# scHelper_cell_type_order <- c('extra_embryonic', 'early_non_neural', 'non_neural', 'early_NNE', 'early_PPR', 'early_aPPR', 'aPPR', 'iPPR',
#                               'early_pPPR', 'pPPR', 'early_border', 'early_NPB', 'NPB', 'early_pNPB', 'pNPB', 'early_aNPB', 'aNPB', 'early_neural',
#                               'early_neural_plate', 'early_caudal_neural', 'neural_progenitors', 'a_neural_progenitors', 'early_forebrain', 'forebrain',
#                               'early_midbrain', 'midbrain', 'p_neural_progenitors', 'early_hindbrain', 'hindbrain', 'NC', 'delaminating_NC', 'node')
# 
# 
# ### create cell_order object based on ordering from metadata
# class_order = c('extra_embryonic', 'early_non_neural', 'non_neural', 'early_NNE', 'early_PPR', 'early_aPPR', 'aPPR', 'iPPR',
#                 'early_pPPR', 'pPPR', 'early_border', 'early_NPB', 'NPB', 'early_aNPB', 'aNPB', 'early_pNPB', 'pNPB', 'NC',
#                 'delaminating_NC', 'early_neural', 'early_neural_plate', 'early_caudal_neural', 'neural_progenitors', 'p_neural_progenitors',
#                 'early_hindbrain', 'hindbrain', 'early_midbrain', 'midbrain', 'a_neural_progenitors', 'early_forebrain', 'forebrain', 'node')
# seurat_data@meta.data$scHelper_cell_type <- factor(seurat_data@meta.data$scHelper_cell_type, levels = class_order)
# col_ann <- seurat_data@meta.data[,metadata, drop=FALSE] %>% mutate_if(is.character, as.factor)
# col_ann <- col_ann[do.call("order", c(col_ann[metadata], list(decreasing = FALSE))), , drop = FALSE]
# cell_order = rownames(col_ann)

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