#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
library(future)
library(Seurat)
library(pheatmap)
library(tidyverse)
library(Antler)
library(RColorBrewer)
library(scHelper)
library(gridExtra)
library(grid)

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

metadata <- read.csv(list.files(data_path, pattern = "*.csv", full.names = TRUE))
# metadata <- read.csv('./output/NF-downstream_analysis_stacas/filtered_seurat/cellrank/NF-scRNAseq_alignment_out_metadata.csv')

seurat_data <- readRDS(list.files(data_path, pattern = "*.RDS", full.names = TRUE)[!list.files(data_path, pattern = "*.RDS") %>% grepl('antler', .)])
# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/filtered_seurat/seurat/state_classification/rds_files/contamination_cell_state_classification.RDS')

# load antler data
antler_data <- readRDS(list.files(data_path, pattern = "antler_out.RDS", full.names = TRUE))
# antler_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/hh5_splitstage_data/antler/stage_gene_modules/rds_files/antler_out.RDS')



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

# Set gene module order
stage_order <- c("hh4", "hh5", "hh6", "hh7", "ss4", "ss8")
scHelper_cell_type_order <- c('extra_embryonic', 'early_non_neural', 'non_neural', 'early_NNE', 'early_PPR', 'early_aPPR', 'aPPR', 'iPPR',
                              'early_pPPR', 'pPPR', 'early_border', 'early_NPB', 'NPB', 'early_pNPB', 'pNPB', 'early_aNPB', 'aNPB', 'early_neural',
                              'early_neural_plate', 'early_caudal_neural', 'neural_progenitors', 'a_neural_progenitors', 'early_forebrain', 'forebrain',
                              'early_midbrain', 'midbrain', 'p_neural_progenitors', 'early_hindbrain', 'hindbrain', 'NC', 'delaminating_NC', 'node')

# ordered_gms = GMOrder(seurat_obj = seurat_data, gms, metadata_1 = 'stage', order_1 = stage_order,
#                       metadata_2 = 'scHelper_cell_type', order_2 = scHelper_cell_type_order)
# 


plot_data <- data.frame()
for(mod_name in names(gms)){
  temp <- GetAssayData(seurat_data, assay = 'RNA', slot = 'scale.data')[gms[[mod_name]],]
  
  temp <- merge(t(temp), seurat_data@meta.data[,c('latent_time', 'lineage_NC_probability', 'lineage_neural_probability', 'lineage_placodal_probability'), drop=FALSE], by=0)
  plot_data <- temp %>%
    column_to_rownames('Row.names') %>%
    pivot_longer(!c(latent_time, lineage_NC_probability, lineage_neural_probability, lineage_placodal_probability)) %>%
    rename(scaled_expression = value) %>%
    rename(gene = name) %>%
    pivot_longer(cols = !c(latent_time, gene, scaled_expression)) %>%
    rename(lineage = name) %>%
    mutate(module = mod_name) %>%
    bind_rows(plot_data)
}

# set colours
colours <- setNames(colorRampPalette(brewer.pal(9, "Paired"))(length(gms)), names(gms))

for(mod_name in names(gms)){
  plot = ggplot(filter(plot_data, module == mod_name), aes(x = latent_time, y = scaled_expression, color = module)) +
    scale_color_manual(values=colours[[mod_name]]) +
    geom_smooth(method="gam", se=FALSE, mapping = aes(weight = value, linetype = lineage, group=lineage)) +
    xlab("Latent time") + ylab("Scaled expression") +
    # facet_wrap(~lineage, dir = 'v') +
    theme_classic()
  
  png(paste0(plot_path, mod_name, '_lineage_dynamics.png'), width = 25, height = 15, units='cm', res=200)
  print(plot)
  graphics.off()
  
  # Calculate average module expression for plotting multi-feature plot
  seurat_data@meta.data[[mod_name]] <-  colMeans(GetAssayData(seurat_data, assay = 'RNA', slot = 'scale.data')[gms[[mod_name]],])
}

ncol = ceiling((length(names(gms))+1)/8)+1
nrow = ceiling((length(names(gms))+1)/ncol)

png(paste0(plot_path, 'gene_module_feature.png'), width = ncol*10, height = nrow*10, units = "cm", res = 200)
MultiFeaturePlot(seurat_object = seurat_data, gene_list = grep('GM', colnames(seurat_data@meta.data), value = TRUE), plot_stage = TRUE,
                 stage_col = 'stage', plot_clusters = FALSE, n_col = ncol)
graphics.off()



# # Calculate average module expression for contamination gene list
# seurat_data <- AverageGeneModules(seurat_obj = seurat_data, gene_list = gms)

# temp <- GetAssayData(seurat_data, assay = 'RNA', slot = 'scale.data')[gms[[1]],]
# 
# temp <- merge(t(temp), seurat_data@meta.data[,c('latent_time', 'lineage_NC_probability', 'lineage_neural_probability', 'lineage_placodal_probability'), drop=FALSE], by=0)
# 
# plot_data <- temp %>%
#   column_to_rownames('Row.names') %>%
#   pivot_longer(!c(latent_time, lineage_NC_probability, lineage_neural_probability, lineage_placodal_probability)) %>%
#   rename(scaled_expression = value) %>%
#   rename(gene = name) %>%
#   pivot_longer(cols = !c(latent_time, gene, scaled_expression)) %>%
#   rename(lineage = name)









