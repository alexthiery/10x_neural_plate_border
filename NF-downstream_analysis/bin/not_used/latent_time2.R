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

# Multi-core when running from command line
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 16* 1024^3) # 32gb

cat(paste0("script ran with ", ncores, " cores\n"))
dir.create(plot_path, recursive = T)
dir.create(rds_path, recursive = T)

# metadata <- read.csv(list.files(data_path, pattern = "*.csv", full.names = TRUE))
metadata <- read.csv('./output/NF-downstream_analysis_stacas/filtered_seurat/cellrank/NF-scRNAseq_alignment_out_metadata.csv')
metadata <- read.csv('./output/NF-downstream_analysis_stacas/transfer_subset/cellrank/transfer_filter/all_stages_filtered_metadata.csv')

metadata$CellID <- paste0(metadata$CellID, "-1")

# seurat_data <- readRDS(list.files(data_path, pattern = "*.RDS", full.names = TRUE)[!list.files(data_path, pattern = "*.RDS") %>% grepl('antler', .)])
seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/filtered_seurat/seurat/state_classification/rds_files/contamination_cell_state_classification.RDS')
seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/transfer_subset/seurat/transfer_filter/rds_files/transfer_filter.RDS')

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

# load antler data
# antler_data <- readRDS(list.files(data_path, pattern = "antler_out.RDS", full.names = TRUE))
antler_data <- readRDS('./output/NF-downstream_analysis_stacas/filtered_seurat/antler/gene_modules/rds_files/antler_out.RDS')

#####################################################################################################
#                           calculate module scores                   #
#####################################################################################################

# access DE gene modules (batchfilt if there are multiple batches in the dataset)
if(multi_run){gms <- antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content}else{gms <- antler_data$gene_modules$lists$unbiasedGMs_DE$content}
if(is.null(names(gms))){names(gms) = paste0("GM: ", 1:length(gms))}

# Set RNA to default assay for plotting expression data
DefaultAssay(seurat_data) <- "RNA"

# Calculate average module expression for contamination gene list
seurat_data <- AverageGeneModules(seurat_obj = seurat_data, gene_list = gms)

plot_data <- seurat_data@meta.data[,c('latent_time', names(gms))] %>%
  pivot_longer(cols = !latent_time) %>%
  rename(module = name) %>%
  rename(module_score = value)


# Mean and SE summary data
plot_data_summary <- plot_data %>%
  mutate(rank_bin = latent_time - (latent_time %% 0.025)) %>%
  group_by(rank_bin, module) %>% 
  summarise(mn = mean(module_score),
            se = sd(module_score)/sqrt(n()))


# Plot GAM for module score without standard error
png(paste0(plot_path, 'gam_gms_latent_time.png'), height = 18, width = 26, units = 'cm', res = 400)
ggplot(plot_data, aes(x = latent_time, y = module_score, colour = module)) +
  geom_smooth(method="gam", se=FALSE) +
  xlab("Latent time") + ylab("Gene module score") +
  theme_classic()
graphics.off()

# Plot GAM for module score with standard error
png(paste0(plot_path, 'gam_se_gms_latent_time.png'), height = 18, width = 26, units = 'cm', res = 400)
ggplot(plot_data, aes(x = latent_time, y = module_score, colour = module)) +
  geom_errorbar(data = plot_data_summary, aes(x = rank_bin, y = mn, ymax = mn + se, ymin = mn - se), width = 0.01) +
  geom_point(data = plot_data_summary, aes(x = rank_bin, y = mn)) +
  geom_smooth(method="gam", se=FALSE) +
  xlab("Latent time") + ylab("Gene module score") +
  theme_classic()
graphics.off()


ncol = ceiling((length(names(gms))+1)/8)+1
nrow = ceiling((length(names(gms))+1)/ncol)

png(paste0(plot_path, 'multi_feature_module_score.png'), width = ncol*10, height = nrow*10, units = "cm", res = 200)
MultiFeaturePlot(seurat_data, gene_list = names(gms), n_col = ncol, label = 'UMAPs showing gene module scores')
graphics.off()

# 
# 
# ############# Plot latent time for selected custom modules
# 
# tissues <- list("Placodal" = c("GM: 7", "GM: 8", "GM: 9"),
#           "Progenitor" = c("GM: 10", "GM: 11", "GM: 13"),
#           "Neural" = c("GM: 36", "GM: 37", "GM: 39", "GM: 43"),
#           "Neural crest" = c("GM: 40", "GM: 41"))
# 
# colours <- c("Placodal" = "#7D08AC", "Progenitor" = "#32CC28", "Neural" = "#F9E300", "Neural crest" = "#2A50DF")
# 
# temp_plot_data <- plot_data %>%
#   filter(module %in% unlist(tissues)) %>%
#   mutate(tissue = unlist(reverseSplit(tissues)[module])) %>%
#   mutate(colour = colours[tissue]) %>%
#   mutate(tissue = factor(tissue, levels=c("Placodal", "Progenitor", "Neural", "Neural crest")))
# 
# # Extract columns for plotting
# colours = temp_plot_data %>% dplyr::pull(colour, module)
# colours = colours[!duplicated(names(colours))]
# 
# # Mean and SE summary data
# plot_data_summary <- temp_plot_data %>%
#   mutate(rank_bin = latent_time - (latent_time %% 0.025)) %>%
#   group_by(rank_bin, module) %>% 
#   summarise(mn = mean(module_score),
#             se = sd(module_score)/sqrt(n()))
# 
# # Plot GAM for module score without standard error
# png(paste0(plot_path, 'select_mod_latent_time.png'), height = 20, width = 15, units = 'cm', res = 400)
# ggplot(temp_plot_data, aes(x = latent_time, y = module_score, colour = module)) +
#   scale_color_manual(values=colours) +
#   geom_smooth(method="gam", se=FALSE) +
#   xlab("Latent time") + ylab("Gene module score") +
#   facet_grid(rows = vars(tissue)) +
#   theme_classic() +
#   theme(legend.position = "none")
# graphics.off()
# 
# 


test = seurat_data@meta.data[,c('latent_time', 'lineage_dNC_probability', 'lineage_neural_probability', 'lineage_placodal_probability', names(gms))] %>%
  pivot_longer(cols = !c(latent_time, lineage_dNC_probability, lineage_neural_probability, lineage_placodal_probability)) %>%
  rename(module = name) %>%
  rename(module_score = value) %>%
  filter(module == 'GM9') %>%
  pivot_longer(cols = !c(latent_time, module, module_score)) %>%
  rename(lineage = name)

plot_data <- seurat_data@meta.data[,c('latent_time', names(gms))] %>%
  pivot_longer(cols = !latent_time) %>%
  rename(module = name) %>%
  rename(module_score = value) %>%
  filter(module == 'GM26')


grep('lineage', colnames(seurat_data@meta.data), value=TRUE)


library('mgcv')

gam(~latent_time, data = seurat_data@meta.data %>% dplyr::select(latent_time, GM13))


# Plot GAM for module score without standard error
png(paste0(plot_path, 'gam_gms_latent_time.png'), height = 18, width = 26, units = 'cm', res = 400)
ggplot(test, aes(x = latent_time, y = module_score, colour = lineage, group = lineage)) +
  geom_smooth(method="gam", se=FALSE, mapping = aes(weight = value)) +
  xlab("Latent time") + ylab("Gene module score") +
  theme_classic()
graphics.off()







