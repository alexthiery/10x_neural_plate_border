


library(Seurat)
library(scHelper)
library(tidyverse)



seurat_data <- readRDS('~/output/NF-downstream_analysis_stacas/stage_split/hh6_splitstage_data/seurat/stage_state_classification/rds_files/cell_state_classification.RDS')

metadata <- read.csv('~/output/NF-downstream_analysis_stacas/scvelo/hh6_splitstage_data/scvelo_run/hh6_splitstage_data_metadata.csv', row.names = 1)



# edit names based on match between scvelo output and seurat data
rownames(metadata) <- sapply(rownames(metadata), function(x) rownames(seurat_data@meta.data)[grepl(x, rownames(seurat_data@meta.data))] )


# re-order metadata rows based on seurat metadata order
metadata <- metadata[ order(match(rownames(metadata), rownames(seurat_data@meta.data))), ]

# replace seurat metadata with scvelo metadata
seurat_data@meta.data <- metadata


DimPlot(seurat_data, group.by = 'clusters_gradients')

FeaturePlot(seurat_data, features = 'latent_time')

DimPlot(seurat_data, group.by = 'seurat_clusters')

DimPlot(seurat_data, group.by = 'scHelper_cell_type')


# set boolean for whether dataset contains multiple batches
multi_run <- ifelse(length(unique(seurat_data$run)) > 1, TRUE, FALSE)

# load antler data
# antler_data <- readRDS(list.files(data_path, pattern = "antler_out.RDS", full.names = TRUE))
antler_data <- readRDS('./output/NF-downstream_analysis_stacas/stage_split/hh6_splitstage_data/antler/stage_gene_modules/rds_files/antler_out.RDS')




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



seurat_data@meta.data$initial_states

ggplot(seurat_data@meta.data, aes(y=))



plot_data <- seurat_data@meta.data[,c('latent_time', 'GM16', 'clusters_gradients')] %>%
  pivot_longer(cols = !c(latent_time, clusters_gradients)) %>%
  rename(module = name) %>%
  rename(module_score = value)


# Mean and SE summary data
plot_data_summary <- plot_data %>%
  mutate(rank_bin = latent_time - (latent_time %% 0.025)) %>%
  group_by(rank_bin, module, clusters_gradients) %>% 
  summarise(mn = mean(module_score),
            se = sd(module_score)/sqrt(n()))


ggplot(plot_data, aes(x = latent_time, y = module_score, colour = clusters_gradients, group = clusters_gradients)) +
  geom_smooth(method="gam", se=FALSE) +
  xlab("Latent time") + ylab("Gene module score") +
  theme_classic()


seurat_data@meta.data$
  
  