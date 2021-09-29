library(Seurat)
library(scHelper)
library(tidyverse)

install.packages('RcppRoll')
library(RcppRoll)

seurat_data <- readRDS('~/output/NF-downstream_analysis_stacas/filtered_seurat/seurat/state_classification/rds_files/contamination_cell_state_classification.RDS')

metadata <- read.csv('~/output/NF-downstream_analysis_stacas/filtered_seurat/cellrank/NF-scRNAseq_alignment_out_metadata.csv', row.names = 1)

antler_data <- readRDS('~/output/NF-downstream_analysis_stacas/stage_split/ss8_splitstage_data/antler/stage_gene_modules/rds_files/antler_out.RDS')

# edit names based on match between scvelo output and seurat data
rownames(metadata) <- sapply(rownames(metadata), function(x) rownames(seurat_data@meta.data)[grepl(x, rownames(seurat_data@meta.data))] )
# re-order metadata rows based on seurat metadata order
metadata <- metadata[ order(match(rownames(metadata), rownames(seurat_data@meta.data))), ]

# replace seurat metadata with scvelo metadata
seurat_data@meta.data <- metadata

DefaultAssay(seurat_data) <- 'RNA'

data = t(as.matrix(GetAssayData(object = seurat_data)))
# Sort cells based on latent time
latent_time_df = seurat_data@meta.data[, 'latent_time', drop=FALSE]
latent_time_df = latent_time_df[order(latent_time_df$latent_time),, drop=FALSE]

data = data[match(rownames(data), rownames(latent_time_df)),]

# subset data for GM of interest
gm_data = data[,antler_data$gene_modules$lists$unbiasedGMs_DE$content$GM2]

window_size = 200
ticker = 1
roll_cor = c()
latent_time = c()
for(i in window_size : nrow(gm_data)){
  cor <- cor(gm_data[ticker:i,])
  cor <- mean(cor[lower.tri(cor)])
  roll_cor <- c(roll_cor, cor)
  
  latent_time <- c(latent_time, mean(latent_time_df$latent_time[ticker:i]))
  ticker = ticker+1
}

plot_data = data.frame(roll_cor = roll_cor, latent_time = latent_time)

ggplot(plot_data, aes(x = latent_time, y = roll_cor)) + 
  geom_point()+
  geom_smooth(method = 'gam', se = FALSE)





# 
# 
# # install.packages('roll')
# 
# # library(roll)
roll_cor_data = roll_cor(gm_data, width=250)

av_roll_cor_data = c()
for(i in 1:dim(roll_cor_data)[3]){
  mat_slice = roll_cor_data[,,i]
  av_roll_cor_data = c(av_roll_cor_data, mean(mat_slice[lower.tri(mat_slice)]))
}

temp = data.frame(latent_time = latent_time_df, roll_cor = av_roll_cor_data)

ggplot(temp, aes(x = latent_time, y = roll_cor)) +
  geom_point()+
  geom_smooth(method = 'gam', se = FALSE)







temp[,,1] = temp[lower.tri(temp[,,1])]

# Calculate average rolling variance
rolling_variance <- rowMeans(roll_var(data, n = 100))

# Calculate rolling mean for latent time values
latent_time <- roll_mean(latent_time$latent_time, n = 100)

# Plot average rolling variance across latent time
plot_data <- as.data.frame(cbind(rolling_variance, latent_time))

ggplot(plot_data, aes(x = latent_time, y = rolling_variance)) +
  geom_smooth()