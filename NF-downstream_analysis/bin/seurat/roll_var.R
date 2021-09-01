library(Seurat)
library(scHelper)
library(tidyverse)

install.packages('RcppRoll')
library(RcppRoll)

seurat_data <- readRDS('~/output/NF-downstream_analysis_stacas/')

metadata <- read.csv('~/output/NF-downstream_analysis_stacas/', row.names = 1)

# edit names based on match between scvelo output and seurat data
rownames(metadata) <- sapply(rownames(metadata), function(x) rownames(seurat_data@meta.data)[grepl(x, rownames(seurat_data@meta.data))] )
# re-order metadata rows based on seurat metadata order
metadata <- metadata[ order(match(rownames(metadata), rownames(seurat_data@meta.data))), ]

# replace seurat metadata with scvelo metadata
seurat_data@meta.data <- metadata

DefaultAssay(seurat_data) <- 'RNA'

data = t(as.matrix(GetAssayData(object = seurat_data)))
# Sort cells based on latent time
latent_time = seurat_data@meta.data[, 'latent_time', drop=FALSE]
latent_time = latent_time[order(latent_time),, drop=FALSE]

data = data[match(rownames(data), rownames(latent_time)),]

# Calculate average rolling variance
rolling_variance <- rowMeans(roll_var(data, n = 100))

# Calculate rolling mean for latent time values
latent_time <- roll_mean(latent_time$latent_time, n = 100)

# Plot average rolling variance across latent time
plot_data <- as.data.frame(cbind(rolling_variance, latent_time))

ggplot(plot_data, aes(x = latent_time, y = rolling_variance)) +
  geom_smooth()