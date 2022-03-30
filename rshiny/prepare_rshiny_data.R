# Prepare skinny data for shiny
library(Seurat)
library(tidyverse)

# Get RDS objects for all different data subsets - Full data and NPB subset are taken from cellrank output to enable running dynamics plots
data_paths = data.frame(paths = c(
  './output/NF-downstream_analysis_stacas/transfer_labels/gene_modules_subset_latent_time/rds_files/seurat_label_transfer_latent_time.RDS',
  './output/NF-downstream_analysis_stacas/transfer_subset/transfer_ppr_nc_subset/gene_modules_npb_latent_time/rds_files/seurat_npb_subset_latent_time.RDS',
  list.files('./output/NF-downstream_analysis_stacas/stage_split/', recursive = TRUE, pattern = 'cell_state_classification.RDS', full.names = TRUE)
),
row.names = c('Full data', 'NPB subset', 'HH4', 'HH5', 'HH6', 'HH7', 'ss4', 'ss8'))




seurat_list <- lapply(data_paths$paths, function(x) readRDS(x))
names(seurat_list) <- row.names(data_paths)

# Slim seurat object to make it run faster
seurat_list <- lapply(seurat_list, function(x) DietSeurat(x, counts = FALSE, data = TRUE, scale.data = FALSE, dimreducs = 'umap', assays = 'RNA'))

# Specify columns to keep
subset_cols <- c('stage', 'scHelper_cell_type', 'seurat_clusters', 'latent_time', 'terminal_states', 'lineage_NC_probability', 'lineage_neural_probability', 'lineage_placodal_probability')

# rename metadata columns and drop other cols
seurat_list <- lapply(seurat_list, function(x) {x@meta.data <- x@meta.data %>%
  dplyr::select(any_of(subset_cols)) %>%
  rename(
    'Stage' = stage,
    'Cell state' = scHelper_cell_type,
    'Clusters' = seurat_clusters
  )
return(x)
})

# filter for variable features
seurat_list_var_features <- lapply(seurat_list, VariableFeatures) %>% Reduce(c, .) %>% unique()

seurat_list <- lapply(seurat_list, function(x) subset(x, features = seurat_list_var_features))

saveRDS(seurat_list, './output/rshiny_input.RDS', compress = FALSE)
