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
# metadata <- read.csv('./output/NF-downstream_analysis_stacas/transfer_subset/transfer_ppr_nc_subset/cellrank/transfer_ppr_nc_subset_metadata.csv')

seurat_data <- readRDS(list.files(data_path, pattern = "*.RDS", full.names = TRUE)[!list.files(data_path, pattern = "*.RDS") %>% grepl('antler', .)])
seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/transfer_subset/transfer_ppr_nc_subset/seurat/transfer_cluster/rds_files/transfer_clustered_data.RDS')

# load antler data
antler_data <- readRDS(list.files(data_path, pattern = "antler_out.RDS", full.names = TRUE))
antler_data <- readRDS('./output/NF-downstream_analysis_stacas/transfer_subset/transfer_ppr_nc_subset/antler/transfer_gene_modules/rds_files/antler_out.RDS')

metadata$CellID <- paste0(metadata$CellID, "-1")

# Set RNA to default assay for plotting expression data
DefaultAssay(seurat_data) <- "RNA"

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

# names(gms) <- c(paste(rep('neural', 4), 1:4, sep = '-'), paste(rep('placodal', 6), 1:6, sep = '-'), paste(rep('NC', 3), 1:3, sep = '-'))
# if(is.null(names(gms))){names(gms) = paste0("GM: ", 1:length(gms))}
# names(gms) <- paste0("GM", 1:length(gms))

# Set RNA to default assay for plotting expression data
DefaultAssay(seurat_data) <- "RNA"

# Iteratively get expression data for each gene module and bind to tidy dataframe
plot_data <- data.frame()
for(module in names(gms)){
  temp <- GetAssayData(seurat_data, assay = 'RNA', slot = 'scale.data')[gms[[module]],]
  
  # temp <- merge(t(temp), seurat_data@meta.data[,c('latent_time', 'lineage_NC_probability', 'lineage_neural_probability', 'lineage_placodal_probability'), drop=FALSE], by=0)
  meta <- merge(t(temp), seurat_data@meta.data[,c('latent_time', 'lineage_NC_probability', 'lineage_placodal_probability'), drop=FALSE], by=0)
  plot_data <- meta %>%
    column_to_rownames('Row.names') %>%
    # pivot_longer(!c(latent_time, lineage_NC_probability, lineage_neural_probability, lineage_placodal_probability)) %>%
    pivot_longer(!c(latent_time, lineage_NC_probability, lineage_placodal_probability)) %>%
    rename(scaled_expression = value) %>%
    rename(gene = name) %>%
    pivot_longer(cols = !c(latent_time, gene, scaled_expression)) %>%
    mutate(module = module) %>%
    rename(lineage_probability = value) %>%
    rename(lineage = name) %>%
    group_by(lineage) %>%
    mutate(lineage = unlist(strsplit(lineage, '_'))[2]) %>%
    bind_rows(plot_data) %>%
    ungroup()
}

# Identify the point along latent time when no more cells are identified as terminal cells in each lineage - subsewuently terminate lineage inference at this point

# plot_data <- plot_data %>%
#   group_by(lineage) %>%
#   mutate(max_lineage_probability = max(lineage_probability)) %>%
#   mutate(max_latent_time =max(latent_time[lineage_probability == max_lineage_probability])) %>%
#   filter(latent_time < max_latent_time) %>%
#   ungroup()

# plot_data <- plot_data %>% filter(latent_time < 0.75)


# ggplot(plot_data %>% filter(gene == 'DLX6'), aes(x = latent_time, y = scaled_expression)) +
#   geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, group=lineage, colour = lineage)) +
#   geom_point()
# 
# 
# library(biomaRt)
# # Get biomart GO annotations for TFs
# ensembl = useMart("ensembl",dataset="ggallus_gene_ensembl")
# 
# annotations <- read.csv('./output/NF-downstream_analysis_stacas/seurat_filtering/1_preprocessing/seurat_annotations.csv')
# 
# # Get gene annotation for genes in gms
# gms_annotations <- filter(annotations, Gene %in% unname(unlist(gms)))
# 
# gms_annotations <- getBM(attributes=c("ensembl_gene_id", 'go_id', "name_1006"),
#                          filters = 'ensembl_gene_id',
#                          values = gms_annotations$Accession,
#                          mart = ensembl)
# 
# # Filter genes with TF annotations
# gms_tf <- filter(gms_annotations, go_id %in% c('GO:0003700', 'GO:0043565', 'GO:0000981')) %>% dplyr::pull(ensembl_gene_id) %>% unique()
# 
# # Get gene name for TF from annotation file
# gms_tf <- filter(annotations, Accession %in% gms_tf) %>% dplyr::pull(Gene)
# 
# # Filter gms for TFs
# gms_tf <- lapply(gms, function(x) x[x %in% gms_tf])
# 
# 
# # Filter plot data for identified TFs
# plot_data_tf = filter(plot_data, gene %in% unlist(gms_tf))
# plot_data_tf_split <- split(plot_data_tf, plot_data_tf$module)
# 
# 
# 
# 
# 
# library(gridExtra)
# # Plot gms c(40, 42, 43, 44) across NC lineage (as they are identified as NC modules)
# p1 <- ggplot(filter(plot_data_tf_split$GM1, lineage == 'placodal'), aes(x = latent_time, y = scaled_expression, colour = gene)) +
#   geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, group=gene)) +
#   facet_wrap(~module, dir = 'v', scales = 'free') +
#   theme_classic()
# 
# # p2 <- p1 %+% plot_data_tf_split$GM1
# p3 <- p1 %+% plot_data_tf_split$GM12
# p4 <- p1 %+% plot_data_tf_split$GM14
# p5 <- p1 %+% plot_data_tf_split$GM10
# p6 <- p1 %+% plot_data_tf_split$GM13
# 
# grid.arrange(p1,p3,p4,p5,p6)
# 
# plot_data <- seurat_data@meta.data %>%
#   dplyr::select(c(lineage_NC_probability, lineage_placodal_probability, latent_time)) %>%
#   pivot_longer(cols = !latent_time) %>%
#   filter(name == 'lineage_placodal_probability')
# 
# ggplot(plot_data, aes(x = latent_time, y = value, colour = name)) +
#   geom_point()
# 
# # Plot gms c(10, 11, 12, 13, 14) across placodal lineage
# p1 <- ggplot(filter(plot_data_tf_split$GM10, lineage == 'placodal'), aes(x = latent_time, y = scaled_expression, colour = gene)) +
#   geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, group=gene)) +
#   facet_wrap(~module, dir = 'v', scales = 'free') +
#   theme_classic()
# 
# p2 <- p1 %+% filter(plot_data_tf_split$GM12, lineage == 'placodal')
# p3 <- p1 %+% filter(plot_data_tf_split$GM13, lineage == 'placodal')
# p4 <- p1 %+% filter(plot_data_tf_split$GM14, lineage == 'placodal')
# 
# grid.arrange(p1,p2,p3,p4)
# 
# ggplot(plot_data_tf %>% filter(gene %in% c('IRF6', 'NFKB1', 'PITX2', 'NKX2-5', 'NR2F2')), aes(x = latent_time, y = scaled_expression)) +
#   geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, colour = gene, linetype = lineage)) +
#   facet_wrap(~gene, dir = 'v', scales = 'free') +
#   theme_classic()
# 
# FeaturePlot(seurat_data, 'NFKB1')
# 
# 
# FeaturePlot(seurat_data, c('GATA2', 'GATA3', 'TSPAN13', 'IRF6', 'METRNL'))
# 
# FeaturePlot(seurat_data, c('DLX5', 'SHISA2'))
# 
# FeaturePlot(seurat_data, c('CLDN3', 'CSRP2'))
# 
# FeaturePlot(seurat_data, c('ASS1', 'SIX1', 'BASP1', 'NR2F2'))





# 
# 
# # Plot gms c(40, 42, 43, 44) across NC lineage (as they are identified as NC modules)
# p1 <- ggplot(filter(plot_data_tf_split$GM40, lineage == 'NC'), aes(x = latent_time, y = scaled_expression, colour = gene)) +
#   geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, group=gene)) +
#   facet_wrap(~module, dir = 'v', scales = 'free') +
#   theme_classic()
# 
# # p2 <- p1 %+% plot_data_tf_split$GM1
# p2 <- p1 %+% plot_data_tf_split$GM42
# p3 <- p1 %+% plot_data_tf_split$GM44
# p4 <- p1 %+% plot_data_tf_split$GM43
# 
# grid.arrange(p1,p2,p3,p4)
# # 
# 
# ggplot(plot_data %>% filter(gene %in% c('MSX1', 'GADD45A', 'AGTRAP', 'SOX11')), aes(x = latent_time, y = scaled_expression)) +
#   geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, colour = gene, linetype = lineage)) +
#   facet_wrap(~gene, dir = 'v', scales = 'free') +
#   theme_classic()
# 
# FeaturePlot(seurat_data, c('MSX1', 'GADD45A', 'AGTRAP', 'SOX11'))
# 
# FeaturePlot(seurat_data, c('PAX7', 'ZNF423', 'Z-ENC1'))
# 
# FeaturePlot(seurat_data, c('SOX5', 'OLFML3', 'GLIPR2'))
# 
# FeaturePlot(seurat_data, c('SOX10', 'Z-MEF2C', 'RFTN2', 'ETS1', 'PPP1R1C'))





# Running GAM before ggplot
# some random data
set.seed(100)



library(mgcv)
library(viridis)


# Set RNA to default assay for plotting expression data
DefaultAssay(seurat_data) <- "RNA"

# Iteratively get expression data for each gene module and bind to tidy dataframe
heatmap_data <- data.frame()
for(module in names(gms)){
  temp <- GetAssayData(seurat_data, assay = 'RNA', slot = 'data')[gms[[module]],]
  
  # temp <- merge(t(temp), seurat_data@meta.data[,c('latent_time', 'lineage_NC_probability', 'lineage_neural_probability', 'lineage_placodal_probability'), drop=FALSE], by=0)
  meta <- merge(t(temp), seurat_data@meta.data[,c('latent_time', 'lineage_NC_probability', 'lineage_placodal_probability'), drop=FALSE], by=0)
  heatmap_data <- meta %>%
    column_to_rownames('Row.names') %>%
    # pivot_longer(!c(latent_time, lineage_NC_probability, lineage_neural_probability, lineage_placodal_probability)) %>%
    pivot_longer(!c(latent_time, lineage_NC_probability, lineage_placodal_probability)) %>%
    rename(scaled_expression = value) %>%
    rename(gene = name) %>%
    pivot_longer(cols = !c(latent_time, gene, scaled_expression)) %>%
    mutate(module = module) %>%
    rename(lineage_probability = value) %>%
    rename(lineage = name) %>%
    group_by(lineage) %>%
    mutate(lineage = unlist(strsplit(lineage, '_'))[2]) %>%
    bind_rows(heatmap_data) %>%
    ungroup()
}

multi_gam <- function(data, genes, lineage, max_latent_time = 1){
  
  pdat <- tibble(latent_time = seq(0, max_latent_time, length = 1000))
  
  for(gene in genes){
    temp <- plot_data %>% filter(lineage == !!lineage & gene == !!gene)
    # your model
    mod <- gam(scaled_expression ~ s(latent_time, bs = "cs"), weights = lineage_probability, data = temp)
    # predictions
    pdat <- cbind(pdat, predict.gam(mod, newdata = pdat))
    colnames(pdat)[ncol(pdat)]  <- gene
  }
  # return(pdat %>% pivot_longer(!latent_time))
  rownames(pdat) <- NULL
  return(pdat %>% column_to_rownames('latent_time'))
}

# Model NC
annotations <- list(GM1 = c('MSX1', 'GADD45A', 'AGTRAP', 'SOX11'),
                    GM2 = c('PAX7', 'ZNF423', 'Z-ENC1'),
                    GM3 = c('SOX5', 'OLFML3', 'GLIPR2'),
                    GM4 = c('SOX10', 'Z-MEF2C', 'RFTN2', 'ETS1', 'PPP1R1C'))

annotations <- stack(annotations) %>% column_to_rownames('values')
colnames(annotations) <- c('gm')

test <- multi_gam(heatmap_data, genes = rownames(annotations), lineage = 'NC', max_latent_time = 0.75)

test <- scale(test)
test[test < -2] <- -2
test[test > 2] <- 2

gaps_row <- cumsum(rle(as.vector(annotations[['gm']]))[["lengths"]])

png(paste0(plot_path, 'NC_latent_hm.png'), width = 30, height = 10, res = 200, units = 'cm')
pheatmap(t(test), cluster_cols = FALSE, color = viridis(n=100), border_color = NA, show_colnames = FALSE, treeheight_row = FALSE, cluster_rows = FALSE,
         annotation_row = annotations, gaps_row = gaps_row, cellwidth = 0.6)
graphics.off()


for(mod in unique(annotations$gm)){
  png(paste0(plot_path, 'NC_', mod, '_feature_plot.png'), width = 25, height = 8 * ceiling(((length(rownames(annotations)[annotations$gm == mod])+2)/3)), res = 200, units = 'cm')
  MultiFeaturePlot(seurat_data, rownames(annotations)[annotations$gm == mod], plot_stage = TRUE, plot_clusters = FALSE,
                   plot_celltype = TRUE, celltype_col = 'scHelper_cell_type', n_col = 3, label = paste0('NC ', mod))
  graphics.off()
  
  p1 = ggplot(filter(plot_data, gene %in% rownames(annotations)[annotations$gm == mod] & latent_time < 0.75), aes(x = latent_time, y = scaled_expression, colour = gene)) +
    geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, linetype = lineage)) +
    facet_wrap(~gene, dir = 'v', scales = 'free', ncol = 2) +
    theme_classic()
  
  png(paste0(plot_path, 'NC_', mod, '_dynamics.png'), width = 20, height = 6 * (length(rownames(annotations)[annotations$gm == mod])/2), res = 200, units = 'cm')
  print(p1)
  graphics.off()
}



# Model placodes
annotations <- list(GM1 = c('GATA2', 'GATA3', 'TSPAN13', 'IRF6', 'METRNL'),
                    GM2 = c('DLX5', 'SHISA2'),
                    GM3 = c('CLDN3', 'CSRP2'),
                    GM4 = c('ASS1', 'SIX1', 'BASP1', 'NR2F2'))

annotations <- stack(annotations) %>% column_to_rownames('values')
colnames(annotations) <- c('gm')

test <- multi_gam(heatmap_data, genes = rownames(annotations), lineage = 'placodal', max_latent_time = 0.75)

test <- scale(test)
test[test < -2] <- -2
test[test > 2] <- 2

gaps_row <- cumsum(rle(as.vector(annotations[['gm']]))[["lengths"]])

png(paste0(plot_path, 'placodal_latent_hm.png'), width = 30, height = 10, res = 200, units = 'cm')
pheatmap(t(test), cluster_cols = FALSE, color = viridis(n=100), border_color = NA, show_colnames = FALSE, treeheight_row = FALSE, cluster_rows = FALSE,
         annotation_row = annotations, gaps_row = gaps_row, cellwidth = 0.6)
graphics.off()


for(mod in unique(annotations$gm)){
  png(paste0(plot_path, 'placodal_', mod, '_feature_plot.png'), width = 25, height = 8 * ceiling(((length(rownames(annotations)[annotations$gm == mod])+2)/3)), res = 200, units = 'cm')
  MultiFeaturePlot(seurat_data, rownames(annotations)[annotations$gm == mod], plot_stage = TRUE, plot_clusters = FALSE,
                   plot_celltype = TRUE, celltype_col = 'scHelper_cell_type', n_col = 3, label = paste0('placodal ', mod))
  graphics.off()
  
  
  p1 = ggplot(filter(plot_data, gene %in% rownames(annotations)[annotations$gm == mod] & latent_time < 0.75), aes(x = latent_time, y = scaled_expression, colour = gene)) +
    geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, linetype = lineage)) +
    facet_wrap(~gene, dir = 'v', scales = 'free', ncol = 2) +
    theme_classic()
  
  png(paste0(plot_path, 'placodal_', mod, '_dynamics.png'), width = 20, height = 6 * (length(rownames(annotations)[annotations$gm == mod])/2), res = 200, units = 'cm')
  print(p1)
  graphics.off()
}















# 
# library(gridExtra)
# # Plot gms c(40, 42, 43, 44) across NC lineage (as they are identified as NC modules)
# p1 <- ggplot(filter(plot_data_tf_split$GM1, lineage == 'placodal'), aes(x = latent_time, y = scaled_expression, colour = gene)) +
#   geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, group=gene)) +
#   facet_wrap(~module, dir = 'v', scales = 'free') +
#   theme_classic()
# 
# # p2 <- p1 %+% plot_data_tf_split$GM1
# p3 <- p1 %+% plot_data_tf_split$GM12
# p4 <- p1 %+% plot_data_tf_split$GM14
# p5 <- p1 %+% plot_data_tf_split$GM10
# p6 <- p1 %+% plot_data_tf_split$GM13
# 
# grid.arrange(p1,p3,p4,p5,p6)



# 
# 
# 
# temp <- plot_data %>% filter(lineage == 'placodal' & gene == 'PAX7')
# # your model
# mod <- gam(scaled_expression ~ s(latent_time, bs = "cs"), weights = lineage_probability, data = temp)
# 
# # predictions
# pdat <- tibble(latent_time = seq(0,1, length = 100)) %>% 
#   # New data, this can include any observations you want, or in my case just a sequence for a particular part of the x axis
#   mutate(fit = predict.gam(mod, newdata = .))
# 
# # new one with predicted data - looks weird cos of the weights I made I think
# ggplot(plot_data %>% filter(gene == 'PAX7'), aes(x = latent_time, y = scaled_expression, group)) +
#   # geom_point() +
#   geom_line(data = pdat, aes(y = fit))
# 
# pdat <- pdat %>% column_to_rownames('latent_time')
# 
# pheatmap(t(pdat), cluster_cols = F, cluster_rows = F)








ggplot(plot_data_tf %>% filter(module == 'GM5' & lineage == 'NC'), aes(x = latent_time, y = scaled_expression, colour = gene)) +
  geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, group=gene))

ggplot(plot_data_tf %>% filter(module == 'GM7' & lineage == 'NC'), aes(x = latent_time, y = scaled_expression, colour = gene)) +
  geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, group=gene))



ggplot(plot_data %>% filter(gene == 'DRAXIN' & lineage == 'placodal'), aes(x = latent_time, y = scaled_expression, colour = lineage)) +
  geom_point()+
  geom_smooth(method="gam", se=FALSE, mapping = aes( linetype = lineage, group=lineage))






FeaturePlot(seurat_data, features = c('DRAXIN', 'lineage_placodal_probability', 'latent_time'))

DimPlot(seurat_data, group.by = 'scHelper_cell_type')


hist(plot_data %>% filter(gene == 'PAX7' & latent_time > 0.5) %>% dplyr::pull(scaled_expression))





# set colours
colours <- setNames(colorRampPalette(brewer.pal(9, "Paired"))(length(gms)), names(gms))

# Plot dynamics for each module individually
for(mod in names(gms)){
  plot = ggplot(filter(plot_data, module == mod), aes(x = latent_time, y = scaled_expression, color = module)) +
    scale_color_manual(values=colours[[mod]]) +
    geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, linetype = lineage, group=lineage)) +
    xlab("Latent time") + ylab("Scaled expression") +
    # facet_wrap(~lineage, dir = 'v') +
    theme_classic()
  
  png(paste0(plot_path, mod, '_lineage_dynamics.png'), width = 25, height = 15, units='cm', res=200)
  print(plot)
  graphics.off()
  
  # Calculate average module expression for plotting multi-feature plot
  seurat_data@meta.data[[mod]] <-  colMeans(GetAssayData(seurat_data, assay = 'RNA', slot = 'scale.data')[gms[[mod]],])
}

ncol = ceiling((length(names(gms))+1)/8)+1
nrow = ceiling((length(names(gms))+1)/ncol)

png(paste0(plot_path, 'gene_module_feature.png'), width = ncol*10, height = nrow*10, units = "cm", res = 200)
MultiFeaturePlot(seurat_object = seurat_data, gene_list = grep('GM', colnames(seurat_data@meta.data), value = TRUE), plot_stage = TRUE,
                 stage_col = 'stage', plot_clusters = FALSE, n_col = ncol)
graphics.off()



# Soft assign gene modules to lineages and plot multi-gm dynamics plots
gm_class = plot_data %>%
  group_by(lineage) %>%
  mutate(lineage_expression = scaled_expression*lineage_probability) %>%
  group_by(module, lineage) %>%
  summarise(lineage_expression = mean(lineage_expression), .groups = 'keep') %>%
  filter(lineage_expression > 0) %>%
  dplyr::select(!lineage_expression)


gm_class <- gm_class %>%
  group_by(lineage) %>%
  group_split() %>%
  setNames(unique(gm_class$lineage))

gm_class <- lapply(gm_class, function(x) as.data.frame(x) %>% dplyr::pull(module))

for(lineage in names(gm_class)){
  plot = ggplot(plot_data[plot_data$module %in% gm_class[[lineage]],], aes(x = latent_time, y = scaled_expression, group = module, colour = module)) +
    scale_color_manual(values= colours[gm_class[[lineage]]]) +
    geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, group = module)) +
    xlab("Latent time") + ylab("Scaled expression") +
    facet_wrap(~lineage, dir = 'v') +
    theme_classic()
  
  png(paste0(plot_path, lineage, '_multi_module_dynamics.png'), width = 12, height = 25, units='cm', res=200)
  print(plot)
  graphics.off()
}


# Plot multi-mod dynamics for selected modules (which show pan expression across sub-cell-types)
gm_class = list(
  neural = c('GM16', 'GM11', 'GM12', 'GM13'),
  NC = c('GM8', 'GM2', 'GM7', 'GM1'),
  placodal = c('GM6', 'GM4', 'GM5', 'GM8'))

for(tissue in names(gm_class)){
  plot = ggplot(plot_data %>% filter(module %in% gm_class[[tissue]]) %>% filter(lineage == tissue),
                aes(x = latent_time, y = scaled_expression, group = module, colour = module)) +
    scale_color_manual(values= colours[gm_class[[tissue]]]) +
    geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, group = module)) +
    xlab("Latent time") + ylab("Scaled expression") +
    theme_classic()
  
  png(paste0(plot_path, tissue, '_selected_multi_module_dynamics.png'), width = 15, height = 12, units='cm', res=200)
  print(plot)
  graphics.off()
}








PPR_genes <- FindMarkers(seurat_data, ident.1 = c('aPPR', 'pPPR'), ident.2 = c('NC', 'delaminating_NC'), group.by = 'scHelper_cell_type', logfc.threshold = 1, only.pos = TRUE)


ggplot(plot_data %>% filter(module %in% gm_class$placodal) %>% filter(gene %in% rownames(PPR_genes)) %>% filter(lineage == 'placodal'),
       aes(x = latent_time, y = scaled_expression, group = gene, colour = gene)) +
  # scale_color_manual(values= colours[gm_class[[tissue]]]) +
  geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, group = gene)) +
  xlab("Latent time") + ylab("Scaled expression") +
  theme_classic()


test = c('SOX10', 'TFAP2B', 'PAX7', 'TFAP2A')
ggplot(plot_data %>% filter(module %in% gm_class$NC) %>% filter(gene %in% test) %>% filter(lineage == 'NC'),
       aes(x = latent_time, y = scaled_expression, group = gene, colour = gene)) +
  # scale_color_manual(values= colours[gm_class[[tissue]]]) +
  geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, group = gene)) +
  xlab("Latent time") + ylab("Scaled expression") +
  theme_classic()



test = c('SIX1', 'TFAP2B', 'PAX7', 'TFAP2A', 'EYA2')
ggplot(plot_data %>% filter(gene %in% test) %>% filter(lineage == 'placodal'),
       aes(x = latent_time, y = scaled_expression, group = gene, colour = gene)) +
  # scale_color_manual(values= colours[gm_class[[tissue]]]) +
  geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, group = gene)) +
  xlab("Latent time") + ylab("Scaled expression") +
  theme_classic()

