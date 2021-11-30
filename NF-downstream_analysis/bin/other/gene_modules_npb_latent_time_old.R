#!/usr/bin/env Rscript

GeneModulePheatmap <- function (seurat_obj, metadata, col_order = metadata[1], custom_order = NULL, 
                                custom_order_column = NULL, assay = "RNA", slot = "scale.data", 
                                gene_modules, selected_genes = NULL, hide_annotation = NULL, 
                                gaps_row = TRUE, gaps_col = NULL, gm_row_annotation = TRUE, 
                                cell_subset = NULL, use_seurat_colours = TRUE, colour_scheme = c("PRGn", 
                                                                                                 "RdYlBu", "Greys"), col_ann_order = rev(metadata), show_colnames = FALSE, 
                                show_rownames = TRUE, cluster_rows = FALSE, cluster_cols = FALSE, 
                                order_genes = TRUE, annotation_names_row = FALSE, ..., return = 'plot') 
{
  if (!is.null(cell_subset)) {
    seurat_obj <- subset(seurat_obj, cells = cell_subset)
  }
  seurat_obj@meta.data <- seurat_obj@meta.data[, metadata, 
                                               drop = FALSE] %>% mutate_if(is.character, as.factor)
  seurat_obj@meta.data[] <- lapply(seurat_obj@meta.data, function(x) if (is.factor(x)) 
    factor(x)
    else x)
  col_ann <- seurat_obj@meta.data
  col_ann <- col_ann[do.call("order", c(col_ann[col_order], 
                                        list(decreasing = FALSE))), , drop = FALSE]
  if (!is.null(custom_order)) {
    if (is.null(custom_order_column)) {
      "custom_order column must be specified \n"
    }
    if (!setequal(custom_order, unique(col_ann[[custom_order_column]]))) {
      stop("custom_order factors missing from custom_order_column \n\n")
    }
    levels(col_ann[[custom_order_column]]) <- custom_order
    col_ann <- col_ann[order(col_ann[[custom_order_column]]), 
                       , drop = FALSE]
  }
  if (!is.null(gaps_col)) {
    ifelse(class(gaps_col) != "character", stop("gaps_col must be a metadata column name"), 
           gaps_col <- cumsum(rle(as.vector(col_ann[[gaps_col]]))[["lengths"]]))
  }
  if (!is.null(hide_annotation)) {
    col_ann[, hide_annotation] <- NULL
  }
  if (use_seurat_colours == FALSE) {
    ann_colours <- lapply(1:ncol(col_ann), function(x) setNames(colorRampPalette(brewer.pal(9, 
                                                                                            colour_scheme[x])[2:9])(length(unique(col_ann[, x]))), 
                                                                unique(col_ann[, x])))
  }
  else {
    ann_colours <- lapply(colnames(col_ann), function(x) {
      temp <- setNames(ggPlotColours(n = length(levels(col_ann[, 
                                                               x]))), levels(col_ann[, x]))
      temp[match(levels(col_ann[[x]]), names(temp))]
    })
  }
  names(ann_colours) <- colnames(col_ann)
  if (!is.null(selected_genes)) {
    selected_GM <- SubsetGeneModules(gm = gene_modules, selected_genes = selected_genes)
  }
  else {
    if (is.null(names(gene_modules))) {
      names(gene_modules) <- paste0("GM: ", 1:length(gene_modules))
    }
    selected_GM <- gene_modules
  }
  if (gm_row_annotation == TRUE) {
    row_ann <- stack(selected_GM) %>% rename(`Gene Modules` = ind) %>% 
      column_to_rownames("values")
  }
  else {
    row_ann <- NA
  }
  if (gaps_row == TRUE) {
    row_ann <- droplevels(row_ann)
    gaps_row = cumsum(summary(row_ann[["Gene Modules"]], 
                              maxsum = max(lengths(lapply(row_ann, unique)))))
  }
  else {
    gaps_row = NULL
  }
  ann_colours[["Gene Modules"]] <- setNames(colorRampPalette(brewer.pal(9, 
                                                                        "Paired"))(length(unique(row_ann[["Gene Modules"]]))), 
                                            unique(row_ann[["Gene Modules"]]))
  if (order_genes == TRUE) {
    dir.create("scHelper_log/gene_hclustering/", recursive = TRUE, 
               showWarnings = FALSE)
    GMs_ordered_genes <- c()
    for (i in names(selected_GM)) {
      mat <- GetAssayData(object = seurat_obj, assay = assay, 
                          slot = slot)
      dist_mat <- dist(mat[selected_GM[[i]], ], method = "euclidean")
      hclust_avg <- hclust(dist_mat, method = "average")
      ordered_genes <- hclust_avg$labels[c(hclust_avg$order)]
      GMs_ordered_genes[[i]] <- ordered_genes
    }
    selected_GM <- GMs_ordered_genes
  }
  plot_data <- t(as.matrix(x = GetAssayData(object = seurat_obj, 
                                            assay = assay, slot = slot)[unlist(selected_GM), rownames(col_ann), 
                                                                        drop = FALSE]))
  if (!is.null(cell_subset)) {
    cat("rescaling data as cells have been subset \n")
    plot_data <- t(scale(t(plot_data)))
  }
  plot_data <- replace(plot_data, plot_data >= 2, 2)
  plot_data <- replace(plot_data, plot_data <= -2, -2)
  
  if(return == 'plot'){
    return(pheatmap(t(plot_data), color = PurpleAndYellow(), 
                    annotation_col = col_ann[, rev(col_ann_order), drop = FALSE], 
                    annotation_row = row_ann, annotation_colors = ann_colours, 
                    cluster_rows = cluster_rows, cluster_cols = cluster_cols, 
                    show_colnames = show_colnames, show_rownames = show_rownames, 
                    gaps_col = gaps_col, gaps_row = gaps_row, annotation_names_row = annotation_names_row, 
                    ...))
  } else if (return == 'plot_data'){
    return(list(plot_data = plot_data,
                row_ann = row_ann,
                col_ann = col_ann,
                ann_colours = ann_colours))
  } else {
    stop('return must be either "plot" or "plot_data"')
  }
  
}


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
library(devtools)
library(mgcv)
library(ComplexHeatmap) # Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. DOI: 10.1093/bioinformatics/btw313

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


# Params
lineage_colours = c('placodal' = '#DE4D00', 'NC' = '#10E0E8')

#####################################################################################################
#                           Read in data and combine into seurat object                  #
#####################################################################################################

#metadata <- read.csv(list.files(data_path, pattern = "*.csv", full.names = TRUE))
metadata <- read.csv('./output/NF-downstream_analysis_stacas/transfer_subset/transfer_ppr_nc_subset/cellrank/transfer_ppr_nc_subset_metadata.csv')

#seurat_data <- readRDS(list.files(data_path, pattern = "*.RDS", full.names = TRUE)[!list.files(data_path, pattern = "*.RDS") %>% grepl('antler', .)])
seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/transfer_subset/transfer_ppr_nc_subset/seurat/transfer_cluster/rds_files/transfer_clustered_data.RDS')

# load antler data
#antler_data <- readRDS(list.files(data_path, pattern = "antler_out.RDS", full.names = TRUE))
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


#############################################################################
#                           NC Plac Neural gm subset                        #
#############################################################################

####### Select modules of interest and plot new heatmap #######


# Set RNA to default assay for plotting expression data
DefaultAssay(seurat_data) <- "RNA"


scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR', 'eNPB', 'NPB',
                              'aNPB', 'pNPB', 'NC', 'dNC', 'eN', 'eCN', 'NP', 'pNP',
                              'HB', 'iNP', 'MB', 'aNP', 'FB', 'vFB', 'node', 'streak')

seurat_data@meta.data$scHelper_cell_type <- factor(seurat_data@meta.data$scHelper_cell_type, levels = scHelper_cell_type_order)

gms <- antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content

PPR_gms <- c("GM12", "GM14", "GM13")
NC_gms <- c("GM40", "GM42", "GM44", "GM43")

gms_sub <- gms[unlist(c(PPR_gms, NC_gms))]

# Get plot data from GeneModulePheatmap to plot with Complex Heatmap for extra functionality
plot_data <- GeneModulePheatmap(seurat_obj = seurat_data,  metadata = c('stage', 'scHelper_cell_type'), gene_modules = gms_sub,
                                col_order = c('stage', 'scHelper_cell_type'), col_ann_order = c('stage', 'scHelper_cell_type'), return = 'plot_data')

goi <- which(rownames(plot_data$row_ann) %in% c('FAM184B', 'TSPAN13', 'GATA2', 'GATA3',
                                                'DLX5', 'DLX6', 'SIX3', 'PAX6', 'NFKB1', 
                                                'SIX1', 'EYA2', 'ASS1',
                                                'MSX1', 'Z-ENC1',
                                                'PAX7', 'SNAI2',
                                                'OLFML3', 'WNT1',
                                                'SOX10', 'RFTN2'
                                                ))

png(paste0(plot_path, 'subsetGMs.png'), width = 100, height = 60, res = 800, units = 'cm')
Heatmap(t(plot_data$plot_data), col = PurpleAndYellow(), cluster_columns = FALSE, cluster_rows = FALSE,
        show_column_names = FALSE, column_title = NULL, show_row_names = FALSE, row_title_gp = gpar(fontsize = 45), row_title_rot = 0,
        row_split = plot_data$row_ann$`Gene Modules`, column_split = plot_data$col_ann$stage, 
        heatmap_legend_param = list(
          title = "Scaled expression", at = c(-2, 0, 2), 
          labels = c(-2, 0, 2),
          legend_height = unit(11, "cm"),
          grid_width = unit(1.5, "cm"),
          title_gp = gpar(fontsize = 35, fontface = 'bold'),
          labels_gp = gpar(fontsize = 30, fontface = 'bold'),
          title_position = 'lefttop-rot',
          x = unit(13, "npc")
        ),
        
        top_annotation = HeatmapAnnotation(stage = anno_block(gp = gpar(fill = plot_data$ann_colours$stage),
                                                              labels = levels(plot_data$col_ann$stage),
                                                              labels_gp = gpar(col = "white", fontsize = 50, fontface='bold'))),
        
        bottom_annotation = HeatmapAnnotation(scHelper_cell_type = anno_simple(x = as.character(plot_data$col_ann$scHelper_cell_type),
                                                                               col = plot_data$ann_colours$scHelper_cell_type, height = unit(1, "cm")), show_annotation_name = FALSE,
                                              labels = anno_mark(at = cumsum(rle(as.character(plot_data$col_ann$scHelper_cell_type))$lengths) - floor(rle(as.character(plot_data$col_ann$scHelper_cell_type))$lengths/2),
                                                                 labels = rle(as.character(plot_data$col_ann$scHelper_cell_type))$values,
                                                                 which = "column", side = 'bottom',
                                                                 labels_gp = gpar(fontsize = 40), lines_gp = gpar(lwd=8))),
        
        right_annotation = rowAnnotation(foo = anno_mark(at = goi,
                                                         labels = rownames(plot_data$row_ann)[goi],
                                                         labels_gp = gpar(fontsize = 35), lines_gp = gpar(lwd=8))),
        
        raster_quality = 8
)
graphics.off()



# Iteratively get expression data for each gene module and bind to tidy dataframe
plot_data <- data.frame()
for(module in names(gms)){
  temp <- GetAssayData(seurat_data, assay = 'RNA', slot = 'scale.data')[gms[[module]],]
  
  # temp <- t(scale(t(temp)))
  
  temp <- merge(t(temp), seurat_data@meta.data[,c('latent_time', 'lineage_NC_probability', 'lineage_placodal_probability'), drop=FALSE], by=0)
  plot_data <- temp %>%
    column_to_rownames('Row.names') %>%
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


# Generate gams for each group (lineage, module)
gams <- plot_data %>%
  group_by(module, lineage) %>%
  do(gams = gam(scaled_expression ~ s(latent_time, bs = "cs"), weights = lineage_probability, data = .))

# Generate latent time values to predict gams on
pdat <- tibble(latent_time = seq(0, 0.75, length = 1000))

# Generate predicted values in long format
plot_data <- data.frame()
for(row in 1:nrow(gams)){
  plot_data <- rbind(plot_data, data.frame(module = gams[[row, 'module']],
                                           lineage = gams[[row, 'lineage']],
                                           scaled_expression = predict.gam(gams[[row,'gams']][[1]], newdata = pdat),
                                           pdat))
}

for(module in unique(plot_data$module)){
  plot <- ggplot(filter(plot_data, module == !!module), aes(x = latent_time, y = scaled_expression, colour = lineage)) +
    geom_point(size = 0.25) +
    geom_line() +
    scale_colour_manual(values=lineage_colours) +
    theme_classic()
  
  png(paste0(plot_path, module, '_lineage_dynamics.png'), width = 15, height = 10, units='cm', res=200)
  print(plot)
  graphics.off()
}

# Calculate average module expression for plotting feature plots
for(module in names(gms)){
  seurat_data@meta.data[[module]] <-  colMeans(GetAssayData(seurat_data, assay = 'RNA', slot = 'scale.data')[gms[[module]],])
  
  png(paste0(plot_path, module, '_feature_plot.png'), width = 15, height = 15, units='cm', res=200)
  print(FeaturePlot(seurat_data, features = module, pt.size = 1) +
          theme_void() +
          theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))) 
  graphics.off()
}






#####################################################################################################
#                    plot latent time heatmaps using GAM data for GMs of interest                  #
#####################################################################################################


# Are they properly row scaled?
# Gene modules need to be arranged by specific order
# Latent time value not derived from plot data -> max(latent_time) is not always 1!

for(module in unique(plot_data$module)){
  plot <- ggplot(filter(plot_data, module == !!module), aes(x = latent_time, y = scaled_expression, colour = lineage)) +
    geom_point(size = 0.25) +
    geom_line() +
    scale_colour_manual(values=lineage_colours) +
    theme_classic()
  
  png(paste0(plot_path, module, '_lineage_dynamics.png'), width = 15, height = 10, units='cm', res=200)
  print(plot)
  graphics.off()
}

# Filter NC modules for plotting heatmap
NC_latent_time <- filter(plot_data, lineage == 'NC') %>%
  filter(module %in% NC_gms) %>%
  pivot_wider(values_from = scaled_expression, names_from = module) %>%
  dplyr::select(!c(lineage, latent_time))

# NC_latent_time = t(scale(t(NC_latent_time)))

library(circlize)
col_fun = colorRamp2(c(0, 1000), c("purple", "yellow"))
column_ha = HeatmapAnnotation(labels = anno_mark(at = c(2, 500, 999), labels = c("0", "Latent Time", "1"), which = "column", side = "top", 
                                                 labels_gp = gpar(fontsize = 22), 
                                                 labels_rot = 0, link_gp = gpar(lty = 0),
                                                 link_height = unit(0.5, "mm")),
                              Latent_time = 1:(ncol(t(NC_latent_time))), col = list(Latent_time = col_fun), show_legend = FALSE, 
                              annotation_label = "Latent time", show_annotation_name = FALSE)

library(viridis)
Heatmap(t(NC_latent_time), cluster_rows = FALSE, cluster_columns = FALSE,
        show_column_names = FALSE, show_row_names = FALSE,
        col = viridis(n=100), row_split = colnames(NC_latent_time),
        row_title_gp = gpar(fontsize = 22),
        top_annotation = column_ha,
        heatmap_legend_param = list(
          title = "Scaled Expression",
          title_gp = gpar(fontsize = 18),
          legend_height = unit(8, "cm"),
          grid_width = unit(1, "cm"),
          title_position = "leftcenter-rot",
          labels_gp = gpar(fontsize = 16)),
        raster_quality = 4)


##################    PPR GMs on PPR lineage
# Filter NC modules for plotting heatmap
PPR_latent_time <- filter(plot_data, lineage == 'placodal') %>%
  filter(module %in% PPR_gms) %>%
  pivot_wider(values_from = scaled_expression, names_from = module) %>%
  dplyr::select(!c(lineage, latent_time))

col_fun = colorRamp2(c(0, 1000), c("purple", "yellow"))
column_ha = HeatmapAnnotation(labels = anno_mark(at = c(2, 500, 999), labels = c("0", "Latent Time", "1"), which = "column", side = "top", 
                                                 labels_gp = gpar(fontsize = 22), 
                                                 labels_rot = 0, link_gp = gpar(lty = 0),
                                                 link_height = unit(0.5, "mm")),
                              Latent_time = 1:(ncol(t(PPR_latent_time))), col = list(Latent_time = col_fun), show_legend = FALSE, 
                              annotation_label = "Latent time", show_annotation_name = FALSE)

# PPR_latent_time = scale(PPR_latent_time)

Heatmap(t(PPR_latent_time), cluster_rows = FALSE, cluster_columns = FALSE,
        show_column_names = FALSE, show_row_names = FALSE,
        col = viridis(n=100), row_split = colnames(PPR_latent_time),
        row_title_gp = gpar(fontsize = 22),
        top_annotation = column_ha,
        heatmap_legend_param = list(
          title = "Scaled Expression",
          title_gp = gpar(fontsize = 18),
          legend_height = unit(8, "cm"),
          grid_width = unit(1, "cm"),
          title_position = "leftcenter-rot",
          labels_gp = gpar(fontsize = 16)),
        raster_quality = 4)









#####################################################################################################
#                plot latent time line graphs using GAM data for genes of interest                  #
#####################################################################################################

##################    NC genes on NC lineage
annotations <- list(GM1 = c('MSX1', 'GADD45A', 'AGTRAP', 'SOX11'),
                    GM2 = c('PAX7', 'ZNF423', 'Z-ENC1'),
                    GM3 = c('OLFML3'),
                    GM4 = c('SOX10', 'Z-MEF2C', 'RFTN2', 'ETS1', 'PPP1R1C'))

annotations <- list(GM1 = c('MSX1', 'GADD45A', 'AGTRAP', 'SOX11'))

annotations <- stack(annotations) %>% column_to_rownames('values')
colnames(annotations) <- c('gm')

NC_genes_NC_lineage <- multi_gam(expression_data, variable = "gene", values = rownames(annotations), lineage = 'NC', max_latent_time = 0.75)

gaps_row <- cumsum(rle(as.vector(annotations[['gm']]))[["lengths"]])

png(paste0(plot_path, 'NC_genes_latent_hm.png'), width = 30, height = 10, res = 200, units = 'cm')
pheatmap(t(NC_genes_NC_lineage), cluster_cols = FALSE, color = viridis(n=100), border_color = NA, show_colnames = FALSE, treeheight_row = FALSE, cluster_rows = FALSE,
         annotation_row = annotations, gaps_row = gaps_row, cellwidth = 0.6)
graphics.off()


NC_genes_PPR_lineage <- multi_gam(expression_data, variable = "gene", values = rownames(annotations), lineage = 'placodal', max_latent_time = 0.75)


NC_genes_NC_lineage <- NC_genes_NC_lineage %>% as.data.frame %>% mutate(lineage = "NC") %>% rownames_to_column("latent_time")
NC_genes_PPR_lineage <- NC_genes_PPR_lineage %>% as.data.frame %>% mutate(lineage = "placodal") %>% rownames_to_column("latent_time")

NC_genes_both_lineages <- rbind(NC_genes_NC_lineage, NC_genes_PPR_lineage)
NC_genes_both_lineages <- NC_genes_both_lineages %>% pivot_longer(cols = !c(latent_time, lineage)) %>% rename(gene = name)


#####     NEED TO WORK HERE TO SET MAX LATENT TIME??
png(paste0(plot_path, 'NC_genes_latent_line.png'), width = 30, height = 20, res = 200, units = 'cm')
ggplot(NC_genes_both_lineages, aes(x = latent_time, y = value, colour = lineage)) +
  geom_point(size = 1) +
  geom_line() +
  facet_wrap(~gene, dir = 'v', scales = 'free', ncol = 2) +
  theme_classic()
graphics.off()

# make sure latent time read in as continuous variable
# change the colours 

ggplot(filter(expression_data, gene %in% rownames(annotations) & latent_time < 0.75), aes(x = latent_time, y = scaled_expression, colour = gene)) +
  geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, linetype = lineage)) +
  facet_wrap(~gene, dir = 'v', scales = 'free', ncol = 2) +
  theme_classic()








########################################################       OLD        ###################################################
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


ggplot(plot_data %>% filter(gene == 'DLX6'), aes(x = latent_time, y = scaled_expression)) +
  geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, group=lineage, colour = lineage)) +
  geom_point()
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

test <- multi_gam(expression_data, variable = "gene", values = rownames(annotations), lineage = 'NC', max_latent_time = 0.75)

test <- scale(test)
test[test < -2] <- -2
test[test > 2] <- 2

gaps_row <- cumsum(rle(as.vector(annotations[['gm']]))[["lengths"]])

png(paste0(plot_path, 'NC_latent_hm.png'), width = 30, height = 10, res = 200, units = 'cm')
pheatmap(t(test), cluster_cols = FALSE, color = viridis(n=100), border_color = NA, show_colnames = FALSE, treeheight_row = FALSE, cluster_rows = FALSE,
         annotation_row = annotations, gaps_row = gaps_row, cellwidth = 0.6)
graphics.off()

pheatmap(t(test), cluster_cols = FALSE, color = viridis(n=100), border_color = NA, show_colnames = FALSE, treeheight_row = FALSE, cluster_rows = FALSE,
         cellwidth = 0.6)



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

