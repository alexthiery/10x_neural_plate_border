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
install_github("jokergoo/ComplexHeatmap") # Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. DOI: 10.1093/bioinformatics/btw313
library(ComplexHeatmap)

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
# metadata <- read.csv('./output/NF-downstream_analysis_stacas/transfer_labels/cellrank1/all_stages_filtered_metadata.csv')

seurat_data <- readRDS(list.files(data_path, pattern = "*.RDS", full.names = TRUE)[!list.files(data_path, pattern = "*.RDS") %>% grepl('antler', .)])
# seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/transfer_labels/seurat/rds_files/seurat_label_transfer.RDS')

# load antler data
antler_data <- readRDS(list.files(data_path, pattern = "antler_out.RDS", full.names = TRUE))
# antler_data <- readRDS('./output/NF-downstream_analysis_stacas/transfer_labels/antler/gene_modules/rds_files/antler_out.RDS')

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


debug(GeneModulePheatmap)

scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR', 'eNPB', 'NPB',
                              'aNPB', 'pNPB', 'NC', 'dNC', 'eN', 'eCN', 'NP', 'pNP',
                              'HB', 'iNP', 'MB', 'aNP', 'FB', 'vFB', 'node', 'streak')

seurat_data@meta.data$scHelper_cell_type <- factor(seurat_data@meta.data$scHelper_cell_type, levels = scHelper_cell_type_order)



gms <- antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content
gms_sub <- gms[paste0('GM', c(24, 21, 23, 13, 9, 7, 11))]

# Get plot data from GeneModulePheatmap to plot with Complex Heatmap for extra functionality
plot_data <- GeneModulePheatmap(seurat_obj = seurat_data,  metadata = c('stage', 'scHelper_cell_type'), gene_modules = gms_sub,
                  col_order = c('stage', 'scHelper_cell_type'), col_ann_order = c('stage', 'scHelper_cell_type'), return = 'plot_data')

goi <- which(rownames(plot_data$row_ann) %in% c('EPCAM', 'SALL4',
                                                'HOMER2', 'TFAP2C',
                                                'SIX1', 'EYA2', 'DLX5', 'DLX6', 'GATA2', 'GATA3',
                                                'PAX7', 'SNAI2', 'SOX10',
                                                'LMO1', 'SOX21',
                                                'SIX3', 'PAX6',
                                                'OTX2', 'PAX2', 'WNT4'))

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
  
  temp <- merge(t(temp), seurat_data@meta.data[,c('latent_time', 'lineage_NC_probability', 'lineage_neural_probability', 'lineage_placodal_probability'), drop=FALSE], by=0)
  plot_data <- temp %>%
    column_to_rownames('Row.names') %>%
    pivot_longer(!c(latent_time, lineage_NC_probability, lineage_neural_probability, lineage_placodal_probability)) %>%
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


# set colours
colours <- setNames(colorRampPalette(brewer.pal(9, "Paired"))(length(gms)), names(gms))

# Plot dynamics for each module individually
for(module in names(gms)){
  plot = ggplot(filter(plot_data, module == module), aes(x = latent_time, y = scaled_expression, color = module)) +
    scale_color_manual(values=colours[[module]]) +
    geom_smooth(method="gam", se=FALSE, mapping = aes(weight = value, linetype = lineage, group=lineage)) +
    xlab("Latent time") + ylab("Scaled expression") +
    # facet_wrap(~lineage, dir = 'v') +
    theme_classic()
  
  png(paste0(plot_path, module, '_lineage_dynamics.png'), width = 25, height = 15, units='cm', res=200)
  print(plot)
  graphics.off()
  
  # Calculate average module expression for plotting multi-feature plot
  seurat_data@meta.data[[module]] <-  colMeans(GetAssayData(seurat_data, assay = 'RNA', slot = 'scale.data')[gms[[module]],])
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








PPR_genes <- FindMarkers(seurat_data, ident.1 = c('aPPR', 'pPPR'), ident.2 = c('NC', 'dNC'), group.by = 'scHelper_cell_type', logfc.threshold = 1, only.pos = TRUE)


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

