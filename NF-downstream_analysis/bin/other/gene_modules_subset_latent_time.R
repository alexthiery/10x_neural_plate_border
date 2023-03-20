#!/usr/bin/env Rscript


############################################################################################################################################################
# Functions
############################################################################################################################################################

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
    # re-order row annotations according to hclust output
    row_ann <- row_ann[match(unlist(selected_GM), rownames(row_ann)),,drop=FALSE]
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

ModelMultiLatentExpression <- function(object, goi, latent_col, lineage_col, assay = 'RNA', slot = 'data', cell_sub = colnames(object), padj_method = 'BH'){
  # Get metadata
  metadata <- object@meta.data[cell_sub, c(latent_col, lineage_col)]
  
  # Get expression data
  expression_data <- GetAssayData(object, slot = slot, assay = assay)[goi, rownames(metadata)] %>% as.matrix() %>% t()
  
  # Make data long and run models
  gam_dat <- cbind(metadata[latent_col], metadata[lineage_col], expression_data) %>%
    as.data.frame() %>%
    pivot_longer(!c(!!sym(latent_col), !!sym(lineage_col)), names_to = 'gene', values_to = 'expression') %>%
    group_by(gene) %>%
    do(model = mgcv::gam(expression ~ s(!!sym(latent_col)) + !!sym(lineage_col), data =., family = nb(link='log'), method = 'REML', control = gam.control(maxit = 10000)))
  
  # Get summary output from gams
  gam_dat <- gam_dat %>%
    mutate(slope = summary(model)$p.table[lineage_col, 'Estimate']) %>%
    mutate(p_lineage = summary(model)$p.table[lineage_col, 'Pr(>|z|)']) %>%
    mutate(p_time = as.data.frame(summary(model)$s.table)[['p-value']])
  
  gam_dat$padj_lineage <- p.adjust(gam_dat$p_lineage, padj_method)
  gam_dat$padj_time <- p.adjust(gam_dat$p_time, padj_method)
  
  return(gam_dat)
}

CalcLatentTimeCutoff <- function(latent_time, lineage_probability, top_frac = 0.2, return = 'intercept', verbose=FALSE){
  data = data.frame(latent_time, lineage_probability)
  
  model_data <- data %>%
    filter(lineage_probability > 0 & lineage_probability < 1) %>%
    filter(lineage_probability > quantile(lineage_probability, 1-top_frac)) %>%
    # If cells remaining in top frac have less than 0.5 probability of giving rise to the lineage then remove (required for short lineages)
    filter(lineage_probability > 0.5)
  
  # Fit model to top frac of data
  fit = lm(lineage_probability ~ latent_time, data = model_data)
  
  # Inverse equation to find X for Y = 1
  x <- (1-coef(fit)[1])/coef(fit)[2]
  
  # Identify max latent time value based on cells which have reached lineage probability == 1
  max_latent_time <- data %>% filter(lineage_probability == 1) %>% filter(latent_time == max(latent_time)) %>% dplyr::pull(latent_time)
  
  if(return == 'plot'){
    p = ggplot(data, aes(latent_time, lineage_probability)) + 
      geom_point(size = 0.1) +
      geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2])
    
    return(p)
  }else if(return == 'intercept'){
    if(x > max_latent_time & verbose == TRUE){
      cat('Predicted latent time is later than any cell observed that has reached full lineage absorbtion. Max latent time is set based on oldest cell at lineage_probability == 1')
      return(max_latent_time)
    }else{
      return(unname(x))
    }
  }else{
    stop('return must be one of intercept or plot')
  }
}

PrepareLineageGamData <- function(object, goi, slot = 'data', assay = 'RNA', lineage, latent_col = 'latent_time'){
  
  # Get metadata
  metadata <- object@meta.data[, c('latent_time', lineage)]
  colnames(metadata) <- c(latent_col, 'lineage')
  
  # Calculate latent_time_cutoff
  latent_time_cutoff <- CalcLatentTimeCutoff(latent_time = object@meta.data[[latent_col]], lineage_probability = object@meta.data[[lineage]])
  
  # Get expression data
  expression_data <- GetAssayData(object, slot = slot, assay = assay)[goi, rownames(metadata), drop = FALSE] %>% as.matrix() %>% t()
  
  if(length(goi) > 1){
    expression_data <- data.frame(multi_gene = rowMeans(expression_data))
  }
  
  # Filter cells below max_latent_time
  lineage_expression_data <- cbind(metadata, expression_data) %>%
    as.data.frame() %>%
    filter(!!sym(latent_col) < latent_time_cutoff)
  
  return(lineage_expression_data)
}


LineageGam <- function(lineage_expression_data, goi, latent_col = 'latent_time'){
  return(mgcv::gam(data = lineage_expression_data, formula = as.formula(paste0(goi, " ~ s(", latent_col,", bs = 'cs', k = 5)")),
                   weights = lineage, family = nb(link='log'), method = 'REML', control = gam.control(maxit = 10000)))
}

CalcGamConfidence <- function(gam){
  # Line of 'best fit'
  fit_gam <- predict.gam(gam, type = 'response')
  
  fit_gam_se <- predict.gam(gam, type = 'response', se.fit = TRUE)
  
  # 1.96*se
  fit_gam_ci <- fit_gam_se %>% as.data.frame() %>%
    mutate(upper = fit + (1.96*se.fit)) %>%
    mutate(lower = fit - (1.96*se.fit))
  
  return(fit_gam_ci)
}


RunLineageGamConfidence <- function(object, goi, slot = 'data', assay = 'RNA', lineage, latent_col = 'latent_time'){
  lineage_expression_data <- PrepareLineageGamData(object = object, goi = goi, slot = slot, assay = assay, lineage = lineage, latent_col = latent_col)
  gam <- LineageGam(lineage_expression_data, goi = ifelse(length(goi) > 1, 'multi_gene', goi), latent_col = latent_col)
  gam_ci <- CalcGamConfidence(gam)
  gam_ci <- cbind(gam_ci, lineage_expression_data[,latent_col, drop = FALSE])
  gam_ci <- gam_ci %>% mutate(!!sym(latent_col) := round(!!sym(latent_col), digits = 2))
  return(gam_ci)
}


MultiRunLineageGamConfidence <- function(object, goi, slot = 'data', assay = 'RNA', lineage_cols='auto', latent_col = 'latent_time'){
  if(lineage_cols == 'auto'){
    lineage_cols <- grep('lineage_', colnames(object@meta.data), value = TRUE)
  } else if (!all(lineage_cols %in% colnames(object@meta.data))){
    stop('lineage_cols missing from object metadata - please check metadata column names and rerun')
  }
  
  if(length(goi) > 1){cat('A group of genes has been provided to plot gene dynamics, gene expression will be averaged for the group of genes. \n')}
  
  gams <- list()
  for(lineage in lineage_cols){
    gams[[lineage]] <- RunLineageGamConfidence(object = object, goi = goi, lineage = lineage, slot = slot, assay = assay) %>%
      group_by(!!sym(latent_col)) %>% summarise(upper = max(upper), lower = min(lower), fit = mean(fit))
  }
  return(gams)
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
library(viridis)
set.seed(100)

########################       CELL STATE COLOURS    ########################################
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", "#53A651", "#6D8470",
                                "#87638F", "#A5548D", "#C96555", "#ED761C", "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30",
                                "#CC9F2C", "#AD6428", "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 'MB', 
                                       'vFB', 'aNP', 'node', 'FB', 'pEpi')
########################       STAGE COLOURS     ###########################################
stage_colours = c("#E78AC3", "#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- c("HH4", "HH5", "HH6", "HH7", "ss4", "ss8")
############################################################################################

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
lineage_colours = c('placodal' = '#3F918C', 'NC' = '#DE4D00', 'neural' = '#8000FF')

metadata <- read.csv(list.files(data_path, pattern = "*.csv", full.names = TRUE))
# metadata <- read.csv('./output/NF-downstream_analysis/transfer_labels/cellrank/all_stages_filtered_metadata.csv')

seurat_data <- readRDS(list.files(data_path, pattern = "*.RDS", full.names = TRUE)[!list.files(data_path, pattern = "*.RDS") %>% grepl('antler', .)])
# seurat_data <- readRDS('./output/NF-downstream_analysis/transfer_labels/seurat/rds_files/seurat_label_transfer.RDS')

# load antler data
antler_data <- readRDS(list.files(data_path, pattern = "antler_out.RDS", full.names = TRUE))
# antler_data <- readRDS('./output/NF-downstream_analysis/transfer_labels/antler/gene_modules/rds_files/antler_out.RDS')

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

# Save seurat object with latent time labels added to metadata
saveRDS(seurat_data, paste0(rds_path, "seurat_label_transfer_latent_time.RDS"), compress = FALSE)

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

# Subset gms of interest 
gms <- antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content
ss8 <- subset(x = seurat_data, subset = stage == "ss8")
lineages <- plyr::mapvalues(ss8@meta.data$scHelper_cell_type, c("HB", "MB", "aNP", "FB", "vFB", "aPPR", "pPPR", "dNC"), 
                      c(rep("neural", 5), "PPR", "PPR", "NC"))

ss8 <- AddMetaData(ss8, metadata = lineages, col.name = "lineage")
gms_sub <- DEGeneModules(ss8, gene_modules = gms, active_ident = "lineage", logfc = 0.25)

# Plot unbiasedly selected GMs
plot_data <- GeneModulePheatmap(seurat_obj = seurat_data,  metadata = c('stage', 'scHelper_cell_type'), gene_modules = gms_sub,
                                col_order = c('stage', 'scHelper_cell_type'), col_ann_order = c('stage', 'scHelper_cell_type'), return = 'plot_data')

plot_data$ann_colours$scHelper_cell_type <- scHelper_cell_type_colours[names(plot_data$ann_colours$scHelper_cell_type)]
plot_data$ann_colours$stage <- stage_colours[names(plot_data$ann_colours$stage)]

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
        
        raster_quality = 8
)
graphics.off()

# Plot unbiasedly selected GMs
gms_sub <- gms[c("GM21", "GM24", "GM23", "GM9", "GM13")]

plot_data <- GeneModulePheatmap(seurat_obj = seurat_data,  metadata = c('stage', 'scHelper_cell_type'), gene_modules = gms_sub,
                                col_order = c('stage', 'scHelper_cell_type'), col_ann_order = c('stage', 'scHelper_cell_type'), return = 'plot_data')

plot_data$ann_colours$scHelper_cell_type <- scHelper_cell_type_colours[names(plot_data$ann_colours$scHelper_cell_type)]
plot_data$ann_colours$stage <- stage_colours[names(plot_data$ann_colours$stage)]

goi <- which(rownames(plot_data$row_ann) %in% c('EPCAM', 'SALL4', 'TSPAN13',
                                                'HOMER2', 'TFAP2C', 'BAMBI', 'CITED4',
                                                'SIX1', 'EYA2', 'DLX5', 'DLX6', 'GATA2', 'GATA3', 'TFAP2A',
                                                'PAX7', 'SNAI2', 'SOX10', 'TFAP2B',
                                                'LMO1', 'SOX21'))

png(paste0(plot_path, 'subsetGMs_manual.png'), width = 100, height = 60, res = 800, units = 'cm')
Heatmap(t(plot_data$plot_data), col = PurpleAndYellow(), cluster_columns = FALSE, cluster_rows = FALSE,
        show_column_names = FALSE, column_title = NULL, show_row_names = FALSE, row_title_gp = gpar(fontsize = 45), row_title_rot = 90,
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



################## Calculate maximum latent time values ################## 
png(paste0(plot_path, 'max_NC_latent_time.png'), width = 15, height = 10, res = 200, units = 'cm')
CalcLatentTimeCutoff(latent_time = seurat_data@meta.data[['latent_time']],
                        lineage_probability = seurat_data@meta.data[['lineage_NC_probability']],

                        return = 'plot')
graphics.off()

png(paste0(plot_path, 'max_placodal_latent_time.png'), width = 15, height = 10, res = 200, units = 'cm')
CalcLatentTimeCutoff(latent_time = seurat_data@meta.data[['latent_time']],
                        lineage_probability = seurat_data@meta.data[['lineage_placodal_probability']],
                        return = 'plot')
graphics.off()

png(paste0(plot_path, 'max_neural_latent_time.png'), width = 15, height = 10, res = 200, units = 'cm')
CalcLatentTimeCutoff(latent_time = seurat_data@meta.data[['latent_time']],
                        lineage_probability = seurat_data@meta.data[['lineage_neural_probability']],
                        return = 'plot')
graphics.off()

max_NC <- CalcLatentTimeCutoff(latent_time = seurat_data@meta.data[['latent_time']], lineage_probability = seurat_data@meta.data[['lineage_NC_probability']])
max_placodal <- CalcLatentTimeCutoff(latent_time = seurat_data@meta.data[['latent_time']], lineage_probability = seurat_data@meta.data[['lineage_placodal_probability']])
max_neural <- CalcLatentTimeCutoff(latent_time = seurat_data@meta.data[['latent_time']], lineage_probability = seurat_data@meta.data[['lineage_neural_probability']])

# Iteratively get expression data for each gene module and bind to tidy dataframe
lineage_expression_data <- data.frame()
for(module in names(gms)){
  temp <- GetAssayData(seurat_data, assay = 'RNA', slot = 'scale.data')[gms[[module]],]
  
  temp <- merge(t(temp), seurat_data@meta.data[,c('latent_time', 'lineage_NC_probability', 'lineage_neural_probability', 'lineage_placodal_probability'), drop=FALSE], by=0)
  lineage_expression_data <- temp %>%
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
    bind_rows(lineage_expression_data) %>%
    ungroup()
}


########## Generate a GAM per lineage per gene for plotting ##########


lineage_names = c('lineage_NC_probability' = 'Neural crest', 'lineage_neural_probability' = 'Neural', 'lineage_placodal_probability' = 'Placodal')

# Plot module dynamics for every module
curr_plot_path <- paste0(plot_path, 'gm_lineage_dynamics/')
dir.create(curr_plot_path, recursive = TRUE, showWarnings = FALSE)


for(module in names(gms)){
  gms_plot_data <- MultiRunLineageGamConfidence(seurat_data, goi = gms[[module]], slot = 'data') %>% dplyr::bind_rows(.id = 'id')
  
  plot <-   ggplot(gms_plot_data, aes(latent_time, fit, colour = id, fill = id)) +
    geom_ribbon(aes(ymax = upper, ymin = lower), alpha=0.4, colour = NA) +
    geom_line() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key.size = unit(0.8,"cm"), 
          legend.title = element_blank(),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),  
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.text=element_text(size=20)) +
    ylab(paste0('Average Normalised Expression')) +
    xlab('Latent Time') +
    labs(color = "Lineage", 
         fill = "Lineage") +
    scale_fill_viridis(discrete = TRUE, end = 0.98,
                       labels=lineage_names) +
    scale_color_viridis(discrete = TRUE, end = 0.98,
                        labels=lineage_names)
  
  png(paste0(curr_plot_path, module, '.png'), width = 18, height = 12, res = 200, units = 'cm')
  print(plot)
  graphics.off()
}


# Calculate average module expression for plotting feature plots
curr_plot_path <- paste0(plot_path, 'gm_feature_plots/')
dir.create(curr_plot_path)

for(module in names(gms)){
  seurat_data@meta.data[[module]] <-  colMeans(GetAssayData(seurat_data, assay = 'RNA', slot = 'data')[gms[[module]],])
  
  png(paste0(curr_plot_path, module, '_feature_plot.png'), width = 15, height = 15, units='cm', res=200)
  print(
    FeaturePlot(seurat_data, features = module, pt.size = 1.4) +
      theme_void() +
      theme(plot.title = element_blank(),
            legend.text = element_text(size=16),
            legend.key.size = unit(1, 'cm'))
  ) 
  graphics.off()
}

