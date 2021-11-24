#!/usr/bin/env Rscript

# Define arguments for Rscript
library(optparse)
library(future)
library(Seurat)
library(pheatmap)
library(tidyverse)
library(Antler)
library(RColorBrewer)
library(scHelper)

### TEMP - overwrite scHelper GeneModuleOrder function
GeneModuleOrder <- function (seurat_obj, gene_modules, metadata_1 = NULL, order_1 = NULL, 
          metadata_2 = NULL, order_2 = NULL, rename_modules = NULL, 
          plot_path = "scHelper_log/GM_classification/") 
{
  classified_gms_1 <- GeneModuleClassify(seurat_obj, gene_modules, 
                                         metadata = metadata_1, plot_path = plot_path)
  classified_gms_1 <- classified_gms_1 %>% arrange(match(!!sym(metadata_1), 
                                                         order_1)) %>% group_by(!!sym(metadata_1)) %>% mutate(pos = 1:n()) %>% 
    mutate(new_name = paste(!!sym(metadata_1), pos, sep = "-")) %>% 
    dplyr::select(gene_module, !!sym(metadata_1), new_name)
  ordered_gms <- gene_modules[order(match(names(gene_modules), 
                                          classified_gms_1$gene_module))] # ordered on metadata_1 but not renamed
  if (is.null(metadata_2) || is.na(metadata_2) || is.nan(metadata_2)) {
    print(paste0("Gene modules ordered only on ", metadata_1))
  }
  else {
    print(paste0("Gene modules ordered on ", metadata_1, 
                 " AND ", metadata_2))
    temp_seurat <- SplitObject(seurat_data, split.by = metadata_1)
    classified_gms_2 <- c()
    for (i in order_1) {
      print(i)
      subset_gms <- gene_modules[classified_gms_1 %>% filter(!!sym(metadata_1) == 
                                                               i) %>% dplyr::pull(gene_module)]
      if (length(subset_gms) != 0) {
        temp <- GeneModuleClassify(seurat_data, subset_gms, 
                                   metadata = metadata_2, plot_path = paste0(plot_path, 
                                                                             i, "/"))
        classified_gms_2 <- rbind(classified_gms_2, temp)
      }
    }
    
    classified_gms_2 <- classified_gms_2 %>% add_column(!!sym(metadata_1) := classified_gms_1[[metadata_1]])
    
    classified_gms_2 <- classified_gms_2 %>% arrange(!!sym(metadata_1), (match(!!sym(metadata_2), order_2))) %>% 
      mutate(pos = 1:n()) %>% mutate(new_name = paste(!!sym(metadata_2), pos, sep = "-")) %>% 
      dplyr::select(gene_module, !!sym(metadata_2), !!sym(metadata_1), new_name)
    
    ordered_gms <- gene_modules[order(match(names(gene_modules), classified_gms_2$gene_module))]
    
    if (!is.null(rename_modules) && rename_modules == metadata_2) {
      names(ordered_gms) <- classified_gms_2$new_name
    }
  }
  if (!is.null(rename_modules) && rename_modules == metadata_1) {
    names(ordered_gms) <- classified_gms_1$new_name
  }
  return(ordered_gms)
}

### TEMP - overwrite pheatmap function to return dfs to pass to complex heatmap
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
                    gaps_col = gaps_col, gaps_row = gaps_row, annotation_names_row = annotation_names_row))
  } else if (return == 'plot_data'){
    return(list(plot_data = plot_data,
                row_ann = row_ann,
                col_ann = col_ann,
                ann_colours = ann_colours))
  } else {
    stop('return must be either "plot" or "plot_data"')
  }
  
}



# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-m", "--meta_col"), action = "store", type = "character", help = "Name of metadata column containing grouping information", default = 'scHelper_cell_type'),
  make_option(c("-f", "--force_order"), action = "store", type = "character", help = "Comma separated values specifying metadata columns used for GeneModuleOrdering", default = NULL),
  make_option(c("", "--verbose"), action = "store_true", type = "logical", help = "Verbose", default = FALSE)
)

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    plot_path = "./output/NF-downstream_analysis_stacas/run_split/2_splitrun_data/antler/plots/"
    rds_path = "./output/NF-downstream_analysis_stacas/run_split/2_splitrun_data/antler/rds_files/"
    gm_path = "./output/NF-downstream_analysis_stacas/run_split/2_splitrun_data/antler/gene_module_lists/"
    antler_path = "./output/NF-downstream_analysis_stacas/run_split/2_splitrun_data/antler/antler_data/"
    data_path = "./output/NF-downstream_analysis_stacas/run_split/2_splitrun_data/seurat/run_state_classification/rds_files/"
    
    ncores = 8
    meta_col = 'scHelper_cell_type'
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    gm_path = "./gene_module_lists/"
    antler_path = "./antler_data/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    meta_col = opt$meta_col
    
    # Multi-core when running from command line
    plan("multiprocess", workers = ncores)
    options(future.globals.maxSize = 16* 1024^3) # 32gb
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
  dir.create(gm_path, recursive = T)
  dir.create(antler_path, recursive = T)
}
seurat_data <- readRDS(list.files(data_path, full.names = TRUE))

########################################################################################################################################################
##################################################   Calculating gene modules using Antler:   ##########################################################

# switch to RNA assay for viewing expression data
DefaultAssay(seurat_data) <- "RNA"

seurat_data <- DietSeurat(seurat_data, counts = TRUE, data = TRUE, scale.data = TRUE, assays = 'RNA')
# seurat_data <- subset(seurat_data, cells = colnames(seurat_data)[1:2500])

# strip end of cell names as this is incorrectly reformated in Antler
seurat_data <- RenameCells(seurat_data, new.names = gsub('-', '_', colnames(seurat_data)))

antler_data <- data.frame(row.names = colnames(seurat_data),
                          "timepoint" = as.numeric(substring(colnames(seurat_data), 3, 3)),
                          "treatment" = rep("null", ncol(seurat_data)),
                          "replicate_id" = rep(1, ncol(seurat_data))
)

# save pheno data
write.table(antler_data, file = paste0(antler_path, "phenoData.csv"), row.names = T, sep = "\t", col.names = T)

# save count data
write.table(GetAssayData(seurat_data, assay = "RNA", slot = "counts"), file = paste0(antler_path, "assayData.csv"), row.names = T, sep = "\t", col.names = T, quote = F)

########################################################################################################
#                            Load Antler data and generate correlation matrix                          #
########################################################################################################

# Change plot path
antler_data <- Antler$new(output_folder = plot_path, num_cores = ncores)
antler_data$load_dataset(folder_path = antler_path)
antler_data$exclude_unexpressed_genes(min_cells=10, min_level=1, verbose=T, data_status='Raw')

antler_data$normalize(method = 'MR')

###########################
# Calculate GMs unbiasedly
antler_data$gene_modules$identify(
  name                  = "unbiasedGMs",
  corr_t                = 0.3,  # the Spearman correlation treshold
  corr_min              = 3,    # min. number of genes a gene must correlate with
  mod_min_cell          = 10,   # min. number of cells expressing the module
  mod_consistency_thres = 0.4,  # ratio of expressed genes among "positive" cells
  process_plots         = TRUE)


if (is.null(names(antler_data$gene_modules$lists$unbiasedGMs$content))) {
  names(antler_data$gene_modules$lists$unbiasedGMs$content) <- paste0("GM", 1:length(antler_data$gene_modules$lists$unbiasedGMs$content))
}

########## DE GMs ##############
# Plot gene modules with at least 50% of genes DE > 0.25 logFC & FDR < 0.001
gms <- DEGeneModules(seurat_data, antler_data$gene_modules$get("unbiasedGMs"), logfc = 0.5, pval = 0.001, selected_gene_proportion = 0.5, active_ident = meta_col)

# save unbiasedGMs_DE in antler object
antler_data$gene_modules$set(name= "unbiasedGMs_DE", content = gms)

########## DE batch filter GMs ##############
# Filter gene modules which are deferentially expressed across batches - first filter stages which have multiple runs to test for DEA
if(length(unique(seurat_data$run)) > 1){
  temp_seurat <- subset(seurat_data, cells = seurat_data@meta.data %>% filter(stage %in% c('HH6', 'ss4')) %>% rownames())
    
  batch_gms <- DEGeneModules(temp_seurat, antler_data$gene_modules$get("unbiasedGMs"), logfc = 0.25, pval = 0.001, selected_gene_proportion = 0.5, active_ident = 'run')
  gms <- antler_data$gene_modules$lists$unbiasedGMs_DE$content[!names(antler_data$gene_modules$lists$unbiasedGMs_DE$content) %in% names(batch_gms)]
  
  # save unbiasedGMs_DE in antler object
  antler_data$gene_modules$set(name= "unbiasedGMs_DE_batchfilt", content = gms)
}


########################################################################################################################################################
##################################################   Setting metadata and colours for heatmaps:   ################################################

# Set metadata and order cells in heatmap
metadata <- if(length(unique(seurat_data@meta.data$stage)) == 1){meta_col
}else{c("stage", meta_col)}
metadata <- if(length(unique(seurat_data@meta.data$run)) == 1){metadata
}else{c(metadata, "run")}

# Allow manual setting of metadata using --force_order command line arg
if(!is.null(opt$force_order)){metadata <- unlist(str_split(opt$force_order, pattern = ','))}


## ADD STAGE AND SCHELPER COLOURS HERE AND AT TOP OF SCRIPT


########################################################################################################################################################
##################################################   Plotting heatmaps with ordered GMs:   #############################################################

# Extract ordering of gms from metadata
labels <- c("stage", "scHelper_cell_type", "seurat_clusters")
stage_order <- levels(seurat_data@meta.data$stage)
scHelper_cell_type_order <- levels(seurat_data@meta.data$scHelper_cell_type)

metadata_1 <- NULL
order_1 <- NULL
metadata_2 <- NULL
order_2 <- NULL

if(sum(labels %in% metadata) !=0){
  metadata_1 <- labels[labels %in% metadata][1]
  order_1 <- get(paste0(metadata_1, '_order'))
  
  if(sum(labels %in% metadata) == 2){
    metadata_2 <- labels[labels %in% metadata][2]
    order_2 <- get(paste0(metadata_2, '_order'))
  }
} else {
  print(paste(c(labels, 'not found in metadata. GMs will not be ordered'), collapse = ' '))
}

##########################################################################################
# plot all gene modules (unbiasedGMs)

# Order gms
if (!is.null(metadata_1)){
  antler_data$gene_modules$lists$unbiasedGMs$content <- GeneModuleOrder(seurat_obj = seurat_data, gene_modules = antler_data$gene_modules$lists$unbiasedGMs$content,
                                                                   metadata_1 = metadata_1, order_1 = order_1,
                                                                   metadata_2 = metadata_2, order_2 = order_2,
                                                                   plot_path = "scHelper_log/GM_classification/unbiasedGMs/")
}

png(paste0(plot_path, 'unbiasedGMs.png'), height = 150, width = 75, units = 'cm', res = 400)
GeneModulePheatmap(seurat_obj = seurat_data, metadata = metadata, gene_modules = antler_data$gene_modules$lists$unbiasedGMs$content,
                   show_rownames = FALSE, col_order = metadata, col_ann_order = metadata, gaps_col = ifelse('stage' %in% metadata, 'stage', meta_col), fontsize = 15)
graphics.off()

##########################################################################################
# plot gene modules with at least 50% of genes DE > 0.25 logFC & FDR < 0.001 (unbiasedGMs_DE)

# Order gms
if (!is.null(metadata_1)){
  antler_data$gene_modules$lists$unbiasedGMs_DE$content <- GeneModuleOrder(seurat_obj = seurat_data, gene_modules = antler_data$gene_modules$lists$unbiasedGMs_DE$content,
                         metadata_1 = metadata_1, order_1 = order_1,
                         metadata_2 = metadata_2, order_2 = order_2,
                         plot_path = "scHelper_log/GM_classification/unbiasedGMs_DE/")
}

ngene = length(unlist(antler_data$gene_modules$lists$unbiasedGMs_DE$content))

# Plot heatmaps
png(paste0(plot_path, 'unbiasedGMs_DE_rownames.png'), height = min(c(150, round(ngene/3))), width = 75, units = 'cm', res = 200)
GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = antler_data$gene_modules$lists$unbiasedGMs_DE$content,
                   show_rownames = TRUE, col_order = metadata, col_ann_order = metadata, gaps_col = ifelse('stage' %in% metadata, 'stage', meta_col),
                   fontsize = 15, fontsize_row = 10)
graphics.off()

png(paste0(plot_path, 'unbiasedGMs_DE.png'), height = min(c(150, ifelse(round(ngene/8) < 20, 20, round(ngene/8)))), width = 75, units = 'cm', res = 600)
GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = antler_data$gene_modules$lists$unbiasedGMs_DE$content,
                   show_rownames = FALSE, col_order = metadata, col_ann_order = metadata, gaps_col = ifelse('stage' %in% metadata, 'stage', meta_col), fontsize = 15, fontsize_row = 10)
graphics.off()


##########################################################################################
# plot gene modules which have been filtered to remove those deferentially expressed across batches
if(length(unique(seurat_data$run)) > 1){
  
  # unbiasedGMs_DE_batchfilt
  ngene = length(unlist(antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content))
  
  # Order gms
  if (metadata_1 %in% c("stage", "scHelper_cell_type", "seurat_clusters")){
    antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content <- GeneModuleOrder(seurat_obj = seurat_data, gene_modules = antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content,
                           metadata_1 = metadata_1, order_1 = order_1,
                           metadata_2 = metadata_2, order_2 = order_2,
                           plot_path = "scHelper_log/GM_classification/unbiasedGMs_DE_batchfilt/")
  }
  
  # Plot heatmaps
  png(paste0(plot_path, 'unbiasedGMs_DE_batchfilt_rownames.png'), height = min(c(150, round(ngene/3))), width = 75, units = 'cm', res = 400)
  GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content,
                     show_rownames = TRUE, col_order = metadata, col_ann_order = metadata, gaps_col = ifelse('stage' %in% metadata, 'stage', meta_col), fontsize = 15, fontsize_row = 10)
  graphics.off()
  
  png(paste0(plot_path, 'unbiasedGMs_DE_batchfilt.png'), height = min(c(150, ifelse(round(ngene/8) < 20, 20, round(ngene/8)))), width = 60, units = 'cm', res = 400)
  GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = antler_data$gene_modules$lists$unbiasedGMs_DE_batchfilt$content,
                     show_rownames = FALSE, col_order = metadata, col_ann_order = metadata, gaps_col = ifelse('stage' %in% metadata, 'stage', meta_col), fontsize = 15, fontsize_row = 10)
  graphics.off()
}


########## Write GMs ##############
export_antler_modules <- function(antler_object, publish_dir, names_list){
  for(gm_list in names_list){
    mods = antler_data$gene_modules$lists[[gm_list]]$content
    for (i in seq(length(mods))) {
      modname = base::names(mods)[i]
      if (is.null(modname)) {
        modname = paste0("GM: ", i)
      }
      write(paste0(modname, "; ", paste0(mods[[i]], collapse = ", ")), file = paste0(publish_dir, '/', gm_list, '.txt'), append = TRUE)
    }
  }
}

export_antler_modules(antler_data, publish_dir = gm_path, names_list = c('unbiasedGMs', 'unbiasedGMs_DE', 'unbiasedGMs_DE_batchfilt'))

########## Save Antler object ##############
saveRDS(antler_data, paste0(rds_path, 'antler_out.RDS'))
