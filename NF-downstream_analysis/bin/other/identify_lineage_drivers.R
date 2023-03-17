
#!/usr/bin/env Rscript

# Load packages
library(optparse)
library(mgcv)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(biomaRt)
library(viridis)

# install.packages('tidymv')
# library(tidymv)

# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("", "--verbose"), action = "store_true", type = "logical", help = "Verbose", default = FALSE),
  make_option(c("-t", "--label"), action = "store", type = "character", help = "Label to name output files")
)

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

ncores = opt$cores

# Set paths and load data
plot_path = "./plots/"
rds_path = "./rds_files/"
data_path = "./input"

dir.create(plot_path, recursive = T)
dir.create(rds_path, recursive = T)

# Set output label
if(!is.null(opt$label)){
  label <- opt$label
} else {
  # Retrieve seurat object label
  label <- sub('_.*', '', list.files(data_path, pattern = "*.RDS"))
}


############################################################################################################################################################
# Functions
############################################################################################################################################################

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


PlotLineageVolcano <- function(model_out, ymax = max(-log10(model_out$padj_lineage)), goi_label = NULL, xmin = min(model_out$slope), xmax = max(model_out$slope),
                               padj_lineage_cutoff = 0.05, padj_time_cutoff = 1, sig_shape = 'triangle', up_col = "#228b22", down_col = "#57228b", neg_col = 'gray70'){
  # Filter genes which do not pass time significance cutoff
  model_out <- model_out %>% filter(padj_time < padj_time_cutoff)
  
  # label significance
  model_out <- model_out %>%
    mutate(sig = case_when((padj_lineage < padj_lineage_cutoff & slope > 0) == 'TRUE' ~ 'upregulated',
                           (padj_lineage < padj_lineage_cutoff & slope < 0) == 'TRUE' ~ 'downregulated',
                           (padj_lineage >= padj_lineage_cutoff) == 'TRUE' ~ 'not sig')) %>%
    arrange(abs(padj_lineage))
  
  
  # label outliers with triangles for volcano plot
  model_out <- model_out %>%
    mutate(shape = ifelse(-log10(padj_lineage) > ymax | slope < xmin | slope > xmax, "triangle", "circle")) %>%
    mutate(`-log10(padj lineage)` = ifelse(-log10(padj_lineage) > ymax, ymax, -log10(padj_lineage))) %>%
    mutate(slope = case_when(slope < xmin ~ xmin,
                             slope > xmax ~ xmax,
                             TRUE ~ slope))
  
  plot <- ggplot(model_out, aes(slope, `-log10(padj lineage)`, shape=shape, label = gene)) +
    geom_point(aes(colour = sig, fill = sig), size = 2) +
    scale_color_manual(breaks = c("not sig", "downregulated", "upregulated"),
                       values= c(neg_col, down_col, up_col)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none", legend.title = element_blank()) +
    ylab(paste0('-log10(p-adjusted lineage')) +
    xlab('Coefficient')
  
  if(!is.null(goi_label)){
    plot <- plot + geom_label_repel(data = subset(model_out, gene %in% goi_label), min.segment.length = 0, segment.size  = 0.6, segment.color = "black", max.overlaps = Inf)
  }
  
  return(plot)
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



PrepareLineageGamData <- function(object, gene, slot = 'data', assay = 'RNA', lineage, latent_col = 'latent_time'){
  
  # Get metadata
  metadata <- object@meta.data[, c('latent_time', lineage)]
  colnames(metadata) <- c(latent_col, 'lineage')
  
  # Calculate latent_time_cutoff
  latent_time_cutoff <- CalcLatentTimeCutoff(latent_time = object@meta.data[[latent_col]], lineage_probability = object@meta.data[[lineage]])
  
  # Get expression data
  expression_data <- GetAssayData(object, slot = slot, assay = assay)[gene, rownames(metadata), drop = FALSE] %>% as.matrix() %>% t()
  
  # Filter cells below max_latent_time
  lineage_expression_data <- cbind(metadata, expression_data) %>%
    as.data.frame() %>%
    filter(!!sym(latent_col) < latent_time_cutoff)
  
  return(lineage_expression_data)
}


LineageGam <- function(lineage_expression_data, gene, latent_col = 'latent_time'){
  return(mgcv::gam(data = lineage_expression_data, formula = as.formula(paste0(gene, " ~ s(", latent_col,", bs = 'cs', k = 5)")),
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


RunLineageGamConfidence <- function(object, gene, slot = 'data', assay = 'RNA', lineage, latent_col = 'latent_time'){
  lineage_expression_data <- PrepareLineageGamData(object = object, gene = gene, slot = slot, assay = assay, lineage = lineage, latent_col = latent_col)
  gam <- LineageGam(lineage_expression_data, gene = gene, latent_col = latent_col)
  gam_ci <- CalcGamConfidence(gam)
  gam_ci <- cbind(gam_ci, lineage_expression_data[,latent_col, drop = FALSE])
  gam_ci <- gam_ci %>% mutate(!!sym(latent_col) := round(!!sym(latent_col), digits = 2))
  return(gam_ci)
}


MultiRunLineageGamConfidence <- function(object, gene, slot = 'data', assay = 'RNA', lineage_cols='auto', latent_col = 'latent_time'){
  if(lineage_cols == 'auto'){
    lineage_cols <- grep('lineage_', colnames(object@meta.data), value = TRUE)
  } else if (!all(lineage_cols %in% colnames(object@meta.data))){
    stop('lineage_cols missing from object metadata - please check metadata column names and rerun')
  }
  
  gams <- list()
  for(lineage in lineage_cols){
    gams[[lineage]] <- RunLineageGamConfidence(object = object, gene = gene, lineage = lineage) %>%
      group_by(!!sym(latent_col)) %>% summarise(upper = max(upper), lower = min(lower), fit = mean(fit))
  }
  return(gams)
}

FindMultiGamBreakawayPoint <- function(multi_gam_simulations, latent_time_col = 'latent_time', target_lineage_col, strip_terminal_overlap = TRUE){
  min_latent_time <- lapply(multi_gam_simulations, function(x) max(x[[latent_time_col]])) %>% unlist() %>% min()
  
  multi_gam_simulations <- lapply(multi_gam_simulations, function(x) x %>% filter(!!sym(latent_time_col) <= min_latent_time))
  
  # rearrange gam predictions, setting target lineage first
  multi_gam_simulations <- multi_gam_simulations[c(target_lineage_col, names(multi_gam_simulations)[names(multi_gam_simulations) != target_lineage_col])]
  
  max_breakpoint <- 0
  for(i in 2:length(multi_gam_simulations)){
    merged_df <- merge(multi_gam_simulations[[1]], multi_gam_simulations[[i]], by = latent_time_col)
    
    # Identify which bins of latent time are equivalent are overlapping between two predicted gam CIs
    merged_df <- mutate(merged_df, xy_overlap = ifelse((upper.x <= upper.y & upper.x >= lower.y) | (lower.x <= upper.y & lower.x >= lower.y), 1, 0))
    
    # If there is always an overlap then force to 1
    if (all(merged_df$xy_overlap == 1)){
      breakpoint <- max(merged_df[[latent_time_col]])
    } else if (all(merged_df$xy_overlap == 0)){
      breakpoint <- 0
    } else {
      # Check if lineages overlap again at the end - if they do then the gene is likely involved in secondary process, so strip the trailing 1's
      if(max(which(merged_df$xy_overlap == 1)) == length(merged_df$xy_overlap) & strip_terminal_overlap){
        merged_df <- merged_df[1:max(which(merged_df$xy_overlap != 1)),]
      }
      
      # If there is never an overlap then force to 0
      if(all(merged_df$xy_overlap != 1)){
        breakpoint <- 0
        # If there is always an overlap then force to max latent time
      } else if (all(merged_df$xy_overlap == 1)){
        breakpoint <- max(merged_df[[latent_time_col]])
      } else {
        breakpoint <- which(merged_df$xy_overlap == 1) %>% max()
        breakpoint <- merged_df[[latent_time_col]][breakpoint]
      }
      
      if(breakpoint > max_breakpoint){
        max_breakpoint <- breakpoint
      }
    }
  }
  return(max_breakpoint)
}


###############################################################################################################################################
# Analysis
###############################################################################################################################################



# Load seurat data
# seurat_data <- readRDS('./output/NF-downstream_analysis/transfer_labels/gene_modules_subset_latent_time/rds_files/seurat_label_transfer_latent_time.RDS')
seurat_data <- readRDS(list.files(data_path, pattern = "*.RDS", full.names = TRUE))

# annotations <- read.csv('./output/NF-downstream_analysis/seurat_filtering/preprocess/seurat_annotations.csv')
annotations <- read.csv(list.files(data_path, pattern = "*.csv", full.names = TRUE))

# Get biomart GO annotations for TFs
ensembl <- useEnsembl(biomart = 'genes', 
                      dataset = 'ggallus_gene_ensembl',
                      version = 106)

var_genes <- seurat_data@assays$RNA@var.features
lineages <- grep('lineage_', colnames(seurat_data@meta.data), value = TRUE)

lineage_colours = c('lineage_NC_probability' = '#FFAE49', 'lineage_neural_probability' = '#44A5C2', 'lineage_placodal_probability' = '#024B7A')
lineage_names = c('lineage_NC_probability' = 'neural crest', 'lineage_neural_probability' = 'neural', 'lineage_placodal_probability' = 'placodal')

# var_genes <- c('TFAP2B', 'PAX7', 'MSX1', 'CSRNP1', 'EYA2', 'SIX1')

# Run analysis
# Calculate lineage markers and plot volcanos
for(lineage in lineages){
  lineage_drivers <- ModelMultiLatentExpression(seurat_data, goi = var_genes, latent_col = 'latent_time', lineage_col = lineage)
  # Add gene id from annotations file
  lineage_drivers$gene_id <- annotations[match(lineage_drivers$gene, annotations$Gene), 'Accession']
  
  write.csv(dplyr::select(lineage_drivers, !model), paste0('./', lineage, '_drivers.csv'), quote = FALSE, row.names = FALSE)
  
  # Pull top 10 down and top 20 upregulated
  goi = lineage_drivers %>% filter(!grepl('ENS', gene))
  goi = rbind(filter(goi, slope < 0) %>% arrange(abs(padj_lineage)) %>% head(10), 
        filter(goi, slope > 0) %>% arrange(abs(padj_lineage)) %>% head(20)) %>% pull(gene)
  
  ymax = sort(-log10(lineage_drivers$padj_lineage), decreasing = TRUE)[5]
  xmax = ceiling(min(quantile(lineage_drivers$slope, probs = 0.99), abs(quantile(lineage_drivers$slope, probs = 0.01))))
  xmin = -xmax
  
  png(paste0(plot_path, lineage, '_volcano.png'), width = 20, height = 13, units = 'cm', res = 400)
  print(PlotLineageVolcano(model_out = lineage_drivers, ymax = Inf, xmin = xmin, xmax = xmax, goi_label = c(goi, 'SIX1', 'EYA2', 'ASS1', 'IRF6'), padj_time_cutoff = 0.05))
  graphics.off()
  
  # Filter TFs and make volcano
  TF_subset <- getBM(attributes=c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003"),
                     filters = 'ensembl_gene_id',
                     values = lineage_drivers$gene_id,
                     mart = ensembl)
  
  # subset genes based on transcription factor GO terms
  TF_subset <- TF_subset$ensembl_gene_id[TF_subset$go_id %in% c('GO:0003700', 'GO:0043565', 'GO:0000981', 'GO:0003677', 'GO:0001228')]
  lineage_drivers_TFs <- dplyr::filter(lineage_drivers, gene_id %in% TF_subset)
  
  write.csv(dplyr::select(lineage_drivers_TFs, !model), paste0('./', lineage, '_drivers_TFs.csv'), quote = FALSE, row.names = FALSE)
  
  ymax = sort(-log10(lineage_drivers_TFs$padj_lineage), decreasing = TRUE)[5]
  xmax = ceiling(min(quantile(lineage_drivers_TFs$slope, probs = 0.99), abs(quantile(lineage_drivers_TFs$slope, probs = 0.01))))
  xmin = -xmax
  
  # Pull top 10 down and top 20 upregulated
  goi = lineage_drivers_TFs %>% filter(!grepl('ENS', gene))
  goi = filter(goi, padj_lineage < 0.01 & abs(slope) > 1)
  
  png(paste0(plot_path, lineage, '_volcano_TFs.png'), width = 20, height = 13, units = 'cm', res = 400)
  print(PlotLineageVolcano(model_out = lineage_drivers_TFs, ymax = Inf, xmin = xmin, xmax = xmax, goi_label = goi$gene, padj_time_cutoff = 0.05))
  graphics.off()
  
  saveRDS(lineage_drivers, paste0(rds_path, lineage, "_model_out.RDS"), compress = FALSE)

  ###########################################################################################
  # For each lineage driver, calculate the order  in which they breakaway from other lineages
  
  lineage_drivers_TFs <- lineage_drivers_TFs %>% filter(slope > 0)
  
  simulated_gam_data <- lapply(lineage_drivers_TFs$gene, function(x) MultiRunLineageGamConfidence(seurat_data, x))
  
  names(simulated_gam_data) <- lineage_drivers_TFs$gene
  
  breakpoints <- lapply(simulated_gam_data, FindMultiGamBreakawayPoint, target_lineage_col = lineage, strip_terminal_overlap = TRUE) %>% unlist()
  
  lineage_drivers_TFs[['breakpoints']] <- breakpoints
  lineage_drivers_TFs <- lineage_drivers_TFs %>% mutate(`-log10(padj lineage)` = -log10(padj_lineage))
  
  # Plot lineage drivers ordered by breakpoint
  png(paste0(plot_path, lineage, '_ordered_breakpoint_TFs.png'), width = 20, height = 13, units = 'cm', res = 400)
  ggplot(lineage_drivers_TFs, aes(y = `-log10(padj lineage)`, x = breakpoints, label = gene)) +
    geom_point() +
    geom_label_repel(data = lineage_drivers_TFs, min.segment.length = 0, segment.size  = 0.6, segment.color = "black", max.overlaps = Inf) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none", legend.title = element_blank()) +
    ylab(paste0('-log10(p-adjusted lineage')) +
    xlab('Breakpoint')
  graphics.off()
  
  write.csv(dplyr::select(lineage_drivers_TFs, !model), paste0('./', lineage, '_drivers_TFs_breakpoint.csv'), quote = FALSE, row.names = FALSE)

  
  plot_data <- lapply(simulated_gam_data, function(x) dplyr::bind_rows(x, .id = 'id'))
  
  curr_plot_path <- paste0(plot_path, lineage, "_TF_breakpoints/")
  dir.create(curr_plot_path, recursive = TRUE)
  
  # Plot lineage dynamics with breakpoint for each lineage driver
  for(gene in names(plot_data)){
    plot <- ggplot(plot_data[[gene]], aes(latent_time, fit, colour = id, fill = id)) +
      geom_ribbon(aes(ymax = upper, ymin = lower), alpha=0.4, colour = NA) +
      geom_line() +
      geom_vline(xintercept=breakpoints[gene], linetype="dashed", color = "gray50") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.title = element_text(size=7),
            legend.text = element_text(size=7)) +
      ylab(paste0('Scaled Expression')) +
      xlab('Latent Time') +
      labs(color = "Lineage", 
           fill = "Lineage") +
      scale_fill_viridis(discrete = TRUE, end = 0.98,
                        labels=lineage_names[lineages]) +
      scale_color_viridis(discrete = TRUE, end = 0.98,
                          labels=lineage_names[lineages])
    
    png(paste0(curr_plot_path, gene, ".png"), width = 16, height = 9, units = 'cm', res = 400)
    print(plot)
    graphics.off()

    # Plot UMAP for cells at breakpoint (+-0.005 latent time)
    png(paste0(curr_plot_path, gene, "_breakpoint_umap.png"), width = 16, height = 9, units = 'cm', res = 400)
    DimPlot(seurat_data, cells.highlight = seurat_data@meta.data %>% filter((latent_time < breakpoints[gene] + 0.005) & (latent_time > breakpoints[gene] - 0.005)) %>% rownames(), cols.highlight = 'purple') +
      theme(legend.position = 'none')
    graphics.off()
  }
}

