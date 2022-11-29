
#!/usr/bin/env Rscript

# Load packages
library(optparse)
library(mgcv)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggrepel)

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

# Load seurat data
seurat_data <- readRDS(list.files(data_path, pattern = "*.RDS", full.names = TRUE))

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
    plot <- plot + geom_text_repel(data = subset(model_out, gene %in% goi_label), min.segment.length = 0, segment.size  = 0.6, segment.color = "black")
  }
  
  return(plot)
}

# Calculate lineage markers and plot volcanos
var_genes <- seurat_data@assays$RNA@var.features
lineages <- grep('lineage_', colnames(seurat_data@meta.data), value = TRUE)

for(lineage in lineages){
  lineage_drivers <- ModelMultiLatentExpression(seurat_data, goi = var_genes[1:10], latent_col = 'latent_time', lineage_col = lineage)
  goi = lineage_drivers %>% filter(!grepl('ENS', gene)) %>% arrange(abs(padj_lineage)) %>% head(10) %>% pull(gene)
  
  write.csv(dplyr::select(lineage_drivers, !model), paste0('./', lineage, '_drivers.csv'))
  
  png(paste0(plot_path, lineage, '_temp.png'), width = 20, height = 13, units = 'cm', res = 400)
  print(PlotLineageVolcano(model_out = lineage_drivers, ymax = 20, xmin = -20, xmax = 20, goi_label = goi, padj_time_cutoff = 0.05))
  graphics.off()

  saveRDS(lineage_drivers, paste0(rds_path, lineage, "_model_out.RDS"), compress = FALSE)
}