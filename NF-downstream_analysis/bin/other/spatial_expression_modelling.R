library(optparse)
library(tidyverse)
library(ggplot2)

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
    
    ncores = 8
    plot_path = "./plots/"
    data_path = "./work/77/81a9e2a039bf09a187fe75bf317ac2/input"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input"
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
}

# Read in seurat stage data
seurat_data <- readRDS(list.files(data_path, full.names = TRUE))

# Select PC of interest for visualisation
pc_oi = 'PC_1'

# Select genes of interest for spatial modelling
goi = list(
  c("PAX7"="#ffd700", "TFAP2A"="#83f52c", "SIX1"="magenta"),
  c("PAX7"="#ffd700", "SIX1"="#83f52c", "MSX1"="magenta"),
  c("PAX7"="#ffd700", "SIX1"="#83f52c", "DLX6"="magenta")
)

# Group cell states into parent cell states to visualise association of M-L patterning along PCs
parent_cell_states <- c('EE'='EE', 'NNE'='NNE', 'pEpi'='NNE', 'PPR'='PPR', 'aPPR'='PPR', 'pPPR'='PPR',
                              'eNPB'='NPB', 'NPB'='NPB', 'aNPB'='NPB', 'pNPB'='NPB','NC'='NC', 'dNC'='NC',
                              'eN'='Neural', 'eCN'='Neural', 'NP'='Neural', 'pNP'='Neural', 'HB'='Neural', 'iNP'='Neural', 'MB'='Neural', 
                              'aNP'='Neural', 'FB'='Neural', 'vFB'='Neural', 'node'='node', 'streak'='streak')

seurat_data@meta.data$parent_cell_states = parent_cell_states[as.character(seurat_data@meta.data$scHelper_cell_type)]

# Subset embeddings for PC of interest
pc_embeddings <- seurat_data@reductions$pca@cell.embeddings[,pc_oi]

# Prepare plot data
plot_data <- data.frame(cell_type = as.character(seurat_data@meta.data[names(pc_embeddings), 'parent_cell_states']), pc_embeddings = pc_embeddings)

# Plot histogram of cell state distributions along PC of interest
png(paste0(plot_path, pc_oi, "_cell_states.png"), width = 11, height = 6, units = 'cm', res = 720)
ggplot(plot_data, aes(x = pc_embeddings, colour = cell_type, fill = cell_type)) +
  geom_density(alpha = 0.2) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  xlab(sub('_', ' ', pc_oi)) +
  ylab('Density') +
  theme_classic() +
  theme(legend.title=element_blank()) +
  scale_x_reverse()
graphics.off()


for(combination in goi){
  expression_subset <- cbind(plot_data, t(seurat_data@assays$RNA@scale.data[names(combination), rownames(plot_data)])) %>%
    pivot_longer(cols = !c(cell_type, pc_embeddings), names_to = "Gene", values_to = "Scaled Expression")
  
  plot <- ggplot(expression_subset, aes(x = pc_embeddings, y = `Scaled Expression`, group = Gene, colour = Gene)) +
    geom_smooth(se = FALSE, method='gam') +
    scale_color_manual(values=combination) +
    theme_classic() +
    xlab(sub('_', ' ', pc_oi)) +
    theme(legend.title=element_blank()) +
    scale_x_reverse()
  
  png(paste0(plot_path, pc_oi, "_", paste(names(combination), collapse = "_"), "_dynamics.png"), width = 11, height = 6, units = 'cm', res = 720)
  print(plot)
  graphics.off()
}




