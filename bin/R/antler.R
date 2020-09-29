#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer",
  'custom_functions', 'm', 2, "character",
  'networkGenes', 'd', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set run location
if(length(commandArgs(trailingOnly = TRUE)) == 0){
  cat('No command line arguments provided, user defaults paths are set for running interactively in Rstudio on docker\n')
  opt$runtype = "user"
} else {
  if(is.null(opt$runtype)){
    stop("--runtype must be either 'user' or 'nextflow'")
  }
  if(tolower(opt$runtype) != "user" & tolower(opt$runtype) != "nextflow"){
    stop("--runtype must be either 'user' or 'nextflow'")
  }
  if(tolower(opt$runtype) == "nextflow"){
    if(is.null(opt$custom_functions) | opt$custom_functions == "null"){
      stop("--custom_functions path must be specified in process params config")
    }
  }
}

# Set paths and load data
{
  if (opt$runtype == "user"){
    sapply(list.files('./bin/R/custom_functions/', full.names = T), source)
    plot.path = "./output/plots/seurat_STACAS/"
    rds.path = "./output/RDS.files/seurat_STACAS/"
    antler.path = "./output/antler/"
    data.path = "./output/RDS.files/seurat_STACAS/"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    sapply(list.files(opt$custom_functions, full.names = T), source)
    plot.path = "./plots/"
    rds.path = "./RDS.files/"
    antler.path = "./antler/"
    data.path = "./"
    
    ncores = opt$cores
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  
  dir.create(plot.path, recursive = T)
  dir.create(rds.path, recursive = T)
  dir.create(antler.path, recursive = T)
  
  # Load packages - packages are stored within renv in the repository
  reticulate::use_python('/usr/bin/python3.7')
  library(Seurat)
  
  library(future)
  library(dplyr)
  library(Antler)
  library(cowplot)
  library(clustree)
  library(gridExtra)
  library(grid)
  library(pheatmap)
  library(biomaRt)
  library(RColorBrewer)
}

# read in seurat data
seurat_out <- readRDS(paste0(data.path, "seurat_out_hh4filt.RDS"))

# strip end of cell names as this is incorrectly reformated in Antler
seurat_out <- RenameCells(seurat_out, new.names = sub('-', '_', colnames(seurat_out)))

seurat_out <- RenameCells(seurat_out, new.names = sub('-.*', '', colnames(seurat_out)))

antler_data <- data.frame(row.names = colnames(seurat_out),
                          "timepoint" = as.numeric(substring(colnames(seurat_out), 3, 3)),
                          "treatment" = rep("null", ncol(seurat_out)),
                          "replicate_id" = rep(1, ncol(seurat_out))
)

# save pheno data
write.table(antler_data, file = paste0(antler.path, "phenoData.csv"), row.names = T, sep = "\t", col.names = T)

# save count data
write.table(GetAssayData(seurat_out, assay = "RNA", slot = "counts"), file = paste0(antler.path, "assayData.csv"), row.names = T, sep = "\t", col.names = T, quote = F)

########################################################################################################
#                            Load Antler data and generate correlation matrix                          #
########################################################################################################

# Change plot path
curr.plot.path <- paste0(plot.path, 'antler/')
dir.create(curr.plot.path)

antler <- Antler$new(output_folder = curr.plot.path, num_cores = ncores)
antler$load_dataset(folder_path = antler.path)
antler$exclude_unexpressed_genes(min_cells=10, min_level=1, verbose=T, data_status='Raw')

antler$normalize(method = 'MR')

antler$gene_modules$identify(
  name                  = "unbiasedGMs",
  corr_t                = 0.3,  # the Spearman correlation treshold
  corr_min              = 3,    # min. number of genes a gene must correlate with
  mod_min_cell          = 10,   # min. number of cells expressing the module
  mod_consistency_thres = 0.4,  # ratio of expressed genes among "positive" cells
  process_plots         = TRUE)

saveRDS(antler, paste0(rds.path, "antler_all.RDS"))

# get automated cluster order based on percentage of cells in adjacent stages
cluster.order = order.cell.stage.clust(seurat_object = seurat_out, col.to.sort = seurat_clusters, sort.by = stage)

# plot all gene modules
png(paste0(curr.plot.path, 'allmodules.png'), height = 100, width = 80, units = 'cm', res = 400)
GM.plot(data = seurat_out, metadata = c("seurat_clusters", "stage"), gene_modules = antler$gene_modules$lists$unbiasedGMs$content,
        show_rownames = F, custom_order = cluster.order, custom_order_column = "seurat_clusters")
graphics.off()

# Plot gene modules with at least 50% of genes DE > 0.25 logFC & FDR < 0.001
DEgenes <- FindAllMarkers(seurat_out, only.pos = T, logfc.threshold = 0.25) %>% filter(p_val_adj < 0.001)
gms <- subset.gm(antler$gene_modules$lists$unbiasedGMs$content, selected_genes = DEgenes$gene, keep_mod_ID = T, selected_gene_ratio = 0.5)

png(paste0(curr.plot.path, 'DE.GM.png'), height = 120, width = 80, units = 'cm', res = 400)
GM.plot(data = seurat_out, metadata = c("seurat_clusters", "stage"), gene_modules = gms, gaps_col = "seurat_clusters",
        show_rownames = T, custom_order = cluster.order, custom_order_column = "seurat_clusters")
graphics.off()

