#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer",
  'custom_functions', 'm', 2, "character"
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
    plot_path = "./output/NF-downstream_analysis/1_seurat_integrate/plots/"
    rds_path = "./output/downstream_analysis/1_seurat_integrate/rds_files/"
    data_path = "./output/cellranger/count/filtered_feature_bc_matrix"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    sapply(list.files(opt$custom_functions, full.names = T), source)
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "."
    
    ncores = opt$cores
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
  
  # Load packages - packages are stored within renv in the repository
  reticulate::use_python('/usr/bin/python3.7')
  library(Seurat)
  
  library(future)
  library(dplyr)
  library(cowplot)
  library(clustree)
  library(gridExtra)
  library(grid)
  library(pheatmap)
  library(RColorBrewer)
  library(STACAS)
  library(tidyverse)
}

# make dataframe with stage and replicate taken from path
files <- list.files(data_path, recursive = T, full.names = T)
meta <- str_split(basename(files), pattern = "_", simplify = T)[,1:2]
files_meta <- data.frame(stage = meta[,1], rep = meta[,2], path = files)


 


sample.paths <- data.frame(row.names = sample, sample = sample, stage = names(matches), path = matches, run = gsub(".*-", "", sample))

seurat <- apply(sample.paths, 1, function(x) CreateSeuratObject(counts= Read10X(data.dir = x[["path"]]), project = x[["sample"]]))
seurat <- merge(x = seurat[[1]], y=seurat[-1], add.cell.ids = names(seurat), project = "chick.10x")

# store mitochondrial percentage in object meta data
seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")

