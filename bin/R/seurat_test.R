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
    plot.path = "./results/plots/"
    rds.path = "./results/RDS.files/"
    data.path = "./alignment_out"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    sapply(list.files(opt$custom_functions, full.names = T), source)
    plot.path = "plots/"
    rds.path = "RDS.files/"
    data.path = "."
    
    ncores = opt$cores
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  
  dir.create(plot.path, recursive = T)
  dir.create(rds.path, recursive = T)
  
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
}


# read all files from dir
files <- list.files(data.path, recursive = T, full.names = T)
# remove file suffix
file.path <- dirname(files)[!duplicated(dirname(files))]
# make dataframe with tissue matching directory
sample = c("THI300A1" = "hh4_1", "THI300A3" = "ss4_1", "THI300A4" = "ss8_1", "THI300A6" = "hh6_1",
           "THI725A1" = "hh5_2", "THI725A2" = "hh6_2", "THI725A3" = "hh7_2", "THI725A4" = "ss4_2")
matches <- sapply(names(sample), function(x) file.path[grep(pattern = x, x = file.path)])
sample.paths <- data.frame(row.names = sample, sample = sample, tissue = names(matches), path = matches, run = gsub(".*_", "", sample))

# Make Seurat objects for each of the different samples and then merge
seurat_data <- apply(sample.paths, 1, function(x) CreateSeuratObject(counts= Read10X(data.dir = x[["path"]]), project = x[["sample"]]))
seurat_data <- merge(x = seurat_data[[1]], y=seurat_data[-1], add.cell.ids = gsub("_.*", "", names(seurat_data)), project = "chick.10x")

# Remove genes expressed in fewer than 3 cells
seurat_data <- DietSeurat(seurat_data, features = names(which(Matrix::rowSums(GetAssayData(seurat_data) > 0) >=3)))

# Store mitochondrial percentage in object meta data
seurat_data <- PercentageFeatureSet(seurat_data, pattern = "^MT-", col.name = "percent.mt")

# Remove data which do not pass filter threshold
seurat_data <- subset(seurat_data, subset = c(nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 15))


#####################################################################################################
#                           Integrate data from different 10x runs                                  #
#####################################################################################################

# Add metadata col for seq run
seurat_data@meta.data[["seq_run"]] <- gsub(".*_", "", as.character(seurat_data@meta.data$orig.ident))

# Split object by run and find integration points
seurat_data_integrated <- SplitObject(seurat_data, split.by = "seq_run")

# Log normalize data and find variable features
seurat_data_integrated <- lapply(seurat_data_integrated, function(x) NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000))
seurat_data_integrated <- lapply(seurat_data_integrated, function(x) FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000))

plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 3000 * 1024^2)
seurat_data_integrated <- FindIntegrationAnchors(object.list = seurat_data_integrated, dims = 1:30)
seurat_data_integrated <- IntegrateData(anchorset = seurat_data_integrated, dims = 1:30)

# set inegrated count data as default
DefaultAssay(seurat_data_integrated) <- "integrated"

# Scale data and regress out MT content
# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)
seurat_data_integrated <- ScaleData(seurat_data_integrated, features = rownames(seurat_data_integrated), vars.to.regress = "percent.mt")

# Save RDS after scaling as this step takes time
saveRDS(seurat_data_integrated, paste0(rds.path, "seurat_data_integrated.RDS"))




