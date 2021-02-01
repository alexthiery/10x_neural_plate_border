# !/usr/bin/env Rscript

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
if (opt$runtype == "user"){
  
  # load custom functions
  sapply(list.files('./NF-downstream_analysis/bin/custom_functions/', full.names = T), source) 
  output_path = "./output/NF-downstream_analysis/test/" #we would want this to match the final ouput path that is made when we run in nf
  
  # set cores
  ncores = 8
  
} else if (opt$runtype == "nextflow"){
  cat('pipeline running through nextflow\n')
  
  # load custom functions
  sapply(list.files(opt$custom_functions, full.names = T), source)
  
  # set cores
  ncores = opt$cores
}

library(Seurat)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
#library(rlist) #dont have this package
library(dplyr)
library(pheatmap)
library(gridExtra)
library(grid)
library(reshape2)
library(viridis)
#data.path = ("./input")

setwd("/home/rstudio/NF-downstream_analysis")
data.path = "../alignment_out/10x_scRNAseq"


files <- list.files(data.path, recursive = T, full.names = T)
# remove file suffix
file.path <- dirname(files)[!duplicated(dirname(files))]
# make dataframe with stage matching directory 
### NEED TO MAKE THIS BIT GENERIC
sample = c("THI300A1" = "hh4-1", "THI300A3" = "ss4-1", "THI300A4" = "ss8-1", "THI300A6" = "hh6-1",
           "THI725A1" = "hh5-2", "THI725A2" = "hh6-2", "THI725A3" = "hh7-2", "THI725A4" = "ss4-2")
matches <- sapply(names(sample), function(x) file.path[grep(pattern = x, x = file.path)])

sample.paths <- data.frame(row.names = sample, sample = sample, stage = names(matches), path = matches, run = gsub(".*-", "", sample))

seurat <- apply(sample.paths, 1, function(x) CreateSeuratObject(counts= Read10X(data.dir = x[["path"]]), project = x[["sample"]]))
seurat <- merge(x = seurat[[1]], y=seurat[-1], add.cell.ids = names(seurat), project = "chick.10x")

# store mitochondrial percentage in object meta data
seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")


############# QC plots - Seurat Guided Clustering Tutorial ############
# Visualize QC metrics as a violin plot
png("VlnPlots.png", width=90, height=40, units = 'cm', res = 200)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
graphics.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
png("VlnPlots.png", width=90, height=40, units = 'cm', res = 200)
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
graphics.off()

