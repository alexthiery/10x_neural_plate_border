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
  

# t <- matrix(1:9, nrow = 3, ncol = 3)
# write.csv(t, "test_1_output.csv")

data.path = ("./input/")

reticulate::use_python('/usr/bin/python3.7')
library(Seurat)

# read all files from dir
files <- list.files(data.path, recursive = T, full.names = T)
# remove file suffix
file.path <- dirname(files)[!duplicated(dirname(files))]
# make dataframe with stage matching directory
sample = c("THI300A1" = "hh4-1", "THI300A3" = "ss4-1", "THI300A4" = "ss8-1", "THI300A6" = "hh6-1",
           "THI725A1" = "hh5-2", "THI725A2" = "hh6-2", "THI725A3" = "hh7-2", "THI725A4" = "ss4-2")
matches <- sapply(names(sample), function(x) file.path[grep(pattern = x, x = file.path)])

sample.paths <- data.frame(row.names = sample, sample = sample, stage = names(matches), path = matches, run = gsub(".*-", "", sample))


# Make Seurat objects for each of the different samples and then merge
seurat_data <- apply(sample.paths, 1, function(x) CreateSeuratObject(counts= Read10X(data.dir = x[["path"]]), project = x[["sample"]]))
seurat_data <- merge(x = seurat_data[[1]], y=seurat_data[-1], add.cell.ids = names(seurat_data), project = "chick.10x")

saveRDS(seurat_data, "test_1_output.RDS")
