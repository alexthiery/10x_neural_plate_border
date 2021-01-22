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
  

t <- matrix(1:9, nrow = 3, ncol = 3)

write.csv(t, "test_1_output.csv")