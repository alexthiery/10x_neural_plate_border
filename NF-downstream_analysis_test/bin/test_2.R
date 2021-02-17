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
}

# Set paths and load data
{
  if (opt$runtype == "user"){
    
    input_path = "./NF-downstream_analysis_test/test_data/"
    output_path = "./output/NF-downstream_analysis_test/test_2/"
    dir.create(output_path, recursive = T)
    
    sapply(list.files('./NF-downstream_analysis_test/bin/custom_functions/', full.names = T), source)
    
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    input_path = "./input/"
    output_path = "./"
    
    sapply(list.files(opt$custom_functions, full.names = T), source)
    
  }
  
}


test_2_data <- readRDS(paste0(input_path, 'test_1.RDS'))

write.csv(test_2_data, paste0(output_path, 'test_2_out.csv'))

