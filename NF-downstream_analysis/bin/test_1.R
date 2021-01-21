# !/usr/bin/env Rscript

# # Define arguments for Rscript
# library(getopt)
# spec = matrix(c(
#   'runtype', 'l', 2, "character",
#   'cores'   , 'c', 2, "integer",
#   'custom_functions', 'm', 2, "character"
# ), byrow=TRUE, ncol=4)
# opt = getopt(spec)

# # Set run location
# if(length(commandArgs(trailingOnly = TRUE)) == 0){
#   cat('No command line arguments provided, user defaults paths are set for running interactively in Rstudio on docker\n')
#   opt$runtype = "user"
# } else {
#   if(is.null(opt$runtype)){
#     stop("--runtype must be either 'user' or 'nextflow'")
#   }
#   if(tolower(opt$runtype) != "user" & tolower(opt$runtype) != "nextflow"){
#     stop("--runtype must be either 'user' or 'nextflow'")
#   }
#   if(tolower(opt$runtype) == "nextflow"){
#     if(is.null(opt$custom_functions) | opt$custom_functions == "null"){
#       stop("--custom_functions path must be specified in process params config")
#     }
#   }
# }

# # Set paths and load data
# {
#   if (opt$runtype == "user"){

#     # load custom functions
#     sapply(list.files('./NF-downstream_analysis/bin/custom_functions/', full.names = T), source)
#     output_path = "./output/NF-downstream_analysis/test/"
#     plot_path = "./output/NF-downstream_analysis/test/plots"
    
#     # set cores
#     ncores = 8
    
#   } else if (opt$runtype == "nextflow"){
#     cat('pipeline running through nextflow\n')
    
#     # load custom functions
#     sapply(list.files(opt$custom_functions, full.names = T), source)
#     output_path = "./output/test"
#     plot_path = "./output/test/plots/"
    
#     # set cores
#     ncores = opt$cores
#   }
  
# dir.create(output_path, recursive = T)
# dir.create(plot_path, recursive = T)


output_path = "./output/"
dir.create(output_path, recursive = T)

t <- matrix(1:9, nrow = 3, ncol = 3)

write.csv(t, paste0(output_path, "test_1_output.csv"))