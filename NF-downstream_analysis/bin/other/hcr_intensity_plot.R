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
    data_path = "./input/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/"
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
}

# Set colour map for different HCR combinations
combinations = list("SIX1_TFAP2A_PAX7" = c("SIX1"="magenta", "TFAP2A"="#83f52c", "PAX7"="#ffd700"),
                    "MSX1_SIX1_PAX7" = c("MSX1"="magenta", "SIX1"="#83f52c", "PAX7"="#ffd700"),
                    "DLX6_SIX1_PAX7" = c("DLX6"="magenta", "SIX1"="#83f52c", "PAX7"="#ffd700"))

file_paths <- list.files(data_path, pattern = "*.csv", full.names = TRUE)

# Get HCR combination and genes of interest from path names
hcr_combination <- sub(".*/", "", file_paths[1])
goi <- combinations[[hcr_combination]]

files <- data.frame(name = str_match(file_paths, paste0(hcr_combination, "/\\s*(.*?)\\s*_intensity"))[,2],
                    path = file_paths)

# Read in intensity data
intensity_data <- lapply(files$path, function(x) read.csv(x, col.names = c("M.L.position", names(goi))))
names(intensity_data) <- files$name

# Bind intensity data from different axial levels
intensity_data <- bind_rows(intensity_data, .id = "A.P.position")

# Z-score intensity score for each gene across all A-P regions
intensity_data[,3:5] <- apply(intensity_data[,3:5], 2, scale)

# Transform data to long format for facet plotting
intensity_data <- intensity_data %>% pivot_longer(cols = names(goi), names_to = "gene", values_to = "Scaled Intensity")

plot <- ggplot(intensity_data, aes(y = `Scaled Intensity`, x = `M.L.position`, colour = gene)) +
  geom_line(alpha = 0.4) +
  scale_color_manual(values=goi) +
  facet_wrap(~`A.P.position`, ncol = 1) +
  geom_smooth(span=0.3, se = FALSE) +
  theme_classic() +
  xlab(paste0("M-L position (um)"))

png(paste0(plot_path, hcr_combination,  "_intensity_plot.png"), width = 11, height = 12, units = 'cm', res = 720)
print(plot)
graphics.off()

