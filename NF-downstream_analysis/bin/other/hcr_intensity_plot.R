
library(tidyverse)
library(ggplot2)

data_path = '~/Desktop/intensity_plots'
plot_path = '~/Desktop/'

# Set colour map for different HCR combinations
combinations = list("SIX1_TFAP2A_PAX7" = c("SIX1"="magenta", "TFAP2A"="#83f52c", "PAX7"="#ffd700"),
                    "MSX1_SIX1_PAX7" = c("MSX1"="magenta", "SIX1"="#83f52c", "PAX7"="#ffd700"),
                    "DLX6_SIX1_PAX7" = c("DLX6"="magenta", "SIX1"="#83f52c", "PAX7"="#ffd700"))

dir_paths <- list.dirs(data_path, recursive = FALSE)

# dir_paths = dir_paths[1]
# Loop through HCR dirs
for(dir in dir_paths){
  file_paths <- list.files(dir, pattern = "*.csv", full.names = TRUE)
  
  # Get HCR combination and genes of interest from path names
  hcr_combination <- sub(".*/", "", dir)
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
}

