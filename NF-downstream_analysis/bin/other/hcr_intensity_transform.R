library(tidyverse)
data_path = '~/Desktop/intensity_plots'

# Set colour map for different HCR combinations
combinations = list("SIX1_TFAP2A_PAX7" = c("SIX1"="magenta", "TFAP2A"="#83f52c", "PAX7"="#ffd700"),
                    "MSX1_SIX1_PAX7" = c("MSX1"="magenta", "SIX1"="#83f52c", "PAX7"="#ffd700"),
                    "DLX6_SIX1_PAX7" = c("DLX6"="magenta", "SIX1"="#83f52c", "PAX7"="#ffd700"))

dir_paths <- list.dirs(data_path, recursive = FALSE)

for(dir in dir_paths){
  file_paths <- list.files(dir_paths, pattern = "*.csv", full.names = TRUE)
  
  # Get HCR combination and genes of interest from path names
  hcr_combination <- sub(".*/", "", dir_paths)
  goi <- combinations[[hcr_combination]]
  
  files <- data.frame(name = str_match(file_paths, paste0(hcr_combination, "/\\s*(.*?)\\s*_intensity"))[,2],
                      path = file_paths)
  
  # Read in intensity data
  intensity_data <- lapply(files$path, function(x) read.csv(x, col.names = c("M.L.position", names(goi))))
  names(intensity_data) <- files$name
  # Bind intensity data from different axial levels
  intensity_data <- bind_rows(intensity_data, .id = "A.P.position")  
  
  write.csv(intensity_data, paste0("~/Desktop/", hcr_combination, "_intensity.csv"), row.names = FALSE)
}
