output_path = "./output/"
dir.create(output_path, recursive = T)

t <- matrix(1:4, nrow = 2, ncol = 2)

write.csv(t, paste0(output_path, "test_2_output.csv"))