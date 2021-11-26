dir.create('csv')
dir.create('rds')

t <- matrix(1:9, nrow = 3, ncol = 3)

write.csv(t, "csv/test_1_output.csv")
write.csv(t, "rds/test_1_output.csv")

