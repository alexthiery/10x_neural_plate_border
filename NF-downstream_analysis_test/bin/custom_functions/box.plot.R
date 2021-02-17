
#### box plot for meta data stats ####
box.plot <- function(dat, y_col, group_by, y_lab, x_lab){
  ggplot(dat, aes(x = dat[[group_by]], y = dat[[y_col]], fill = dat[[group_by]])) +
    geom_boxplot() +
    xlab(x_lab) +
    ylab(y_lab) +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.background = element_rect(colour = "white", fill = "white"),
      strip.text = element_text(size = 18),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18))
}
