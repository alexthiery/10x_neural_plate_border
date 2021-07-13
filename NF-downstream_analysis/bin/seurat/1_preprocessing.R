# make dataframe with different filtering parameters which can be put into a loop for carrying out downstream analysis
filter_thresholds <- data.frame(gene_min = c(0, 750, 1000, 1500), gene_max = c(Inf, 7000, 6000, 5000), MT_max = c(Inf, 15, 15, 15), row.names = c("unfilt", "low", "med", "high"))

# Plot remaining cells following different filter thresholds
filter_qc <- lapply(rownames(filter_thresholds), function(condition){
  seurat_all@meta.data %>%
    filter(nFeature_RNA > filter_thresholds[condition,'gene_min']) %>%
    filter(nFeature_RNA < filter_thresholds[condition,'gene_max']) %>%
    filter(percent.mt < filter_thresholds[condition,'MT_max']) %>%
    group_by(orig.ident) %>%
    tally() %>%
    rename(!!condition := n)
})

filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc)

png(paste0(plot_path, 'remaining_cell_table.png'), height = 10, width = 18, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(filter_qc, rows=NULL, theme = ttheme_minimal()))
graphics.off()

png(paste0(plot_path, 'remaining_cell_bar.png'), height = 15, width = 21, units = 'cm', res = 400)
ggplot(filter_qc %>% reshape2::melt(), aes(x=variable, y=value, fill=orig.ident)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  xlab("Filter Condition") +
  ylab("Cell Count") +
  ggtitle("Cell count after filtering") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
graphics.off()


# Plot median gene count per cell following different filter thresholds
filter_qc <- lapply(rownames(filter_thresholds), function(condition){
  seurat_all@meta.data %>%
    filter(nFeature_RNA > filter_thresholds[condition,'gene_min']) %>%
    filter(nFeature_RNA < filter_thresholds[condition,'gene_max']) %>%
    filter(percent.mt < filter_thresholds[condition,'MT_max']) %>%
    group_by(orig.ident) %>%
    summarise(median = median(nFeature_RNA, na.rm = TRUE)) %>%
    mutate(median = as.integer(median)) %>%
    rename(!!condition := median)
})

filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc)

png(paste0(plot_path, 'median_gene_count_table.png'), height = 10, width = 18, units = 'cm', res = 400)
grid.arrange(top=textGrob("Median Gene Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(filter_qc, rows=NULL, theme = ttheme_minimal()))
graphics.off()

png(paste0(plot_path, 'median_gene_count_bar.png'), height = 15, width = 21, units = 'cm', res = 400)
ggplot(filter_qc %>% reshape2::melt(), aes(x=variable, y=value, group=orig.ident)) +
  geom_line(aes(colour = orig.ident)) +
  xlab("Filter Condition") +
  ylab("Median Gene Count") +
  ggtitle("Median gene count per cell after filtering") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
graphics.off()


# Plot median gene count on simulated minimum filter thresholds
tests <- data.frame(gene_cutoff = seq(from = 0, to = 3000, by = 10))

filter_qc <- lapply(seq(from = 0, to = 5000, by = 10), function(cutoff){
  seurat_all@meta.data %>%
    filter(nFeature_RNA > cutoff) %>%
    group_by(orig.ident) %>%
    summarise(median = median(nFeature_RNA, na.rm = TRUE)) %>%
    mutate(median = as.integer(median)) %>%
    rename(!! paste(cutoff) := median)
})

filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc) %>% reshape2::melt() %>% mutate(variable = as.integer(variable)*10)

png(paste0(plot_path, 'median_gene_count_bar.png'), height = 15, width = 21, units = 'cm', res = 400)
ggplot(filter_qc, aes(x=variable, y=value, group=orig.ident)) +
  geom_line(aes(colour = orig.ident)) +
  xlab("Lower Gene Threshold") +
  ylab("Median Gene Count") +
  ggtitle("Median gene counts at simulated minimum filter thresholds") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
graphics.off()
