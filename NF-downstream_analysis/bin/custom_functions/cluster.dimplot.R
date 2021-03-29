# Plot cell subset
cluster.dimplot <- function(seurat_object, group.by = 'seurat_clusters', clusters, cluster_col = 'seurat_clusters', pt.size = 0.2, xlim = c(-15,15), ylim = c(-15,15), nrow = 1, shuffle = TRUE, seed = 123){
  
  plots = list()
  plots[['all cells']] <- DimPlot(seurat_object, group.by = group.by, pt.size = pt.size, shuffle = shuffle, seed = seed) + xlim(xlim) + ylim(ylim) + labs(title = 'all cells')

  # set colours to the same as when plotting all cells
  colours <- factor(seurat_object[[cluster_col, drop = T]])
  
  dat <- filter(seurat_object@meta.data, (!!sym(cluster_col)) %in% clusters)
  plots[['cluster_subset']] <- DimPlot(seurat_object, cells = rownames(dat), pt.size = pt.size, shuffle = shuffle, seed = seed,
          cols = setNames(ggplotColours(n = length(levels(colours))), levels(colours))) + xlim(xlim) + ylim(ylim) + labs(title = 'cluster subset')
  
  do.call("grid.arrange", c(plots, nrow=nrow))
}
