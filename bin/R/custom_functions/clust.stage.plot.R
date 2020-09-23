# Simple function to plot UMAP and dev.stage
clust.stage.plot <- function(data, cluster.col = "seurat_clusters"){
  clust.plot <- DimPlot(data, group.by = cluster.col) + 
    ggtitle(paste("Clusters")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  stage.plot <- DimPlot(data, group.by =  "orig.ident") + 
    ggtitle("Developmental Stage") +
    theme(plot.title = element_text(hjust = 0.5))
  
  gridExtra::grid.arrange(stage.plot, clust.plot, nrow = 1)
}
