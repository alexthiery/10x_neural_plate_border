#### plot multiple feature plots at once ####
multi.feature.plot <- function(seurat.obj, stage.name = NULL, gene.list, plot.clusters = T, plot.stage = F, cluster.col = "seurat_clusters",
                               n.col = 4, label = "UMAP plots for GOI on normalised filtered data", legend.pos = "right", plot.celltype = F,
                               celltype.col = NA){
  if(plot.stage == F & plot.clusters == F){
    tot.len = length(gene.list)
  } else if (plot.stage == T & plot.clusters == T){
    tot.len = length(gene.list) + 2
  } else {
    tot.len = length(gene.list) + 1
  }
  if(plot.celltype == T){
    celltype.plot <- list(DimPlot(seurat.obj, group.by = celltype.col, label = T) + 
                              ggtitle(paste("Cell Types")) +
                              theme(plot.title = element_text(hjust = 0.5)) +
                              NoLegend())
  } else {
    celltype.plot <-NULL
    }
  if(plot.clusters == T){
    clust.plot <- list(DimPlot(seurat.obj, group.by = cluster.col) + 
                           ggtitle(paste("Clusters")) +
                           theme(plot.title = element_text(hjust = 0.5)))
  } else {
    clust.plot <- NULL
    }
  if(plot.stage == T){
    stage.plot <- list(DimPlot(seurat.obj, group.by =  "orig.ident") + 
                           ggtitle("Developmental Stage") +
                           theme(plot.title = element_text(hjust = 0.5)))
    } else {
      stage.plot <- NULL
      }
  plots <- lapply(gene.list, function(x) FeaturePlot(seurat.obj, features = x))
  plots <- c(stage.plot, celltype.plot, clust.plot, plots)
  print(gridExtra::grid.arrange(grobs = plots, ncol = n.col, top = textGrob(label = label, gp=gpar(fontsize=20, font = 2))))
}


