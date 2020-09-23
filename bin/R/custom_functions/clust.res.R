
#### cluster at different resolutions and plot clustree
clust.res <- function(seurat.obj, by = 0.1, starting_res = 0){
  plots <- list()
  resolutions <- c(seq(starting_res, starting_res+9*by, by=by))
  if(length(seurat.obj@reductions) == 0){
    stop("Carry out dimensionality reduction (PCA) before clustering")
  }
  seurat.obj@meta.data <- seurat.obj@meta.data[,!grepl("RNA_snn_res.", colnames(seurat.obj@meta.data))]
  seurat.obj <- FindClusters(seurat.obj, resolution = resolutions, verbose = F)
  plots[["clustree"]] <- clustree(seurat.obj@meta.data, prefix = "RNA_snn_res.")
  for(res in resolutions[2:length(resolutions)]){
    plots[[paste(res)]] <- DimPlot(seurat.obj, group.by =  paste0("RNA_snn_res.", res)) +
      ggtitle(paste("resolution = ", res))
  }
  lay <- rbind(c(1,1,1,2,3,4),
               c(1,1,1,5,6,7),
               c(1,1,1,8,9,10))
  plots2 <- gridExtra::arrangeGrob(grobs = plots, layout_matrix = lay)
  return(gridExtra::grid.arrange(plots2))
}
