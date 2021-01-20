# plot feature plots for a list of genes
umap.gene.list <- function(data, gene.list, plot.path){
  plot.path = plot.path
  dir.create(plot.path, recursive = TRUE)
  
  missing.genes = c()
  for(g in gene.list){
    if(!g %in% rownames(data)){
      missing.genes = c(g, missing.genes)
      next
    }else{
      print(g)
      pdf(paste0(plot.path, g, "_UMAP.pdf"), height = 5, width = 10)
      plot(gridExtra::grid.arrange(grobs = c(list(DimPlot(data) +
                                                    ggtitle("Seurat clusters") +
                                                    theme(plot.title = element_text(hjust = 0.5))),
                                             list(FeaturePlot(data, g))), ncol = 2))
      dev.off()
    }
  }
  cat('\nFollowing genes not expressed in dataset:', missing.genes, '\n\n')
  # system(paste0("zip -rj ", dirname(plot.path), "/", basename(plot.path), ".zip ", plot.path))
  # unlink(plot.path, recursive=TRUE, force=TRUE)
}