# plot dimplots for every stage from different datasets
check.integration = function(seurat_object, group.by = 'orig.ident', pt.size = 0.2, xlim = c(-15,15), ylim = c(-15,15), nrow = 1){
  plots = list()
  plots[['all cells']] <- DimPlot(seurat_object, group.by = group.by, pt.size = 0.2) + xlim(xlim) + ylim(ylim)
  
  stages = unique(gsub('_.*', '', seurat_object[[group.by, drop=T]]))
  
  for(i in stages){
    dat <- filter(seurat_object@meta.data, grepl(pattern = i, orig.ident))
    if(length(unique(dat[[group.by]])) > 1){
      plots[[i]] <- DimPlot(seurat_object, cells = rownames(dat), group.by = group.by, pt.size = 0.2) + xlim(xlim) + ylim(ylim) + labs(title = i)
    }
  }
  do.call("grid.arrange", c(plots, nrow=nrow))
}
