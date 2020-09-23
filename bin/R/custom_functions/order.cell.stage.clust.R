
# this function orders cell clusters by a column in the seurat metadata. i.e. sort seurat clusters by their developmental stage, with clusters
# shared between stages ordered in between those stages

# col.to.sort (seurat_clusters) is the metadata column which needs to be re-arranged based on the levels in sort.by (orig.ident)
# only rearrange cluster if at least 10% of cells in cluster are from a different identity


order.cell.stage.clust = function(seurat_object, col.to.sort, sort.by){
  
  enquo_col.to.sort = enquo(col.to.sort)
  enquo_sort.by = enquo(sort.by)
  
  dat <- seurat_object@meta.data %>%
    group_by(!! enquo_sort.by) %>%
    count(!! enquo_col.to.sort) %>%
    arrange(!! enquo_col.to.sort) %>%
    group_by(!! enquo_col.to.sort) %>%
    mutate(n = n/sum(n))
  
  
  top2 <- dat %>%
    top_n(n = 2) %>%
    arrange(!! enquo_sort.by, desc(n)) %>%
    filter(n > 0.1)
  
  
  top1 <- dat %>%
    top_n(n = 1) %>%
    arrange(!! enquo_sort.by, desc(n)) %>%
    distinct(!! enquo_col.to.sort, .keep_all = TRUE)
  
  sort.by = as.character(substitute(sort.by))
  col.to.sort = as.character(substitute(col.to.sort))
  
  for(i in 1:nrow(top1)){
    
    tempnext = unique(as.integer(top2[[sort.by]][as.character(top2[[sort.by]]) == as.character(top1[[sort.by]])[i]]) + 1)
    tempprev = unique(as.integer(top2[[sort.by]][as.character(top2[[sort.by]]) == as.character(top1[[sort.by]])[i]]) - 1)
    
    if(sum(as.integer(top2[[sort.by]]) == tempnext) == 0 & sum(as.integer(top2[[sort.by]]) == tempprev) == 0){
      next
    } else if (sum(as.integer(top2[[sort.by]]) == tempnext) != 0){
      nextfactor = top2[as.integer(top2[[sort.by]]) == tempnext,]
      
      if(top1[[col.to.sort]][i] %in% nextfactor[[2]]){
        top1[i,"n"] <- top1[i,"n"] - 1
      } else {next}
      
      
    } else if (sum(as.integer(top2[[sort.by]]) == tempprev) != 0){
      prevfactor = top2[as.integer(top2[[sort.by]]) == tempprev,]
      
      if(top1[[col.to.sort]][i] %in% prevfactor[[2]]){
        top1[i,"n"] <- (1- top1[i,"n"]) + 1
      } else {next}
      
    } else {stop("error - factor present in more than two levels")}
  }
  
  top1$n = top1$n + (rev(order(levels(top1[[sort.by]])))[top1[[sort.by]]] * 2)
  
  top1 = top1 %>%
    arrange(!! enquo_sort.by, desc(n))
  
  return(top1[[col.to.sort]])
}



