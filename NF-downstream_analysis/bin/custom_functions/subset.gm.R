
# subset gene modules based on selected gene list provided
# keep_mod_ID sets the names of the module outputs as their initial input position
# selected_gene_ratio is the ratio of selected genes which must be present in order for the module to be kept
subset.gm <- function(gm, selected_genes, keep_mod_ID = F, selected_gene_ratio = 0) {
  outlist <- list()
  if(is.null(names(gm))){
    names(gm) <- paste0("GM: ", 1:length(gm))
  } else {}
  
  # this filters gene modules for which the percentage of differentially expressed genes passes a threshold test
  gm = gm[unlist(lapply(gm, function(x) sum(x %in% selected_genes) >= round(length(x)*selected_gene_ratio)))]
  
  for(mod in names(gm)){
    if(any(selected_genes %in% gm[[mod]])){
      gene.match <- stringr::str_c(selected_genes[which(selected_genes %in% gm[[mod]])], collapse = "; ")
      if(keep_mod_ID == T){
        outlist[[mod]] <- gm[[mod]]
      } else {
        outlist[[gene.match]] <- gm[[mod]]
      }
    }
    else{next}
  }
  return(outlist)
}





