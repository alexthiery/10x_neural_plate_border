
calc_latent_time_cutoff <- function(latent_time, lineage_probability, top_frac = 0.2, return = 'intercept'){
  data = data.frame(latent_time, lineage_probability)
  
  model_data <- data %>%
    filter(lineage_probability > 0 & lineage_probability < 1) %>%
    filter(lineage_probability > quantile(lineage_probability, 1-top_frac)) %>%
    # If cells remaining in top frac have less than 0.5 probability of giving rise to the lineage then remove (required for short lineages)
    filter(lineage_probability > 0.5)
  
  # Fit model to top frac of data
  fit = lm(lineage_probability ~ latent_time, data = model_data)
  
  # Inverse equation to find X for Y = 1
  x <- (1-coef(fit)[1])/coef(fit)[2]
  
  # Identify max latent time value based on cells which have reached lineage probability == 1
  max_latent_time <- data %>% filter(lineage_probability == 1) %>% filter(latent_time == max(latent_time)) %>% dplyr::pull(latent_time)
  
  if(return == 'plot'){
    p = ggplot(data, aes(latent_time, lineage_probability)) + 
      geom_point(size = 0.1) +
      geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2])
    
    return(p)
  }else if(return == 'intercept'){
    if(x > max_latent_time){
      cat('Predicted latent time is later than any cell observed that has reached full lineage absorbtion. Max latent time is set based on oldest cell at lineage_probability == 1')
      return(max_latent_time)
    }else{
      return(unname(x))
    }
  }else{
    stop('return must be one of intercept or plot')
  }
}


# Function for generating dataframe for gene expression dynamics
lineage_gam <- function(gene, seurat_object, assay = 'RNA', slot = 'scale.data', lineages, lineage_cutoffs = NULL){
  
  if(!all(names(lineage_cutoffs) %in% lineages)){
    stop('lineage_cutoffs names do not match selected lineages')
  }
  
  assay_dat <- as.data.frame(GetAssayData(seurat_object, assay = assay, slot = slot)[gene,, drop=FALSE])
  
  subset_lineage_probabilities <- c('Neural' = 'lineage_neural_probability', 'Neural crest' = 'lineage_NC_probability', 'Placodal' = 'lineage_placodal_probability')
  
  cols_to_sub <- c('latent_time', unname(subset_lineage_probabilities[lineages]))
  
  assay_dat <- seurat_object@meta.data %>%
    dplyr::select(cols_to_sub) %>%
    merge(t(assay_dat), ., by = 0)
  
  plot_data <- assay_dat %>%
    column_to_rownames('Row.names') %>%
    pivot_longer(!all_of(cols_to_sub)) %>%
    rename(scaled_expression = value) %>%
    rename(gene = name) %>%
    pivot_longer(cols = !c(latent_time, gene, scaled_expression)) %>%
    rename(lineage_probability = value) %>%
    rename(lineage = name) %>%
    group_by(lineage) %>%
    mutate(lineage = names(subset_lineage_probabilities)[subset_lineage_probabilities %in% lineage])
  
  # GAMs
  gams <- plot_data %>%
    group_by(lineage) %>%
    do(gams = gam(scaled_expression ~ s(latent_time, bs = "cs", k=5), weights = lineage_probability, data = .))
  
  # Filter expression values which are above the predicted lineage cutoffs 
  if(!is.null(lineage_cutoffs)){
    plot_data <- plot_data %>%
      mutate(max_latent_time = lineage_cutoffs[names(lineage_cutoffs) %in% lineage]) %>%
      filter(latent_time < max_latent_time)
  }
  
  # Add module column and max latent time to gam data
  plot_data <- plot_data %>%
    ungroup() %>%
    dplyr::select(gene, lineage, max_latent_time) %>%
    distinct()
  
  gams <- inner_join(gams, plot_data)
  
  # Generate predicted values for each gam in tidy df -> output in long format
  out_data <- data.frame()
  for(row in 1:nrow(gams)){
    # Generate latent time values to predict gams on -> use max latent_time calculated per lineage
    pdat <- tibble(latent_time = seq(0, gams[[row, 'max_latent_time']], length = 100))
    new_data <- predict.gam(gams[[row,'gams']][[1]], newdata = pdat, se=TRUE)
    
    out_data <- rbind(out_data, data.frame(gene = gams[[row, 'gene']],
                                             lineage = gams[[row, 'lineage']],
                                             scaled_expression = new_data[['fit']],
                                             se = new_data[['se.fit']],
                                             pdat))
  }
  
  return(out_data)
}
