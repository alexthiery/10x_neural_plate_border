
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

coexpression_umap <- function(seurat_object, gene_1, gene_2, col.threshold = 0, two.colors = c('#FF0000', '#00ff00'), negative.color = 'gray80',
                                      highlight_cell_size = 1, show_legend = TRUE,
                                      axes_label_size = 12, axes_title_size = 10, axes_tick_size = 0.15){

  dat <- t(as.matrix(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data')[c(gene_1, gene_2),]))
  dat <- apply(dat, 2, function(x) ifelse(!x, 0, as.integer((x - min(x)) / (max(x) - min(x)) * 100)))
  
  
  col_mat = Seurat:::BlendMatrix(n = 101, col.threshold = col.threshold, two.colors =  two.colors, negative.color = negative.color)
  col_mat <- as.data.frame.table(col_mat, responseName = "value") %>% mutate_if(is.factor, as.numeric)
  # Set base col values to 0
  col_mat[1:2] <- col_mat[1:2] - 1
  col_mat[!(col_mat$Var1 > col.threshold*100 & col_mat$Var2 > col.threshold*100), 'value'] <- negative.color
  
  colnames(col_mat) <- c('a', 'b', 'mix')
  
  # Vectorised subset of col mat based on expression values
  cell_cols <- set_names(col_mat[[3]], paste(col_mat[[1]], col_mat[[2]], sep = '_'))
  cell_cols <- data.frame(cell_cols = unname(cell_cols[paste(dat[,1], dat[,2], sep = '_')]), row.names = rownames(dat))
  
  col_mat[,1:2] <- col_mat[,1:2]/100
  key_plot <- ggplot(col_mat, aes(x = a, y = b)) +
    xlab(gene_1) +
    ylab(gene_2) +
    geom_tile(aes(fill = mix)) +
    scale_fill_identity() +
    scale_x_continuous(breaks = c(0, 1), expand = c(0.01, 0.01)) +
    scale_y_continuous(breaks = c(0, 1), expand = c(0.01, 0.01)) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_text(size = axes_title_size),
          axis.text.x = element_text(size = axes_label_size),
          axis.title.y = element_text(size = axes_title_size),
          axis.text.y = element_text(size = axes_label_size),
          axis.ticks.length=unit(axes_tick_size,"cm"))

  plot_data <- as.data.frame(seurat_object[["umap"]]@cell.embeddings)
  
  plot_data$cell_cols = cell_cols$cell_cols[match(rownames(plot_data), rownames(cell_cols))]

  positive_cells = filter(plot_data, cell_cols != negative.color)

  umap_plot <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, colour = rownames(positive_cells))) +
    geom_point(colour = negative.color, size = 2) +
    geom_point(data = positive_cells, size = highlight_cell_size) +
    scale_colour_manual(breaks = rownames(positive_cells), values=positive_cells$cell_cols) +
    theme_classic() +
    NoLegend()


  layout <- '
    BA
    B#
    '
    return(wrap_plots(A = key_plot, B = umap_plot, design = layout, widths = c(4,1), heights = c(1,3)))
}


