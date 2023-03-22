CalcLatentTimeCutoff <- function(latent_time, lineage_probability, top_frac = 0.2, return = 'intercept', verbose=FALSE){
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
    if(x > max_latent_time & verbose == TRUE){
      cat('Predicted latent time is later than any cell observed that has reached full lineage absorbtion. Max latent time is set based on oldest cell at lineage_probability == 1')
      return(max_latent_time)
    }else{
      return(unname(x))
    }
  }else{
    stop('return must be one of intercept or plot')
  }
}



PrepareLineageGamData <- function(object, gene, slot = 'data', assay = 'RNA', lineage, latent_col = 'latent_time'){
  
  # Get metadata
  metadata <- object@meta.data[, c('latent_time', lineage)]
  colnames(metadata) <- c(latent_col, 'lineage')
  
  # Calculate latent_time_cutoff
  latent_time_cutoff <- CalcLatentTimeCutoff(latent_time = object@meta.data[[latent_col]], lineage_probability = object@meta.data[[lineage]])
  
  # Get expression data
  expression_data <- GetAssayData(object, slot = slot, assay = assay)[gene, rownames(metadata), drop = FALSE] %>% as.matrix() %>% t()
  
  # Filter cells below max_latent_time
  lineage_expression_data <- cbind(metadata, expression_data) %>%
    as.data.frame() %>%
    filter(!!sym(latent_col) < latent_time_cutoff)
  
  return(lineage_expression_data)
}


LineageGam <- function(lineage_expression_data, goi, latent_col = 'latent_time'){
  # Remove - from goi as it is a special character in mgcv
  if(grepl("-", goi)){
    names(lineage_expression_data)[names(lineage_expression_data) == goi] <- sub("-", "_", goi)
    goi <- sub("-", "_", goi)
  }
  
  return(mgcv::gam(data = lineage_expression_data, formula = as.formula(paste0(goi, " ~ s(", latent_col,", bs = 'cs', k = 5)")),
                   weights = lineage, family = nb(link='log'), method = 'REML', control = gam.control(maxit = 10000)))
}

CalcGamConfidence <- function(gam){
  # Line of 'best fit'
  fit_gam <- predict.gam(gam, type = 'response')
  
  fit_gam_se <- predict.gam(gam, type = 'response', se.fit = TRUE)
  
  # 1.96*se
  fit_gam_ci <- fit_gam_se %>% as.data.frame() %>%
    mutate(upper = fit + (1.96*se.fit)) %>%
    mutate(lower = fit - (1.96*se.fit))
  
  return(fit_gam_ci)
}


RunLineageGamConfidence <- function(object, gene, slot = 'data', assay = 'RNA', lineage, latent_col = 'latent_time'){
  lineage_expression_data <- PrepareLineageGamData(object = object, gene = gene, slot = slot, assay = assay, lineage = lineage, latent_col = latent_col)
  gam <- LineageGam(lineage_expression_data, gene = gene, latent_col = latent_col)
  gam_ci <- CalcGamConfidence(gam)
  gam_ci <- cbind(gam_ci, lineage_expression_data[,latent_col, drop = FALSE])
  gam_ci <- gam_ci %>% mutate(!!sym(latent_col) := round(!!sym(latent_col), digits = 2))
  return(gam_ci)
}


MultiRunLineageGamConfidence <- function(object, gene, slot = 'data', assay = 'RNA', lineage_cols='auto', latent_col = 'latent_time'){
  if(lineage_cols == 'auto'){
    lineage_cols <- grep('lineage_', colnames(object@meta.data), value = TRUE)
  } else if (!all(lineage_cols %in% colnames(object@meta.data))){
    stop('lineage_cols missing from object metadata - please check metadata column names and rerun')
  }
  
  gams <- list()
  for(lineage in lineage_cols){
    gams[[lineage]] <- RunLineageGamConfidence(object = object, gene = gene, lineage = lineage) %>%
      group_by(!!sym(latent_col)) %>% summarise(upper = max(upper), lower = min(lower), fit = mean(fit))
  }
  return(gams)
}


CoexpressionUMAP <- function(seurat_object, gene_1, gene_2, col.threshold = 0, two.colors = c('#FF0000', '#00ff00'), negative.color = 'gray80',
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


