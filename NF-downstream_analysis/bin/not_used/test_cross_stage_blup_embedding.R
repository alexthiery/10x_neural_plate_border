
CheckNormalisation <- function(seurat_data, assay="RNA") {
  return(!all(seurat_data[[assay]]@counts@x == seurat_data[[assay]]@data@x))
}

ConditionalNormalisation <- function(seurat_data, assay="RNA", normalization.method = "LogNormalize", scale.factor = 10000, ...){
  if(!CheckNormalisation(seurat_data, assay=assay)){
    seurat_data <-  NormalizeData(seurat_data, assay=assay, normalization.method = normalization.method, scale.factor = scale.factor, ... )
  }
  return(seurat_data)
}


# Function to determine features to use for integration - can be either all shared features between seurat objects, shared variable features (determined by nfeatures), or a vector of selected features
GeneralisedIntegrationFeatureSelection <- function(seurat_list, integration_features = "variable", nfeatures = 2000, ...){
  # find all shared features between datasets
  shared_features <- Reduce(intersect, lapply(seurat_list, function(x) rownames(GetAssayData(x))))
  
  if (length(integration_features) > 1){
    
    missing_features <- integration_features[!integration_features %in% shared_features]
    present_features <- integration_features[integration_features %in% shared_features]
    
    if (length(present_features) == 0){
      stop("\nNo selected features are present in all datasets. Please select a valid list of features to integrate, or alternatively integrate 'all' or 'variable' features.\n")
      
    } else if (length(missing_features) > 0){
      cat("\nSelected features which are not present in all datasets will be removed for integration:", paste0(missing_features, collapse = ", "), "\n")
      
    }
    features <- present_features
    
  } else if (integration_features == "variable"){
    if(is.null(nfeatures)){stop("\nnfeatures must be set when selecting shared variable features.\n")}
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = nfeatures, ...)
    
  } else if (integration_features == "all"){
    cat("\n", length(shared_features), "features are shared between datasets. All shared features are being used for integration.\n")
    features <- shared_features
    
  } else {
    stop("\n'integration_features' must be either 'all', 'variable' or a vector of features to integrate.\n")
  }
  return(features)
}



IntegrateSeurat <- function(seurat_list, nfeatures = 2000, dims = 1:20, k.filter = 200, reduction = "cca", integration_features = "variable", verbose = TRUE){
  seurat_list <- lapply(seurat_list, ConditionalNormalisation, verbose = verbose)
  seurat_list <- lapply(seurat_list, FindVariableFeatures, selection.method = "vst", nfeatures = nfeatures, verbose = verbose)
  
  features <- GeneralisedIntegrationFeatureSelection(seurat_list, integration_features = integration_features, nfeatures = nfeatures, verbose = verbose)
  anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features, dims = dims, k.filter = k.filter, reduction = reduction, verbose = verbose)
  
  # this command creates an 'integrated' data assay
  seurat_integrated <- IntegrateData(anchorset = anchors, dims = dims, features = features, verbose = verbose)
  DefaultAssay(object = seurat_integrated) <- "integrated"
  return(seurat_integrated)
}


ScaleClusterSeurat <- function(seurat_data, vars_to_regress = NULL, assay = "RNA", features = rownames(seurat_data), dims = "auto", resolution = 0.5, umap_components = 2, verbose = TRUE){
  
  DefaultAssay(seurat_data) <- assay
  
  seurat_data <- ScaleData(object = seurat_data, vars.to.regress = vars_to_regress, features = features, verbose = verbose)
  seurat_data <- RunPCA(object = seurat_data, verbose = verbose)
  
  # Automatically determine cutoff if desired
  if(length(dims) == 1 && dims == "auto"){
    dims <- 1:scHelper::ElbowCutoff(seurat_data)
  }
  
  seurat_data <- FindNeighbors(object = seurat_data, dims = dims, verbose = verbose)
  seurat_data <- FindClusters(object = seurat_data, resolution = resolution, verbose = verbose)
  seurat_data <- RunUMAP(object = seurat_data, reduction = "pca", dims = dims, n.components = umap_components, verbose = verbose)
  return(seurat_data)
}



RbindSelectedSeuratMeta <- function(object_list, cols=NULL, ident_name = "id"){
  # If unnamed - name objects based on index
  if(!length(names(object_list)) == length(object_list)){
    names(object_list) <- 1:length(object_list)
  }
  
  meta <- lapply(names(seurat_list), function(x) {
    dat <- seurat_list[[x]]@meta.data[,cols ,drop=FALSE]
    dat[[ident_name]] <- x
    return(dat)
  })
  meta <- do.call("rbind", meta)
  return(meta)
}


FindJointNN <- function(object_a, object_b, dims = 'auto', cell_ids = colnames(object_a), k = 5, iterate = 500, frac_subsample = 0.8, verbose = TRUE, ...){
  require(FNN)
  
  # Check that source cell_ids are in all in object_a
  if(!all(cell_ids %in% colnames(object_a))){stop('cell_ids are missing from object_a.')}
  
  # Code to run integration and scaling
  integrated_object <- IntegrateSeurat(list(object_a, object_b), reduction = "rpca", integration_features = "variable", verbose = verbose, ...)
  integrated_object <- ScaleClusterSeurat(integrated_object, assay = "integrated", dims = dims, verbose = verbose)
  
  # Automatically determine cutoff if desired
  if(length(dims) == 1 && dims == "auto"){
    dims <- 1:scHelper::ElbowCutoff(integrated_object)
  }
  
  # Extract joint embeddings from integrated object
  joint_embedding <- as.data.frame(Embeddings(integrated_object, reduction = "pca"))[,dims]
  
  # Only calculate NN using cells of interest (cell_ids) from object_a (in interest of speed)
  if(length(cell_ids) < colnames(object_a) && verbose){
    message('\nCalculating NN for ', sum(cell_ids %in% colnames(object_a)), ' cells_ids from object_a\n')
  }
  
  proportion_iter_match <- IterativeNN(source_ids = cell_ids, target_ids = colnames(object_b), joint_embedding = joint_embedding, k = k, iterate = iterate, frac_subsample = frac_subsample)
  
  return(proportion_iter_match)
}



IterativeNN <- function(source_ids, target_ids, joint_embedding, k = 5, iterate = 500, frac_subsample = 0.8){
  
  # Generate matrix to populate how many observations each target receives on every iteration
  obs_out <- matrix(data = 0, nrow = length(target_ids), ncol = iterate, dimnames = list(target_ids))
  
  for(i in 1:iterate){
    # Get target cells to subset for each iteration
    target_ids_sub <- sample(target_ids, round(length(target_ids)*frac_subsample))
    
    nn <- get.knnx(joint_embedding[target_ids_sub,], joint_embedding[source_ids,], k = k)$nn.index
    
    # Get array of target cell names which are NN for at least 1 source cell
    nn_target_cells <- target_ids_sub[sort(unique(as.vector(nn)))]
    
    # Add 1 for hit
    obs_out[rownames(obs_out) %in% nn_target_cells,i] <- 1
  }
  
  # Calculate the proportion of iterations which a given target cell is matched to at least one source cell
  proportion_iter_match <- apply(obs_out, 1, function(x) sum(x != 0)) / iterate
  
  return(proportion_iter_match) 
}


# Function to re-arrange object based on levels
ReorderObjectList <- function(object_list, levels){
  
  if(length(levels) != length(object_list)){
    stop('levels must be the same length as object_list')
  }
  
  if(is.character(levels)){
    
    if(length(names(object_list)) < length(object_list)){
      stop('levels provided as character vector but object_list contains unnamed elements')
    }
    
    if(!all(levels %in% names(object_list))){
      stop(levels[!levels %in% names(object_list)], ' missing from object list. Please check levels.')
    }
  }
  
  if(is.numeric(levels)){
    if(any(duplicated(levels))){
      stop('levels must not contain duplicated elements: ', levels[duplicated(levels)])
    }
    
    if(max(levels) > length(object_list) | 0 %in% levels){
      stop('Subscript out of bounds. Check levels.')
    }
  }
  
  return(object_list[levels])
}


# Unlist NN matrix output into list of IDs per stage
UnlistNNids <- function(multi_NN_out){
  output <- list(rownames(multi_NN_out[[1]]))
  for(mat in multi_NN_out){
    output <- c(output, list(unname(unlist(as.data.frame(mat)))))
  }
  return(output)
}


# top_perc is the percentage threshold for cells matched cells which are selected for matching with the subsequent stage
MultiJointNN <- function(object_list, cell_ids, levels = NULL, top_perc = 0.5, k = 5, iterate = 500, frac_subsample = 0.8, verbose = TRUE, ...){
  
  if(length(names(object_list)) != length(object_list)){
    message('Names are missing from object list. re-naming objects based on index')
    names(object_list) <- 1:length(object_list)
  }
  
  if(!is.null(levels)){
    object_list <- ReorderObjectList(object_list, levels = levels)
  }
  
  source = 1
  target = 2
  
  output <- list()
  
  # Iterate through pairwise NN for timepoints
  while(target <= length(object_list)){
    # Set output name
    output_name <- paste0(names(object_list)[source], '-', names(object_list)[target], ' NN')
    
    # Find NN
    NN_freq <- FindJointNN(object_a = object_list[[source]], object_b = object_list[[target]], cell_ids = cell_ids,
                           k = k, iterate = iterate, frac_subsample = frac_subsample, verbose = verbose, ...)
    
    output[[output_name]] <- list(source = cell_ids, target_freq = NN_freq)
    
    # Select top n percent of positive matches based on frequency for passing to the next stage
    NN_freq_pos <- NN_freq[NN_freq > 0]
    cell_ids <- names(NN_freq_pos[NN_freq_pos > quantile(NN_freq_pos, prob=1-top_perc)])
    
    # Update ticker
    source = source + 1
    target = target + 1
  }
  
  return(output)
}

NNMapPlot <- function(seurat_list, multi_nn_out, colour_low = "#e1e1e1", colour_high = "#0000d6", reduction = 'umap',
                      label = TRUE, order = TRUE, point_size = 1){
  require(ggplot2)
  require(patchwork)
  
  plots <- list()
  if(length(seurat_list) != (length(multi_nn_out) + 1)){
    stop('Different number of seurat objects in seurat list to multi_nn_out. Please check MultiJoinNN has been run on correct seurat_list')
  }
  
  for(pair_id in names(multi_nn_out)){
    source_id = gsub(pattern = "-.*", "", pair_id)
    target_id = gsub(" .*", "", gsub(pattern = ".*-", "", pair_id))
    
    source = seurat_list[[source_id]]
    target = seurat_list[[target_id]]
    
    # prepare and plot source data
    source = as.data.frame(Embeddings(source, reduction = reduction)[, 1:2])
    source[['source_cells']] = rownames(source) %in% multi_nn_out[[pair_id]]$source
    
    if(order){
      source <- source[order(source$source_cells),]
    }
    
    source_plot <- ggplot(source, aes_string(x = colnames(source)[1], y = colnames(source)[2], colour = "source_cells")) +
      geom_point(size = point_size) +
      scale_color_manual(values = c('FALSE' = colour_low, 'TRUE' = colour_high)) +
      theme_classic() +
      theme(legend.position = "none")
    
    if(label){
      source_plot <- source_plot +
        ggtitle(paste0("Source cells '", source_id, "'")) +
        theme(plot.title = element_text(hjust = 0.5))
    }
    
    # prepare and plot target data
    target = as.data.frame(Embeddings(target, reduction = reduction)[, 1:2])
    target[['target_freq']] = multi_nn_out[[pair_id]]$target_freq[rownames(target)]
    
    if(order){
      target <- target[order(target$target_freq),]
    }
    
    target_plot <- ggplot(target, aes_string(x = colnames(target)[1], y = colnames(target)[2], colour = "target_freq")) +
      geom_point(size = point_size) +
      scale_colour_gradient(low = colour_low, high = colour_high, name="Observed\nFrequency") +
      theme_classic()
    
    if(label){
      target_plot <- target_plot +
        ggtitle(paste0("Target cells '", target_id, "'")) +
        theme(plot.title = element_text(hjust = 0.5))
    }
    plots <- append(plots, list(source_plot))
    plots <- append(plots, list(target_plot))
  }
  return(wrap_plots(plots, ncol = 2))
}


# Cross stage embedding using approach developed by Qiu et al. 2022 (Systematic reconstruction of cellular trajectories across mouse embryogenesis)
library(dplyr)
library(Seurat)
# library(future)
# library(future.apply)
# plan("multiprocess", workers = 4)
# options(future.globals.maxSize = 100000 * 1024^2)

seurat_data <- readRDS('./output/NF-downstream_analysis/transfer_subset/transfer_ppr_nc_subset/seurat/transfer_cluster/rds_files/transfer_clustered_data.RDS')


#################################################################################################################################################
# Extract bins function

extract_bin <- function(seurat_object, gm_1, gm_2, meta_data='scHelper_cell_type', bin_number = 10, bin_extract){
  if(sum(is.null(gm_1), is.null(gm_2)) != 0){
    stop("one or more gene modules lists are empty")
  }
  if(class(gm_1 ) != 'character' | class(gm_1 ) != 'character'){
    stop("one or more gene modules lists are not characters")
  }
  # calculate expression aggregates and products per cell and use that to order them
  x = seurat_object@meta.data[,'scHelper_cell_type']
  gm_1_sum <- t(as.matrix(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data')))[,gm_1] %>% rowSums(.)
  gm_2_sum <- t(as.matrix(GetAssayData(object = seurat_object, assay = 'RNA', slot = 'data')))[,gm_2] %>% rowSums(.)
  plot_data <- data.frame(gm_1_sum = gm_1_sum, gm_2_sum = gm_2_sum, x = x)
  plot_data <- plot_data %>% mutate(ratio = gm_1_sum/(gm_1_sum + gm_2_sum)) %>% arrange(ratio)
  ordered_cells <- rownames(plot_data)
  
  bin_size <- length(ordered_cells)/bin_number
  start = bin_size*(bin_extract[1]-1)
  end = bin_size*(bin_extract[length(bin_extract)])-1
  # end = start + (bin_size-1)
  return(ordered_cells[start:end])
}


antler_data <- readRDS('./output/NF-downstream_analysis/stage_split/ss8_splitstage_data/antler/stage_gene_modules/rds_files/antler_out.RDS')

HH5 <- readRDS('./output/NF-downstream_analysis/stage_split/HH5_splitstage_data/seurat/stage_state_classification/rds_files/HH5_cell_state_classification.RDS')
HH6 <- readRDS('./output/NF-downstream_analysis/stage_split/HH6_splitstage_data/seurat/stage_state_classification/rds_files/HH6_cell_state_classification.RDS')
HH7 <- readRDS('./output/NF-downstream_analysis/stage_split/HH7_splitstage_data/seurat/stage_state_classification/rds_files/HH7_cell_state_classification.RDS')
ss4 <- readRDS('./output/NF-downstream_analysis/stage_split/ss4_splitstage_data/seurat/stage_state_classification/rds_files/ss4_cell_state_classification.RDS')

# Order stages by time (but dont include ss8)
seurat_list <- list(HH5, HH6, HH7, ss4)
stages <- c('HH5', 'HH6', 'HH7', 'ss4')
names(seurat_list) <- stages

# re-scale and cluster individual stages using defined number of PCs
seurat_list_rescaled <- lapply(seurat_list, ScaleClusterSeurat, assay = 'RNA', dims = 'auto', verbose = FALSE)

# subset cells for testing
# seurat_list <- lapply(seurat_list, function(x) subset(x, cells = sample(colnames(x), 500)))

#################################################################################################################################################

# Run coexpression using PPR and NC gene modules from ss8
ppr_gm <- unlist(antler_data$gene_modules$lists$unbiasedGMs_DE$content[c('GM5')])
nc_gm <- unlist(antler_data$gene_modules$lists$unbiasedGMs_DE$content[c('GM2')])

# Extract bin 3 cells
bin_3_cells <- extract_bin(seurat_data, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = 10, bin_extract = 3)

# Filter ss4 bin 3 cells
ss4_blups <- grep('ss4', bin_3_cells, value=TRUE)

multi_NN_output <- MultiJointNN(seurat_list_rescaled,  cell_ids = ss4_blups, levels = rev(stages), dims = 'auto', top_perc = 0.25, iterate = 500, k = 3, frac_subsample = 0.9)

png('./bin3_ss4_backwards_NN.png', width = 20, height = 27, units = 'cm', res = 400)
NNMapPlot(seurat_list = seurat_list, multi_nn_out = multi_NN_output)
graphics.off()


# Extract bin 1 cells
bin_1_cells <- extract_bin(seurat_data, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = 10, bin_extract = 1)

# Filter ss4 bin 1 cells
ss4_NC <- grep('ss4', bin_1_cells, value=TRUE)

multi_NN_output <- MultiJointNN(seurat_list_rescaled,  cell_ids = ss4_NC, levels = rev(stages), dims = 'auto', top_perc = 0.25, iterate = 500, k = 3, frac_subsample = 0.9)

png('./bin1_ss4_backwards_NN.png', width = 20, height = 27, units = 'cm', res = 400)
NNMapPlot(seurat_list = seurat_list, multi_nn_out = multi_NN_output)
graphics.off()


# Extract bin 10 cells
bin_10_cells <- extract_bin(seurat_data, gm_1 = ppr_gm, gm_2 = nc_gm, meta_data = c('scHelper_cell_type'), bin_number = 10, bin_extract = 10)

# Filter ss4 bin 10 cells
ss4_plac <- grep('ss4', bin_10_cells, value=TRUE)

multi_NN_output <- MultiJointNN(seurat_list_rescaled,  cell_ids = ss4_plac, levels = rev(stages), dims = 'auto', top_perc = 0.25, iterate = 500, k = 3, frac_subsample = 0.9)

png('./bin10_ss4_backwards_NN.png', width = 20, height = 27, units = 'cm', res = 400)
NNMapPlot(seurat_list = seurat_list, multi_nn_out = multi_NN_output)
graphics.off()




##################################################################################################################################################################
# Identify NPB and plac descendents


hh5_npb_cells <- seurat_list$HH5@meta.data %>% filter(scHelper_cell_type == "pNPB") %>% rownames()

multi_NN_output <- MultiJointNN(seurat_list_rescaled,  cell_ids = hh5_npb_cells, levels = stages, dims = 'auto', top_perc = 0.1, iterate = 500, k = 3, frac_subsample = 0.9)

png('./npb_hh5_forwards_NN.png', width = 20, height = 27, units = 'cm', res = 400)
NNMapPlot(seurat_list = seurat_list, multi_nn_out = multi_NN_output)
graphics.off()



hh5_ppr_cells <- seurat_list$HH5@meta.data %>% filter(scHelper_cell_type == "aPPR") %>% rownames()

multi_NN_output <- MultiJointNN(seurat_list_rescaled,  cell_ids = hh5_ppr_cells, levels = stages, dims = 'auto', top_perc = 0.1, iterate = 500, k = 3, frac_subsample = 0.9)

png('./ppr_hh5_forwards_NN.png', width = 20, height = 27, units = 'cm', res = 400)
NNMapPlot(seurat_list = seurat_list, multi_nn_out = multi_NN_output)
graphics.off()




# Run from HH6 onwards as HH5 is too heterogeneous
seurat_list_rescaled_sub <- seurat_list_rescaled[-1]

hh6_npb_cells <- seurat_list_rescaled_sub$HH6@meta.data %>% filter(scHelper_cell_type == "pNPB") %>% rownames()

multi_NN_output <- MultiJointNN(seurat_list_rescaled_sub,  cell_ids = hh6_npb_cells, levels = stages[-1], dims = 'auto', top_perc = 0.1, iterate = 500, k = 3, frac_subsample = 0.9)

png('./npb_hh6_forwards_NN.png', width = 20, height = 18, units = 'cm', res = 400)
NNMapPlot(seurat_list = seurat_list[-1], multi_nn_out = multi_NN_output)
graphics.off()



hh6_ppr_cells <- seurat_list_rescaled_sub$HH6@meta.data %>% filter(scHelper_cell_type == "aPPR") %>% rownames()

multi_NN_output <- MultiJointNN(seurat_list_rescaled_sub,  cell_ids = hh6_ppr_cells, levels = stages[-1], dims = 'auto', top_perc = 0.1, iterate = 500, k = 3, frac_subsample = 0.9)

png('./appr_hh6_forwards_NN.png', width = 20, height = 18, units = 'cm', res = 400)
NNMapPlot(seurat_list = seurat_list[-1], multi_nn_out = multi_NN_output)
graphics.off()
