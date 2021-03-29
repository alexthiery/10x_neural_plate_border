# Function for plotting gene module heatmap from seurat 10x object

# col_order is a single element or a vector element. By default cells are ordered by the first element of the metadata variable.
# User can specify col_order as a vector, with cell order prioritised based on the order of the col_order variable.

# A custom_order can be provided, however if so the custom_order_column must also be specified.

# By default cells are coloured based on the default ggplot colours (to match seurat cluster colours).
# This behaviour can be changed by changing use_seurat_colours and colour_scheme (RColourBrewer)

# A list of gene modules must be provided - this can be filtered by a set of selected_genes if desired

# col_ann_order can be specified to change the order in which the column annotations appear on the heatmap

GM.plot <- function(data, metadata, col_order = metadata[1], custom_order = NULL, custom_order_column = NULL, assay = "RNA", slot = "scale.data",
                    gene_modules, selected_genes = NULL, main = '', hide_annotation = NULL, show_rownames = TRUE, annotation_colors = NA,
                    hclust_rows = FALSE, hclust_cols = FALSE, gaps_row = TRUE, gaps_col = NULL, gm_row_annotation = TRUE, cell_subset = NULL,
                    treeheight_row = 0, use_seurat_colours = TRUE,  colour_scheme = c("PRGn", "RdYlBu", "Greys"),
                    col_ann_order = rev(metadata), ...){
  
  # subset data based on vector of cell names
  if(!is.null(cell_subset)){
    data <- subset(data, cells = cell_subset)
  } else {}
  
  
  # reset levels in seurat_clusters metadata to default numerical order as default
  if("seurat_clusters" %in% metadata){
    data@meta.data[,"seurat_clusters"] <-  factor(data@meta.data[,"seurat_clusters"], levels = sort(unique(as.numeric(as.character(data@meta.data[,"seurat_clusters"])))))
  } else {}
  
  # initialise heatmap metadata
  HM.col <- droplevels(data@meta.data[, metadata, drop=FALSE])
  
  # order HM metadata based on col_order variable
  HM.col <- HM.col[do.call('order', c(HM.col[col_order], list(decreasing=FALSE))), , drop = FALSE]
  
  # order HM metadata based on custom order
  if(!is.null(custom_order)){
    if(is.null(custom_order_column)){
      "custom_order column must be specified \n"
    } else {}
    if(!setequal(custom_order, unique(HM.col[[custom_order_column]]))){
      stop("custom_order factors missing from custom_order_column \n\n")
    } else {}
    HM.col[[custom_order_column]] <-  factor(HM.col[[custom_order_column]], levels = custom_order)
    HM.col <- HM.col[order(HM.col[[custom_order_column]]),,drop = FALSE]  
  }
  
  # subset gene modules by selected gene list if desired
  if(!is.null(selected_genes)){
    selected_GM <- subset.gm(gm = gene_modules, selected_genes = selected_genes)
  } else {
    if(is.null(names(gene_modules))){
      names(gene_modules) <- paste0("GM:", 1:length(gene_modules))
      selected_GM <- gene_modules
    } else {}
    selected_GM <- gene_modules
  }
  
  # gaps_col specifies a metadata column which column gaps are calculated from
  if(!is.null(gaps_col)) {
    if(class(gaps_col) != "character"){
      stop("gaps_col must be a metadata column name")
    } else {
      gaps_col = cumsum(rle(as.vector(HM.col[[gaps_col]]))[["lengths"]])
    }
  } else {
  }
  
  # hide as many annotations in metadata as desired with hide_annotation
  if(!is.null(hide_annotation)){
    HM.col[,hide_annotation] <- NULL
  } else{}
  
  # add row annotations for GM
  if(gm_row_annotation == TRUE) {
    row_ann <- stack(selected_GM)
    rownames(row_ann) <- row_ann$values
    colnames(row_ann)[2] <- "Gene Modules"
    row_ann$values <- NULL
  } else {
    row_ann = NA
  }
  
  # gaps_row splits heatmap rows by gene modules
  if(gaps_row == TRUE) {
    row_ann <- droplevels(row_ann)
    gaps_row = cumsum(summary(as.factor(row_ann[["Gene Modules"]]), maxsum = max(lengths(lapply(row_ann, unique)))))
  } else {
    gaps_row = NULL
  }
  
  
  # Set HM annotation colours
  if(use_seurat_colours == FALSE){
    # set colours for metadata
    ann_colours <- list()
    for(tic in 1:ncol(HM.col)){
      ann_colours[[colnames(HM.col[tic])]] <- setNames(colorRampPalette(brewer.pal(9, colour_scheme[tic])[2:9])(length(unique(HM.col[,tic]))),
                                                       unique(HM.col[,tic]))
    }
  } else {
    # set colours ggplot default colours, as in Seurat::DimPlot
    ann_colours <- list()
    for(column in colnames(HM.col)){
      ann_colours[[column]] <- setNames(ggplotColours(n = length(levels(HM.col[,column]))), levels(droplevels(data@meta.data[, column])))
      
      # change levels of HM col so that heatmap annotations are in the same order as plotted
      ann_colours[[column]] <- ann_colours[[column]][match(levels(HM.col[[column]]), names(ann_colours[[column]]))]
    }
  }
  # asign colours for row annotations
  ann_colours[["Gene Modules"]] <- setNames(colorRampPalette(brewer.pal(9, "Paired"))(length(unique(row_ann$`Gene Modules`))),
                                            unique(row_ann$`Gene Modules`))
  
  
  # extract data to plot from seurat object
  new.dat <- t(as.matrix(x = GetAssayData(object = data, assay = assay, slot = slot)[unlist(selected_GM), rownames(HM.col), drop = FALSE]))
  if(!is.null(cell_subset)){
    cat("rescaling data as cells have been subset \n")
    new.dat <- t(scale(t(new.dat)))
  } else {}
  # add limits to expression as in seurat heatmap plots
  new.dat <- replace(new.dat, new.dat >= 2, 2)
  new.dat <- replace(new.dat, new.dat <= -2, -2)
  
  
  # plot heatmap
  print(pheatmap(t(new.dat), color = PurpleAndYellow(),
                 cluster_rows = hclust_rows, cluster_cols = hclust_cols, show_colnames = FALSE,
                 annotation_col = HM.col[,rev(col_ann_order), drop = FALSE], gaps_col = gaps_col,
                 gaps_row = gaps_row, main = main, show_rownames = show_rownames, annotation_row = row_ann, annotation_colors = ann_colours, treeheight_row = treeheight_row, ...))
}
