library(Seurat)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(rlist)
library(dplyr)
library(pheatmap)
library(gridExtra)
library(grid)
library(reshape2)
library(viridis)

source("rscripts/0_custom_functions.R")

dir.create("RDS.files/1_filter_QC/", recursive = T)
plot.path <- 'plots/1_filter_QC/'
dir.create(plot.path, recursive = T)
# make directory for the rest of the QC plots which directly compare the filtering parameters
plot.path.filtcomp <- "plots/1_filter_QC/filter_comparison/"
dir.create(plot.path.filtcomp, recursive = T)


sample.paths<-data.frame(tissue = c("hh4", "hh6", "ss4", "ss8"),
                         path = c("../cellranger/cellranger_output/cellranger_count_hh4/outs/filtered_feature_bc_matrix_chr_edit/",
                                  "../cellranger/cellranger_output/cellranger_count_hh6/outs/filtered_feature_bc_matrix_chr_edit/",
                                  "../cellranger/cellranger_output/cellranger_count_4ss/outs/filtered_feature_bc_matrix_chr_edit/",
                                  "../cellranger/cellranger_output/cellranger_count_8ss/outs/filtered_feature_bc_matrix_chr_edit/"))


# Make Seurat objects for each of the different samples. The raw data for each sample is found in the relative directory assigned above in sample.paths.
for(i in 1:nrow(sample.paths["path"])){
  name<-paste(sample.paths[i,"tissue"])
  assign(name, CreateSeuratObject(counts= Read10X(data.dir = paste(sample.paths[i,"path"])), project = paste(sample.paths[i, "tissue"])))
}

# The four Seurat objects are then merged, before running CreateSeuratObject again on the output in order to apply the min.cells parameter on the final merged dataset.
temp <- merge(hh4, y = c(hh6, ss4, ss8), add.cell.ids = c("hh4", "hh6", "ss4", "ss8"), project = "chick.10x")
merged.data<-CreateSeuratObject(GetAssayData(temp), min.cells = 3, project = "chick.10x.mincells3")

# The original Seurat objects are then removed from the global environment
rm(hh4, hh6, ss4, ss8, sample.paths, temp)

# store mitochondrial percentage in object meta data
merged.data <- PercentageFeatureSet(merged.data, pattern = "^MT-", col.name = "percent.mt")

# make dataframe with different filtering parameters which can be put into a loop for carrying out downstream analysis
filt_crit_df<-data.frame(gene_min = c(0, 750, 1000, 1250), gene_max = c(Inf, 7000, 6500, 6000), MT_max = c(Inf, 15, 15, 15))
rownames(filt_crit_df)<-c("unfilt", "low", "med", "high")

# save image of filter conditions table
filt_tab <- t(filt_crit_df)[c(3,1,2), 2:4]
rownames(filt_tab) <- c("% max MT content", "min gene count", "max gene count")
pdf(paste0(plot.path, "filter_conditions.pdf"), height = 1.5, width = 3.5)
grid.arrange(top=textGrob("Filter Conditions",gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.3, vjust = 1), tableGrob(filt_tab, theme = ttheme_minimal()))
dev.off()
rm(filt_tab)


# Add columns of metadata for each set of filtering parameters in order to label which cells should be kept/removed
for(condition in rownames(filt_crit_df)){
  merged.data <- add.filter.to.meta(data = merged.data, filt.param.df = filt_crit_df, filt.param.name = condition)
}

# Stats on no. genes expressed in each cell at each stage and no. of genes in >3 cells - pre and post filering
filt.stats <- alternate.cols(gene.stats(data=merged.data, min.cell = 3, filter_type = "unfilt"),
                             gene.stats(data=merged.data, min.cell = 3, filter_type = "low"),
                             gene.stats(data=merged.data, min.cell = 3, filter_type = "med"),
                             gene.stats(data=merged.data, min.cell = 3, filter_type = "high"))

############# plot filtering statistics ############

# save image of remaining cells stats after filtering
filt_tab <- sapply(c("unfilt", "low", "med", "high"),
       function(x) filt.stats["no. remaining cells", grep(x, colnames(filt.stats))])
rownames(filt_tab) <- sapply(strsplit(rownames(filt_tab), " "), "[", 1)
# plot table
pdf(paste0(plot.path, "remaining_cells.pdf"), height = 2, width = 3.5)
grid.arrange(top=textGrob("Remaining Cell Count",gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.3, vjust = 0.5), tableGrob(filt_tab, theme = ttheme_minimal()))
dev.off()
# plot bar plot
filt_tab <- melt(filt_tab[1:4,])
g <- ggplot(filt_tab, aes(x = Var2, y = value, fill = Var1)) +
  xlab("Filter Condition") +
  ylab("Cell Count") +
  labs(fill = "Stage") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = viridis(4)) +
  ggtitle("Cell count after filtering") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
pdf(paste0(plot.path, "remaining_cells_bar.pdf"), height = 5, width = 7)
plot(g)
dev.off()
rm(filt_tab)

# save image of how many cells are remaining after filtering
filt_tab <- sapply(c("unfilt", "low", "med", "high"),
                   function(x) filt.stats["median genes per cell", grep(x, colnames(filt.stats))])
rownames(filt_tab) <- sapply(strsplit(rownames(filt_tab), " "), "[", 1)
filt_tab <- melt(filt_tab[1:4,])

g <- ggplot(filt_tab, aes(x = Var2, y = value, group = Var1)) +
  xlab("Filter Condition") +
  ylab("Median Gene Count") +
  geom_line(aes(colour = Var1)) +
  scale_colour_manual(values = viridis(4)) +
  ggtitle("Median gene count per cell after filtering") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
pdf(paste0(plot.path, "median_rem_genes.pdf"), height = 5, width = 7)
plot(g)
dev.off()
rm(filt_tab)


############# QC plots ############

# plot number of remaining cells and median genes per cell at simulated minimum cutoff thresholds
plots <- plot.gene.count(merged.data, sep.stages = T)
pdf(paste0(plot.path.filtcomp, "gene.count.sim.pdf"))
print(plot_grid(plotlist = plots))
dev.off()


# plot histograms of gene number after applying different thresholds
plots <- hist.gene.count(merged.data, filt.param.df = filt_crit_df)
pdf(paste0(plot.path.filtcomp, "gene.count.hist.pdf"))
grid.arrange(plots$unfilt, plots$low, plots$med, plots$high, nrow = 3, 
                        layout_matrix = cbind(c(1,1,2), cbind(c(1,1,3)), c(1,1,4)))
dev.off()
rm(plots)
rel_heights=c(0.1, 1)


# For each of the different filtering columns in the metadata plot the QC plots before filtering
for(condition in rownames(filt_crit_df)){
  temp.filt_crit_df<-filt_crit_df[condition,]
  temp.plot.path <- paste0(plot.path, "filt_param_", condition, "/")
  dir.create(paste(temp.plot.path))
  # violin plots of cell gene no, UMI count and MT content post filtering
  pdf(paste0(temp.plot.path, "filtered_vioplot.pdf"),width=30,height=15)
  p1 = vio.plot(meta.data = merged.data@meta.data, y.dat = "nCount_RNA", x_split = condition, y.lab = "UMI Count")
  p2 = vio.plot(meta.data = merged.data@meta.data, y.dat = "nFeature_RNA", x_split = condition, y.lab = "Gene Count")
  p3 = vio.plot(meta.data = merged.data@meta.data, y.dat = "percent.mt", x_split = condition, y.lab = "Percent Mitochondrial Genes")
  p = plot_grid(p1, p2, p3, ncol = 3)
  title <- ggdraw() + draw_label(paste0("Filtered Cells (<", temp.filt_crit_df[,"gene_min"], " or ", temp.filt_crit_df[,"gene_max"], " genes and >", temp.filt_crit_df[,"MT_max"], "% MT genes)"),
                                 fontface='bold', size = 20)
  print(plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
  dev.off()
  
  # same plot again, however this time the filtering criteria are combined and coloured instead so that you can see where the cells are on the original data
  pdf(paste0(temp.plot.path, "unfiltered_vioplot.pdf"),width=30,height=15)
  p1 = vio.plot(meta.data = merged.data@meta.data, y.dat = "nCount_RNA", dot.col = condition, y.lab = "UMI Count", fill.vio = F)
  p2 = vio.plot(meta.data = merged.data@meta.data, y.dat = "nFeature_RNA", dot.col = condition, y.lab = "Gene Count", fill.vio = F)
  p3 = vio.plot(meta.data = merged.data@meta.data, y.dat = "percent.mt", dot.col = condition, y.lab = "Percent Mitochondrial Genes", fill.vio = F)
  p = plot_grid(p1, p2, p3, ncol = 3)
  title <- ggdraw() + draw_label(paste0("Unfiltered Cells (<", temp.filt_crit_df[,"gene_min"], " or ", temp.filt_crit_df[,"gene_max"], " genes and >", temp.filt_crit_df[,"MT_max"], "% MT genes)"),
                                 fontface='bold', size = 20)
  print(plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
  dev.off()
  
  # plot feature scatter plots before filtering so that I can see which cells will be removed
  pdf(paste0(temp.plot.path, "featscat.prefilt_UMI.MT.pdf"),width=15,height=8)
  plot1 <- FeatureScatter(merged.data, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = condition, pt.size = 0.1)
  plot2 <- FeatureScatter(merged.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = condition, pt.size = 0.1)
  print(CombinePlots(plots = list(plot1, plot2)))
  dev.off()
}

########### Filter subset ##########

# Subset data based on the different set of filtering parameters and plot QC after filtering
filt.object.list <- list()
for(condition in rownames(filt_crit_df)){
  temp.filt_crit_df<-filt_crit_df[condition,]
  temp.plot.path <- (paste0(plot.path, "filt_param_", condition, "/"))
  
  # subset data
  filt.data<-subset(merged.data, subset = c(nFeature_RNA > temp.filt_crit_df[,"gene_min"] & nFeature_RNA < temp.filt_crit_df[,"gene_max"] & percent.mt < temp.filt_crit_df[,"MT_max"]))
  
  # plot only remaining cells
  pdf(paste0(temp.plot.path, "postfiltered_vioplot.pdf"),width=30,height=15)
  p1 = vio.plot(meta.data = filt.data@meta.data, y.dat = "nCount_RNA", y.lab = "UMI Count")
  p2 = vio.plot(meta.data = filt.data@meta.data, y.dat = "nFeature_RNA", y.lab = "Gene Count")
  p3 = vio.plot(meta.data = filt.data@meta.data, y.dat = "percent.mt", y.lab = "Percent Mitochondrial Genes")
  p = plot_grid(p1, p2, p3, ncol = 3)
  title <- ggdraw() + draw_label(paste0("Post Filtering (<", temp.filt_crit_df[,"gene_min"], " or ", temp.filt_crit_df[,"gene_max"], " genes and >", temp.filt_crit_df[,"MT_max"], "% MT genes)"), fontface='bold', size = 20)
  print(plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
  dev.off()
  
  # replot feature scatter plots after filtering
  pdf(paste0(temp.plot.path, "featscat.postfilt_UMI.MT.pdf"),width=15,height=8)
  plot1 <- FeatureScatter(filt.data, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = condition, pt.size = 0.1)
  plot2 <- FeatureScatter(filt.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = condition, pt.size = 0.1)
  print(CombinePlots(plots = list(plot1, plot2)))
  dev.off()
  
  # save object in loop into the environment
  name<-paste0(condition, ".filt.data")
  filt.object.list[[name]]<-filt.data
}


# plot violin plots of all data after applying different filtering parameters - allows cross comparison of filtering levels
pdf(paste0(plot.path.filtcomp, "filter.comparison.vioplot.pdf"),width=30,height=15)
p1 = vio.plot.filt(meta.data = merged.data@meta.data, filt.param.df = filt_crit_df, y.dat = "nCount_RNA", y.lab = "UMI", x.lab = "Filtering Threshold")
p2 = vio.plot.filt(meta.data = merged.data@meta.data, filt.param.df = filt_crit_df, y.dat = "nFeature_RNA", y.lab = "Gene Count", x.lab = "Filtering Threshold")
p3 = vio.plot.filt(meta.data = merged.data@meta.data, filt.param.df = filt_crit_df, y.dat = "percent.mt", y.lab = "Percent Mitochondrial Genes", x.lab = "Filtering Threshold")
p = plot_grid(p1, p2, p3, ncol = 3)
title <- ggdraw() + draw_label(paste0("Filtering Threshold QC Comparison (all stages)"), fontface='bold', size = 20)
print(plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
dev.off()

# remove objects which are not needed to save RAM
rm(p, p1, p2, p3, temp.filt_crit_df, title, plot1, plot2, filt.data, g, filt.stats)


saveRDS(filt.object.list, "RDS.files/1_filter_QC/filt.object.list.RDS")



########### SCTransform ##########

                                              #
# run downstream analysis in order to get the HM plot to show which cells are lost during filtering #
#               then run downstream on only the dataset which will be kept                          #
                                              #

# run SCTransform (normalise the data)
# SCTransform replaces NormaliseData, ScaleData and FindVariableFeatures
# during this step upi can also remove confounding sources of variation, i.e. percent.mt
norm.unfilt <- SCTransform(object = filt.object.list[["unfilt.filt.data"]], vars.to.regress = "percent.mt", verbose = FALSE)
saveRDS(norm.unfilt, "RDS.files/1_filter_QC/norm.unfilt.RDS")
#norm.unfilt <- readRDS("RDS.files/1_filter_QC/norm.unfilt.RDS")

########### PCA, UMAP and clustering ##########

# Run PCA analysis on the each set of data
norm.unfilt <- RunPCA(object = norm.unfilt, verbose = FALSE)

# Seurat's clustering algorithm is based on principle components, so we need to ensure that only the informative PCs are kept!
# in order to determine PCA cutoffs you need to generate different types of plots with sufficient PCA dimensions in order to visualise the level at which including further PCA dimensions no longer provides further information - if 30 PCAs is not sufficient to capture this cutoff then increase this value
# dim heatmap is a heuristic method which will plot heatmaps of most extreme cells and their variable features for each principle component. This is a supervised analysis and is not necessarily the most appropriate to determine which PCs are the most informative
pdf(paste0(plot.path.filtcomp, "unfilt_dimHM.pdf"),width=15,height=25)
DimHeatmap(norm.unfilt, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

# another heuristic method is ElbowPlot which ranks PCs based on the % variance explained by each PC
pdf(paste0(plot.path.filtcomp, "unfilt_elbowplot.pdf"),width=12,height=10)
print(ElbowPlot(norm.unfilt, ndims = 40))
dev.off()

# FindNeighbors and FindClusters carry out the clustering
# FindNeighbors builds the SNN (shared nearest neighbour) graph
norm.unfilt <- FindNeighbors(norm.unfilt, dims = 1:15, verbose = FALSE)
      
# FindClusters runs community detection on the SNN graph - uses PCA data by default in order to identify clusters resolution arg can be altered to increase or reduce no. of clusters
norm.unfilt <- FindClusters(norm.unfilt, resolution = 0.5, verbose = FALSE)
      
# RunUMAP carries out non-linear dimensionality reduction for data visualisation
norm.unfilt <- RunUMAP(object = norm.unfilt, dims = 1:15, verbose = FALSE)

########### UMAP feature plots ##########

# UMAP plots of QC features to compare different filtering thresholds
p1 <- DimPlot(norm.unfilt, label = TRUE)
p2 <- DimPlot(norm.unfilt, group.by = "orig.ident", label = TRUE)
p3 <- FeaturePlot(norm.unfilt, features = "nFeature_RNA") + theme(plot.title = element_blank())
p4 <- FeaturePlot(norm.unfilt, features = "percent.mt") + theme(plot.title = element_blank())
p5 <- FeaturePlot(norm.unfilt, features = "nCount_RNA") + theme(plot.title = element_blank())
p6 <- DimPlot(norm.unfilt, group.by = "low")
p7 <- DimPlot(norm.unfilt, group.by = "med")
p8 <- DimPlot(norm.unfilt, group.by = "high")
lab = c("Clusters", "Developmental Stage", "Gene Count", "Percentage Mitochondrial Content", "UMI Count", "Low Filter Threshold",
        "Medium Filter Threshold", "High Filter Threshold")
p <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4, labels = lab, label_x = 0.1, hjust = 0, label_size = 17)
title <- ggdraw() + draw_label(paste0("UMAP plots of different QC metrics in the unfiltered normalised dataset"), fontface='bold', size = 20)
pdf(paste0(plot.path.filtcomp, "UMAP_filter_QC.pdf"),width=36,height=16)
print(plot_grid(title, p, ncol=1, rel_heights=c(0.2, 1)))
dev.off()
rm(p1, p2, p3, p4, p5, p6, p7, p8, p, lab, title)


# UMAP plot just showing QC features before filtering
p1 <- DimPlot(norm.unfilt, group.by = "orig.ident", label = F)
p2 <- FeaturePlot(norm.unfilt, features = "nFeature_RNA") + theme(plot.title = element_blank())
p3 <- FeaturePlot(norm.unfilt, features = "percent.mt", max.cutoff = 75) + theme(plot.title = element_blank())
p4 <- FeaturePlot(norm.unfilt, features = "nCount_RNA", max.cutoff = 75000) + theme(plot.title = element_blank())
lab = c("Developmental Stage", "Gene Count", "Percentage Mitochondrial Content", "UMI Count")
p <- plot_grid(p1, p2, p3, p4, ncol = 2, labels = lab, label_x = 0.1, hjust = 0, label_size = 17, label_fontface = "plain")
pdf(paste0(plot.path, "UMAP_QC.pdf"),width=16,height=16)
print(p)
dev.off()
rm(p1, p2, p3, p4, p, lab)

#################### DGEA ########################
# Differential expression test
markers <- FindAllMarkers(norm.unfilt, only.pos = T, logfc.threshold = 0.25, test.use = "negbinom")
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
top15 <- unique(subset(top15$gene, top15$gene %in% rownames(x = GetAssayData(object = norm.unfilt, slot = "scale.data"))))
clusters <- norm.unfilt[["seurat_clusters"]][order(norm.unfilt[["seurat_clusters"]]),, drop=FALSE]

# merge two dataframes together (clusters and filter condition)
HMcol <- merge(norm.unfilt[["low"]], clusters, by="row.names", all.x = T)
rownames(HMcol) <- HMcol$Row.names
HMcol$Row.names <- NULL
HMcol <- merge(norm.unfilt[["med"]], HMcol, by="row.names", all.x = T)
rownames(HMcol) <- HMcol$Row.names
HMcol$Row.names <- NULL
HMcol <- merge(norm.unfilt[["high"]], HMcol, by="row.names", all.x = T)
rownames(HMcol) <- HMcol$Row.names
HMcol$Row.names <- NULL
# add stage to metadata
HMcol <- merge(norm.unfilt[["orig.ident"]], HMcol, by="row.names", all.x = T)
rownames(HMcol) <- HMcol$Row.names
HMcol$Row.names <- NULL
# order cells by clusters
HMcol <- HMcol[order(HMcol$seurat_clusters),]
# Specify colors
ann_colors = list(
  low = c("Remaining Cells" = "#2e9112", "Filtered Cells"  = "#c80e0e") ,
  med = c("Remaining Cells" = "#2e9112", "Filtered Cells"  = "#c80e0e") ,
  high = c("Remaining Cells" = "#2e9112", "Filtered Cells"  = "#c80e0e"))

# make matrix from SCT scaled data - selecting only the genes from the topn object, and ordering the data according to topn and clusters which has been ordered above to match the DoHeatmap plot
data <- t(as.matrix(x = GetAssayData(object = norm.unfilt, assay = "SCT", slot = "scale.data")[top15, rownames(clusters), drop = FALSE]))
data <- data[rownames(HMcol),] 
data <- replace(data, data >= 2.5, 2.5)
data <- replace(data, data <= -2.5, -2.5)

# check all rownames in data and heatmap annotation dataframe (HMcol) are the same
setdiff(rownames(data),rownames(HMcol)) 

wid <- 60
ann_col <- HMcol[,1:4]

pdf(paste0(plot.path.filtcomp, "HM.top15.DE.PCA15.pdf"), width=28, height= 25)
pheatmap(t(data), color = colorRampPalette(c("magenta", "gray1", "yellow"))(100),
         cluster_rows = F, cluster_cols = F, show_colnames = F,
         annotation_col = ann_col, annotation_colors = ann_colors,
         fontsize = 22, fontsize_row = 12, gaps_col = cumsum(as.vector(table(clusters))),
         main = paste0("Heatmap of top15 DE genes per cluster in the unfiltered cell dataset"))
dev.off()


# Med filter is best

