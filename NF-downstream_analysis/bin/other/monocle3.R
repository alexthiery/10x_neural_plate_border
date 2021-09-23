library(Seurat)
library(scHelper)
library(monocle3)
library(igraph)


seurat_data <- readRDS('./output/NF-downstream_analysis_stacas/filtered_seurat/seurat/state_classification/rds_files/contamination_cell_state_classification.RDS')



gene_annotation = seurat_data@assays$RNA@meta.features
gene_annotation$gene_short_name = rownames(gene_annotation)

cds <- new_cell_data_set(seurat_data@assays$RNA@counts,
                         cell_metadata = seurat_data@meta.data,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 100)


reducedDims(cds) <- list(UMAP = seurat_data@reductions$umap@cell.embeddings)

cds <- cluster_cells(cds)


colData(cds)$scHelper_cell_type = factor(seurat_data@meta.data$scHelper_cell_type)
cds <- learn_graph(cds, use_partition = T)


plot_cells(cds,
           color_cells_by = "scHelper_cell_type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=T,
           cell_size = 0.5,
           label_cell_groups=F)+NoLegend()


cds <- order_cells(cds)


plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
