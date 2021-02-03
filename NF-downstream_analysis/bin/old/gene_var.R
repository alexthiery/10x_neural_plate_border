#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer",
  'custom_functions', 'm', 2, "character",
  'networkGenes', 'd', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set run location
if(length(commandArgs(trailingOnly = TRUE)) == 0){
  cat('No command line arguments provided, user defaults paths are set for running interactively in Rstudio on docker\n')
  opt$runtype = "user"
} else {
  if(is.null(opt$runtype)){
    stop("--runtype must be either 'user' or 'nextflow'")
  }
  if(tolower(opt$runtype) != "user" & tolower(opt$runtype) != "nextflow"){
    stop("--runtype must be either 'user' or 'nextflow'")
  }
  if(tolower(opt$runtype) == "nextflow"){
    if(is.null(opt$custom_functions) | opt$custom_functions == "null"){
      stop("--custom_functions path must be specified in process params config")
    }
  }
}

# Set paths and load data
{
  if (opt$runtype == "user"){
    sapply(list.files('./bin/R/custom_functions/', full.names = T), source)
    plot.path = "./output/plots/seurat_STACAS/"
    rds.path = "./output/RDS.files/seurat_STACAS/"
    antler.path = "./output/antler/"
    data.path = "./output/RDS.files/seurat_STACAS/"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    sapply(list.files(opt$custom_functions, full.names = T), source)
    plot.path = "./plots/"
    rds.path = "./RDS.files/"
    antler.path = "./antler/"
    data.path = "./"
    
    ncores = opt$cores
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  
  dir.create(plot.path, recursive = T)
  dir.create(rds.path, recursive = T)
  dir.create(antler.path, recursive = T)
  
  # Load packages - packages are stored within renv in the repository
  reticulate::use_python('/usr/bin/python3.7')
  library(Seurat)
  
  library(future)
  library(dplyr)
  library(Antler)
  library(cowplot)
  library(clustree)
  library(gridExtra)
  library(grid)
  library(pheatmap)
  library(biomaRt)
  library(RColorBrewer)
}

# read in seurat data
seurat_out <- readRDS(paste0(data.path, "seurat_out_hh4filt.RDS"))


plot_dat <- GetAssayData(subset(seurat_out), assay = "RNA", slot = "data")
# remove genes with median expression of 0
plot_dat <- plot_dat[apply(plot_dat, 1, function(x) median(x) != 0),]


temp <- lapply(c("hh5", "hh6", "hh7", "ss4", "ss8"), function(x) {
  dat <- plot_dat[,grepl(x, colnames(plot_dat))]
  data.frame(gene_var = apply(dat, 1, var), median_expression = apply(dat, 1, median), cluster = x)
})

temp <- bind_rows(temp)






dat <- dat[apply(dat, 1, function(x) median(x) != 0),]
y = apply(dat, 1, var)
x = apply(dat, 1, median)
hh5 <- as.data.frame(cbind(y, x))

dat <- GetAssayData(subset(seurat_out, cells = rownames(seurat_out@meta.data)[seurat_out@meta.data$stage == "hh6"]), assay = "RNA", slot = "data")
dat <- dat[apply(dat, 1, function(x) median(x) != 0),]
y = apply(dat, 1, var)
x = apply(dat, 1, median)
hh6 <- as.data.frame(cbind(y, x))

dat <- GetAssayData(subset(seurat_out, cells = rownames(seurat_out@meta.data)[seurat_out@meta.data$stage == "hh7"]), assay = "RNA", slot = "data")
dat <- dat[apply(dat, 1, function(x) median(x) != 0),]
y = apply(dat, 1, var)
x = apply(dat, 1, median)
hh7 <- as.data.frame(cbind(y, x))

dat <- GetAssayData(subset(seurat_out, cells = rownames(seurat_out@meta.data)[seurat_out@meta.data$stage == "ss4"]), assay = "RNA", slot = "data")
dat <- dat[apply(dat, 1, function(x) median(x) != 0),]
y = apply(dat, 1, var)
x = apply(dat, 1, median)
ss4 <- as.data.frame(cbind(y, x))

dat <- GetAssayData(subset(seurat_out, cells = rownames(seurat_out@meta.data)[seurat_out@meta.data$stage == "ss8"]), assay = "RNA", slot = "data")
dat <- dat[apply(dat, 1, function(x) median(x) != 0),]
y = apply(dat, 1, var)
x = apply(dat, 1, median)
ss8 <- as.data.frame(cbind(y, x))

AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), cluster = x)
  }))
}

temp <- AppendMe(c("hh5", "hh6", "hh7", "ss4", "ss8"))

ggplot(temp, aes(x = median_expression, y = gene_var, group = cluster, color = cluster)) +
  geom_smooth(aes(group = cluster), method="auto", se=FALSE, fullrange=FALSE, level=0.95) +
  ylim(c(0,0.3)) +
  theme_classic() +
  ylab("gene variance") +
  xlab("median gene expression")



ggplot(temp, aes(x=cluster, y=y, colour = cluster)) +
  geom_violin() +
  ylim(c(0,0.3)) +
  theme_classic() +
  ylab("gene variance")

ggplot(temp, aes(x=cluster, y=gene_var, colour = cluster)) +
  geom_boxplot() +
  ylim(c(0,0.3)) +
  theme_classic() +
  ylab("gene variance")


x = lm(formula = temp$gene_var ~ temp$cluster)
summary(x)

TukeyHSD(x)
TukeyHSD(aov(temp$gene_var ~ temp$cluster))

pairwise.wilcox.test(temp$gene_var, temp$cluster, p.adjust.method = "bonf")

wilcox.test(rnorm(10), rnorm(10, 2), conf.int = TRUE)

wilcox.test(x = temp$x)

