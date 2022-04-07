library(shiny)
library(shinydashboard)
library(Seurat)
library(tidyverse)
library(viridis)
library(mgcv)

source('./custom_functions.R')

# dat <- readRDS('../shiny_test.RDS')
dat_list <- readRDS('../output/rshiny_input.RDS')



# Set parameters
lineage_map = c('Neural' = 'lineage_neural_probability', 'Neural crest' = 'lineage_NC_probability', 'Placodal' = 'lineage_placodal_probability')
scvelo_features = c('latent_time', 'lineage_neural_probability', 'lineage_NC_probability', 'lineage_placodal_probability')



# Get shared var features for optional input
gene_ids <- lapply(dat_list, VariableFeatures) %>% Reduce(c, .) %>% unique()

# Set options for viewing DimPlots
group_options = c('Stage', 'Cell state', 'Clusters')