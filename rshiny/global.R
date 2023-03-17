library(shiny)
library(bs4Dash)
library(Seurat)
library(tidyverse)
library(viridis)
library(mgcv)
library(patchwork)
library(shinycssloaders)

source('./custom_functions.R')

options(scipen = 1)
options(digits = 2)

# dat <- readRDS('../shiny_test.RDS')
dat_list <- readRDS('../output/rshiny_input.RDS')



# Set parameters
lineage_map = c('Neural' = 'lineage_neural_probability', 'Neural crest' = 'lineage_NC_probability', 'Placodal' = 'lineage_placodal_probability')
scvelo_features = c('latent_time', 'lineage_neural_probability', 'lineage_NC_probability', 'lineage_placodal_probability')
lineage_names = c('lineage_NC_probability' = 'Neural crest', 'lineage_neural_probability' = 'Neural', 'lineage_placodal_probability' = 'Placodal')


# Get shared var features for optional input
gene_ids <- lapply(dat_list, VariableFeatures) %>% Reduce(c, .) %>% unique() %>% sort()

# Set options for viewing DimPlots
group_options = c('Stage', 'Cell state', 'Clusters')

my_theme <- theme(axis.text=element_text(size=14),
                  axis.title=element_text(size=16))