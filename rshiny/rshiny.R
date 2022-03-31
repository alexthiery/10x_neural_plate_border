library(shiny)
library(shinydashboard)
library(Seurat)
library(tidyverse)
library(viridis)

# dat <- readRDS('../shiny_test.RDS')
dat_list <- readRDS('./output/rshiny_input.RDS')

# Function for generating dataframe for gene expression dynamics
prepare_dynamics <- function(gene, seurat_object, assay = 'RNA', slot = 'data', lineages){
  assay_dat <- as.data.frame(GetAssayData(seurat_object, assay = 'RNA', slot = 'data')[gene,, drop=FALSE])
  
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
    mutate(lineage = unlist(strsplit(lineage, '_'))[2])
  
  return(plot_data)
}




scvelo_features <- c('latent_time', 'lineage_neural_probability', 'lineage_NC_probability', 'lineage_placodal_probability')



# Get shared var features for optional input
gene_ids <- lapply(dat_list, VariableFeatures) %>% Reduce(c, .) %>% unique()

# Set options for viewing DimPlots
group_options = c('Stage', 'Cell state', 'Clusters')

tab_home          <- tabItem(tabName = "home",
                             h2("Home page <TO FILL>")
)

tab_feature_plots <- tabItem(tabName = "featureplots",
                             fluidRow(
                               column(11,
                                      radioButtons("subset_featureplots", "Select dataset to visualise", names(dat_list), inline = TRUE, selected = 'Full data')
                               )
                             ),
                             fluidRow(
                               column(5,
                                      box(
                                        selectInput("group", "Group by", group_options, width = "200", selected = 'Cell state'),
                                        plotOutput("dimplot"),
                                        width = 12
                                      )
                               ),
                               column(5, offset = 1,
                                      box(
                                        selectInput("gene_id", "Select Gene", gene_ids, width = "200", selected = 'PAX7'),
                                        plotOutput("featureplot"),
                                        width = 12
                                      )
                               )
                             )
)

tab_lineage_dynamics <- tabItem(tabName = "lineage_dynamics",
                                fluidRow(
                                  column(11,
                                         radioButtons("select_lineage_dataset", "Select dataset to visualise", c('Full data', 'NPB subset'), inline = TRUE, selected = 'Full data')
                                  ),
                                  column(6,
                                         box(
                                           selectInput("feature_lineage_dynamics", "Select Feature"),
                                           plotOutput("scvelo_umaps"),
                                           width = 12
                                         )
                                  ),
                                  column(6,
                                         box(
                                           checkboxGroupInput("select_lineage", "Select lineages of interest", inline = TRUE),
                                           selectInput("gene_id_lineage_dynamics", "Select Gene", gene_ids, width = "200", selected = 'PAX7'),
                                           plotOutput("lineage_dynamics"),
                                           width = 12
                                         )
                                  )
                                )
)


ui <- dashboardPage(
  dashboardHeader(),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home"),
      menuItem("Feature Plots", tabName = "featureplots"),
      menuItem("Lineage Dynamics", tabName = "lineage_dynamics")
    )
  ),
  
  dashboardBody(
    tabItems(
      tab_home,
      tab_feature_plots,
      tab_lineage_dynamics
    )
  )
)



server <- function(input, output, session){
  
  ####################################################################
  # Generate seurat umaps (cell state, stage, clusters)
  output$dimplot <- renderPlot(DimPlot(dat_list[[input$subset_featureplots]], group.by = input$group))
  
  
  subset_dataset <- reactive({dat_list[[input$select_lineage_dataset]]})
  ####################################################################
  # Generate gene feature plots
  output$featureplot <- renderPlot(FeaturePlot(dat_list[[input$subset_featureplots]], features = input$gene_id))

  
  ####################################################################
  # Generate scvelo latent time/lineage probability plots 
  # scvelo_umap <- reactive({FeaturePlot(dat_list[[input$select_lineage_dataset]], features = 'lineage_neural_probability', pt.size = 1, order = FALSE) +
  #                            scale_colour_viridis(option="viridis")})
  
  
  # Set lineage options based on dataset selected
  observe({
    
    if(input$select_lineage_dataset == 'Full data'){
      features = c('latent_time', 'lineage_neural_probability', 'lineage_NC_probability', 'lineage_placodal_probability')
    } else if (input$select_lineage_dataset == 'NPB subset'){
      features = c('latent_time', 'lineage_NC_probability', 'lineage_placodal_probability')
    } else {
      features = character(0)
    }
    
    updateSelectInput(session, "feature_lineage_dynamics",
                      choices = features,
                      selected = 'latent_time'
    )
    
  })
  

  output$scvelo_umaps <- renderPlot(FeaturePlot(subset_dataset(), features = input$feature_lineage_dynamics, pt.size = 1, order = FALSE) +
                                      scale_colour_viridis(option="viridis"))

  FeaturePlot(dat_list$`Full data`, features = 'latent_time', pt.size = 2) +
    scale_colour_viridis(option="plasma")


  ####################################################################
  # Generate lineage dynamics plots
  
  # Set lineage options based on dataset selected
  observe({
    
    if(input$select_lineage_dataset == 'Full data'){
      lineages = c('Neural crest', 'Neural', 'Placodal')
    } else if (input$select_lineage_dataset == 'NPB subset'){
      lineages = c('Neural crest', 'Placodal')
    } else {
      lineages = character(0)
    }

    updateCheckboxGroupInput(session, "select_lineage",
                             choices = lineages,
                             selected = 'Neural crest'
    )
  })
  
  plot_data <- reactive({prepare_dynamics(input$gene_id_lineage_dynamics, subset_dataset(), lineages = input$select_lineage)})

  output$lineage_dynamics <- renderPlot(ggplot(plot_data(), aes(x = latent_time, y = scaled_expression, colour = lineage)) +
                                          geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, group=lineage)) +
                                          xlab("Latent time") + ylab("Scaled expression") +
                                          theme_classic())
}

shinyApp(ui, server)

# library(viridis)
# 
# 
# subset_lineage_probabilities <- c('Neural' = 'lineage_neural_probability', 'Neural crest' = 'lineage_NC_probability', 'Placodal' = 'lineage_placodal_probability')
# 

