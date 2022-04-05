library(shiny)
library(shinydashboard)
library(Seurat)
library(tidyverse)
library(viridis)
library(mgcv)

source('./custom_functions.R')

# dat <- readRDS('../shiny_test.RDS')
dat_list <- readRDS('./output/rshiny_input.RDS')



# Set parameters
lineage_map = c('Neural' = 'lineage_neural_probability', 'Neural crest' = 'lineage_NC_probability', 'Placodal' = 'lineage_placodal_probability')
scvelo_features = c('latent_time', 'lineage_neural_probability', 'lineage_NC_probability', 'lineage_placodal_probability')



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
                                           selectInput("feature_lineage_dynamics", "Select Feature", choices = scvelo_features),
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

tab_coexpression_umaps <- tabItem(tabName = "coexpression_umaps",
                             fluidRow(
                               column(11,
                                      radioButtons("coexpression_subset_featureplots", "Select dataset to visualise", names(dat_list), inline = TRUE, selected = 'Full data')
                               )
                             ),
                             fluidRow(
                               column(5, offset = 1,
                                      box(
                                        selectInput("coexpression_gene_id_1", "Select Gene 1", gene_ids, width = "200", selected = 'PAX7'),
                                        plotOutput("coexpression_featureplot_1"),
                                        width = 12
                                      )
                               ),
                               column(5, offset = 1,
                                      box(
                                        selectInput("coexpression_gene_id_2", "Select Gene 2", gene_ids, width = "200", selected = 'SIX1'),
                                        plotOutput("coexpression_featureplot_2"),
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
      menuItem("Lineage Dynamics", tabName = "lineage_dynamics"),
      menuItem("UMAP co-expression", tabName = "coexpression_umaps")
    )
  ),
  
  dashboardBody(
    tabItems(
      tab_home,
      tab_feature_plots,
      tab_lineage_dynamics,
      tab_coexpression_umaps
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
  
  output$scvelo_umaps <- renderPlot(
    if(input$feature_lineage_dynamics == 'latent_time'){
      FeaturePlot(subset_dataset(), features = input$feature_lineage_dynamics, pt.size = 1, order = FALSE) +
        scale_colour_viridis(option="magma")
    } else {
      FeaturePlot(subset_dataset(), features = input$feature_lineage_dynamics, pt.size = 1, order = FALSE) +
        scale_colour_viridis(option="viridis")
      }
    )


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


  subset_lineage_map <- reactive({lineage_map[input$select_lineage]})

  max_cutoffs <- reactive({
    unlist(
      setNames(
          lapply(names(subset_lineage_map()), function(x) {
          calc_latent_time_cutoff(
            latent_time = subset_dataset()@meta.data[['latent_time']],
            lineage_probability = subset_dataset()@meta.data[[subset_lineage_map()[[x]]]]
            )
          }),
          names(subset_lineage_map())
      )
    )
  })
  
  plot_data <- reactive(lineage_gam(input$gene_id_lineage_dynamics, subset_dataset(), slot = 'scale.data', lineages = input$select_lineage, lineage_cutoffs = max_cutoffs()))
  
  output$lineage_dynamics <- renderPlot(ggplot(plot_data(), aes(x = latent_time, y = scaled_expression, colour = lineage, fill = lineage)) +
    geom_line(size = 2) +
    geom_ribbon(aes(ymin=scaled_expression-se, ymax=scaled_expression+se), alpha = .3, colour = NA) +
    labs(x = "Latent time", y = "Scaled Expression") +
    theme_classic() +
    theme(legend.key.size = unit(1,"cm"), 
          legend.title = element_blank(),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),  
          axis.title.x = element_text(size = 24),
          axis.title.y = element_text(size = 24),
          legend.text=element_text(size=20))
  )
  

  ####################################################################
  # Co-expression feature plots
  output$coexpression_featureplot_1 <- renderPlot(FeaturePlot(dat_list[[input$coexpression_subset_featureplots]], features = input$coexpression_gene_id_1))

  output$coexpression_featureplot_2 <- renderPlot(FeaturePlot(dat_list[[input$coexpression_subset_featureplots]], features = input$coexpression_gene_id_2))

}

shinyApp(ui, server)

