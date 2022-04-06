library(shiny)
library(shinydashboard)


# Set parameters
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
                               column(12,
                                      radioButtons("subset_featureplots", "Select dataset to visualise", names(dat_list), inline = TRUE, selected = 'Full data', width = '500')
                               )
                             ),
                             fluidRow(
                               column(5,
                                      box(
                                        selectInput("group", "Group by", group_options, width = "250", selected = 'Cell state'),
                                        plotOutput("dimplot"),
                                        width = 12,
                                        style='height:35vw'
                                      )
                               ),
                               column(5,
                                      box(
                                        selectInput("gene_id", "Select Gene", gene_ids, width = "250", selected = 'PAX7'),
                                        plotOutput("featureplot"),
                                        width = 12,
                                        style='height:35vw'
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
                               column(
                                 4,
                                 fluidRow(
                                   box(
                                     fluidRow(
                                       column(5,
                                              selectInput("coexpression_gene_id_1", "Select Gene 1", gene_ids, width = "200", selected = 'PAX7'),
                                       ),
                                       column(5,
                                              offset = 1,
                                              textInput("gene_1_col", "Colour", value = 'magenta'),
                                              ),
                                       ),
                                     
                                     fluidRow(
                                       column(12,
                                              # offset = 1,
                                         plotOutput("coexpression_featureplot_1"),
                                         )
                                       ),
                                     
                                     width = 12, style='height:23vw'
                                     )
                                  ),
                                 fluidRow(
                                   box(
                                     fluidRow(
                                       column(5,
                                              selectInput("coexpression_gene_id_2", "Select Gene 2", gene_ids, width = "200", selected = 'SIX1'),
                                       ),
                                       column(5,
                                              offset = 1,
                                              textInput("gene_2_col", "Colour", value = 'green'),
                                       ),
                                     ),
                                     
                                     fluidRow(
                                       column(12,
                                              # offset = 1,
                                              plotOutput("coexpression_featureplot_2", width = '80%')
                                              ),
                                       tags$style(
                                         "coexpression_featureplot_2 {border-color: red;}")
                                       ),
                                     width = 12, style='height:23vw'
                                     )
                                   )
                                 ),
                               
                               column(
                                 6,
                                  box(
                                    fluidRow(
                                      column(
                                        6, offset = 3,
                                        h3('Co-expression of gene 1 and gene 2')
                                      )
                                    ),
                                    br(),
                                    fluidRow(
                                      column(
                                        6, offset = 1,
                                        sliderInput("threshold", "Threshold", value = 0.1, min = 0, max = 1, step = 0.05, width = '300'),
                                      )
                                    ),
                                    br(),
                                    fluidRow(
                                      column(
                                        10, offset = 1,
                                        plotOutput("coexpression_umap")
                                      )
                                    ),
                                    width = 12, style='height:46vw'
                                    )
                                  )
                               )
)


ui <- dashboardPage(
  dashboardHeader(title = "Custom font"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home"),
      menuItem("Feature Plots", tabName = "featureplots"),
      menuItem("Lineage Dynamics", tabName = "lineage_dynamics"),
      menuItem("UMAP co-expression", tabName = "coexpression_umaps")
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "www/custom.css")
    ),
    tabItems(
      tab_home,
      tab_feature_plots,
      tab_lineage_dynamics,
      tab_coexpression_umaps
    )
  )
)








