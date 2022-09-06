library(shinythemes)
tab_home          <- tabItem(tabName = "home",
                             fluidRow(
                               column(8, offset = 2,
                                      includeMarkdown("home.md")
                               )
                             )
)

tab_feature_plots <- tabItem(tabName = "featureplots",
                             fluidRow(
                               column(12,
                                      radioButtons("subset_featureplots", "Select dataset to visualise", names(dat_list), inline = TRUE, selected = 'Full data', width = '800')
                               )
                             ),
                             fluidRow(
                               column(5,
                                      box(
                                        selectizeInput("group", "Group by", group_options, width = "250", selected = 'Cell state'),
                                        plotOutput("dimplot"),
                                        width = 12,
                                        style='height:35vw'
                                      )
                               ),
                               column(5,
                                      box(
                                        selectizeInput("gene_id", "Select Gene", choices = NULL, width = "250"),
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
                                           fluidRow(
                                             column(5,
                                                    selectizeInput("feature_lineage_dynamics", "Select Feature", choices = NULL)
                                             )
                                           ),
                                           plotOutput("scvelo_umaps"),
                                           width = 12,
                                           style='height:35vw'
                                         )
                                  ),
                                  column(6,
                                         box(
                                           fluidRow(
                                             column(4,
                                                    selectizeInput("gene_id_lineage_dynamics", "Select Gene", choices = NULL, width = "200"),
                                             ),
                                             column(6,
                                                    offset = 2,
                                                    checkboxGroupInput("select_lineage", "Select lineages of interest", inline = TRUE),
                                             )
                                           ),
                                           br(),
                                           br(),
                                           fluidRow(
                                             column(12,
                                                    plotOutput("lineage_dynamics")
                                             )
                                           ),
                                           width = 12,
                                           style='height:35vw'
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
                                                   selectizeInput("coexpression_gene_id_1", "Select Gene 1", choices = NULL, width = "200"),
                                            ),
                                            column(5,
                                                   offset = 1,
                                                   textInput("gene_1_col", "Colour", value = 'magenta'),
                                            ),
                                          ),
                                          
                                          fluidRow(
                                            column(12,
                                                   # offset = 1,
                                                   plotOutput("coexpression_featureplot_1", width = '90%', height = '90%')
                                            )
                                          ),
                                          
                                          width = 12, style='height:23vw'
                                        )
                                      ),
                                      fluidRow(
                                        box(
                                          fluidRow(
                                            column(5,
                                                   selectizeInput("coexpression_gene_id_2", "Select Gene 2", choices = NULL, width = "200"),
                                            ),
                                            column(5,
                                                   offset = 1,
                                                   textInput("gene_2_col", "Colour", value = 'green'),
                                            ),
                                          ),
                                          
                                          fluidRow(
                                            column(12,
                                                   # offset = 1,
                                                   plotOutput("coexpression_featureplot_2", width = '90%', height = '90%')
                                            )
                                          ),
                                          
                                          width = 12, style='height:23vw'
                                        )
                                      )
                                    ),
                                    
                                    column(
                                      6,
                                      offset = 1,
                                      box(
                                        fluidRow(
                                          column(
                                            8, offset = 2,
                                            uiOutput('coexpression_header')
                                          )
                                        ),
                                        br(),
                                        br(),
                                        fluidRow(
                                          column(
                                            6, offset = 1,
                                            sliderInput("threshold", "Threshold", value = 0.1, min = 0, max = 1, step = 0.05, width = '300'),
                                          )
                                        ),
                                        br(),
                                        br(),
                                        fluidRow(
                                          column(
                                            12,
                                            # offset = 1,
                                            plotOutput("coexpression_umap")
                                          )
                                        ),
                                        width = 12, style='height:46vw'
                                      )
                                    )
                                  )
)

tab_dea <- tabItem(tabName = "dea",
                   fluidRow(
                     column(5,
                            radioButtons("subset_dea", "Select dataset to visualise", names(dat_list), inline = TRUE, selected = 'Full data', width = '800')
                     )
                   ),
                   br(),
                   fluidRow(
                     column(2,
                            selectInput("group_dea", "Group by", group_options, selected = 'Cell state', width = "250"),
                            selectInput('select_a', 'Select A', choices = NULL, multiple = TRUE, width = "250"),
                            selectInput("select_b", "Select B", choices = NULL, multiple = TRUE, width = "250"),
                            textInput("padj_threshold", "Adjusted p-value threshold", value = '0.01', width = "250"),
                            textInput("fc_threshold", "log2 FC threshold", value = '0.25', width = "250"),
                            fluidRow(column(6, actionButton("run_dea", "Run DEA")), column(6, uiOutput("download"))),
                            
                            align = 'left'
                     ),
                     column(6, offset = 1,
                            div(withSpinner(DT::dataTableOutput("dea_table"), color="#0dc5c1"), style = "font-size: 90%")
                     )
                   ),
)

ui <- dashboardPage(
  fullscreen = TRUE,
  header = dashboardHeader(
    title = dashboardBrand(
      title = "10x Neural Plate Border",
      href = "https://github.com/alexthiery/10x_neural_plate_border"
    )
    # skin = "light",
    # status = "white",
    # border = TRUE,
    # sidebarIcon = icon("bars"),
    
  ),
  
  
  
  dashboardSidebar(
    # tags$style("@import url(https://use.fontawesome.com/releases/v5.7.2/css/all.css);"),
    
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon('home')),
      menuItem("Feature Plots", tabName = "featureplots", icon = icon("border-none")),
      menuItem("Lineage Dynamics", tabName = "lineage_dynamics", icon = icon('chart-line')),
      menuItem("UMAP co-expression", tabName = "coexpression_umaps", icon = icon("braille")),
      menuItem("Differential expression", tabName = "dea", icon = icon("arrows-alt"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),
    tabItems(
      tab_home,
      tab_feature_plots,
      tab_lineage_dynamics,
      tab_coexpression_umaps,
      tab_dea
    )
  )
)

