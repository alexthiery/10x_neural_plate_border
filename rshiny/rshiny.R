library(shiny)
library(shinydashboard)
library(Seurat)
library(dplyr)

# dat <- readRDS('../shiny_test.RDS')
dat_list <- readRDS('./output/rshiny_input.RDS')

# ####################################################################
# # Initialise test list of diff stages
# start = c(1, 1:7 * 10 + 1)
# end = start + 9
# 
# dat_list2 = list()
# for(i in 1:8){
#   dat_list2[[i]] <- subset(dat, cells = start[i]:end[i])
# }
# 
# names(dat_list2) <- c('Full dataset', 'NPB subset', 'HH4', 'HH5', 'HH6', 'HH7', 'ss4', 'ss8')
# 
# 
# gene_ids = rownames(dat)
####################################################################


gene_ids <- lapply(dat_list, VariableFeatures) %>% Reduce(c, .) %>% unique()


group_options = c('Stage', 'Cell state', 'Clusters')

tab_home          <- tabItem(tabName = "home",
                             h2("Widgets tab content")
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
                                      radioButtons("subset_lineage_dynamics", "Select dataset to visualise", c('Full dataset', 'NPB subset'), inline = TRUE, selected = 'Full dataset')
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


server <- function(input, output){
  output$dimplot <- renderPlot(DimPlot(dat_list[[input$subset_featureplots]], group.by = input$group))
  output$featureplot <- renderPlot(FeaturePlot(dat_list[[input$subset_featureplots]], features = input$gene_id))
}

shinyApp(ui, server)


