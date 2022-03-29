library(shiny)
# install.packages("shinydashboard")
library(shinydashboard)
library(Seurat)
library(dplyr)

dat <- readRDS('../shiny_test.RDS')

# Slim seurat object to make it run faster
# Remove genes which have 0 expression in all cells


####################################################################
# Initialise test list of diff stages
start = c(1, 1:7 * 10 + 1)
end = start + 9

dat_list = list()
for(i in 1:8){
  dat_list[[i]] <- subset(dat, cells = start[i]:end[i])
}

names(dat_list) <- c('Full dataset', 'NPB subset', 'HH4', 'HH5', 'HH6', 'HH7', 'ss4', 'ss8')


gene_ids = rownames(dat)
####################################################################




# rename metadata columns
dat_list <- lapply(dat_list, function(x) {x@meta.data <- x@meta.data %>%
         rename(
           Stage = stage,
           'Cell state' = scHelper_cell_type,
           Clusters = seurat_clusters
         )
       return(x)
       })



group_options = c('Stage', 'Cell state', 'Clusters')

subset_options = c('Full dataset', 'NPB subset', 'HH4', 'HH5', 'HH6', 'HH7', 'ss4', 'ss8')


tab_home          <- tabItem(tabName = "home",
                             h2("Widgets tab content")
                             )

tab_feature_plots <- tabItem(tabName = "featureplots",
                             fluidRow(
                               column(11,
                                radioButtons("subset", "Select dataset to visualise", subset_options, inline = TRUE, selected = 'Full dataset')
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
                                      radioButtons("subset", "Select dataset to visualise", c('Full dataset', 'NPB subset'), inline = TRUE, selected = 'Full dataset')
                               )
                             )
                             # fluidRow(
                             #   column(5,
                             #          box(
                             #            selectInput("group", "Group by", group_options, width = "200", selected = 'Cell state'),
                             #            plotOutput("dimplot"),
                             #            width = 12
                             #          )
                             #   ),
                             #   column(5, offset = 1,
                             #          box(
                             #            selectInput("gene_id", "Select Gene", gene_ids, width = "200", selected = 'PAX7'),
                             #            plotOutput("featureplot"),
                             #            width = 12
                             #          )
                             #   )
                             # )
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
  output$dimplot <- renderPlot(DimPlot(dat_list[[input$subset]], group.by = input$group))
  output$featureplot <- renderPlot(FeaturePlot(dat_list[[input$subset]], features = input$gene_id))
}

shinyApp(ui, server)

