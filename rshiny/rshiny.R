library(shiny)
library(shinydashboard)
library(Seurat)
library(tidyverse)
library(viridis)
library(mgcv)


# dat <- readRDS('../shiny_test.RDS')
dat_list <- readRDS('./output/rshiny_input.RDS')






calc_latent_time_cutoff <- function(latent_time, lineage_probability, top_frac = 0.2, return = 'intercept'){
  data = data.frame(latent_time, lineage_probability)
  
  model_data <- data %>%
    filter(lineage_probability > 0 & lineage_probability < 1) %>%
    filter(lineage_probability > quantile(lineage_probability, 1-top_frac)) %>%
    # If cells remaining in top frac have less than 0.5 probability of giving rise to the lineage then remove (required for short lineages)
    filter(lineage_probability > 0.5)
  
  # Fit model to top frac of data
  fit = lm(lineage_probability ~ latent_time, data = model_data)
  
  # Inverse equation to find X for Y = 1
  x <- (1-coef(fit)[1])/coef(fit)[2]
  
  # Identify max latent time value based on cells which have reached lineage probability == 1
  max_latent_time <- data %>% filter(lineage_probability == 1) %>% filter(latent_time == max(latent_time)) %>% dplyr::pull(latent_time)
  
  if(return == 'plot'){
    p = ggplot(data, aes(latent_time, lineage_probability)) + 
      geom_point(size = 0.1) +
      geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2])
    
    return(p)
  }else if(return == 'intercept'){
    if(x > max_latent_time){
      cat('Predicted latent time is later than any cell observed that has reached full lineage absorbtion. Max latent time is set based on oldest cell at lineage_probability == 1')
      return(max_latent_time)
    }else{
      return(unname(x))
    }
  }else{
    stop('return must be one of intercept or plot')
  }
}

calc_gams <- function(dynamics_data){
  # GAMs
  gams <- dynamics_data %>%
    group_by(lineage) %>%
    do(gams = gam(scaled_expression ~ s(latent_time, bs = "cs", k=5), weights = lineage_probability, data = .))
  
  # Add module column and max latent time to gam data
  extra_dat <- dynamics_data %>%
    ungroup() %>%
    dplyr::select(gene, lineage, max_latent_time) %>%
    distinct()
  
  gams <- inner_join(gams, extra_dat)
  
  # Generate predicted values for each gam in tidy df -> output in long format
  plot_data <- data.frame()
  for(row in 1:nrow(gams)){
    # Generate latent time values to predict gams on -> use max latent_time calculated per lineage
    pdat <- tibble(latent_time = seq(0, gams[[row, 'max_latent_time']], length = 100))
    new_data <- predict.gam(gams[[row,'gams']][[1]], newdata = pdat, se=TRUE)
    
    plot_data <- rbind(plot_data, data.frame(gene = gams[[row, 'gene']],
                                             lineage = gams[[row, 'lineage']],
                                             scaled_expression = new_data[['fit']],
                                             se = new_data[['se.fit']],
                                             pdat))
  }
  return(plot_data)
}



# Function for generating dataframe for gene expression dynamics
get_lineage_expression_data <- function(gene, seurat_object, assay = 'RNA', slot = 'data', lineages, lineage_cutoffs = NULL){
  
  if(!all(names(lineage_cutoffs) %in% lineages)){
    stop('lineage_cutoffs names do not match selected lineages')
  }
  
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
    mutate(lineage = names(subset_lineage_probabilities)[subset_lineage_probabilities %in% lineage])
  
  # Filter expression values which are above the predicted lineage cutoffs 
  if(!is.null(lineage_cutoffs)){
    plot_data <- plot_data %>%
      mutate(max_latent_time = lineage_cutoffs[names(lineage_cutoffs) %in% lineage]) %>%
      filter(latent_time < max_latent_time)
  }
  
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
                                           # verbatimTextOutput("lineage_dynamics"),
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

lineage_map = c('Neural' = 'lineage_neural_probability', 'Neural crest' = 'lineage_NC_probability', 'Placodal' = 'lineage_placodal_probability')





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
  
  linegae_expression_data <- reactive(get_lineage_expression_data(input$gene_id_lineage_dynamics, subset_dataset(), lineages = input$select_lineage, lineage_cutoffs = max_cutoffs()))

  plot_data <- reactive(calc_gams(linegae_expression_data()))
  
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
  
}

shinyApp(ui, server)

