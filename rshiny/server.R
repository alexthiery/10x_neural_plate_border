
server <- function(input, output, session){
  
  ####################################################################
  # Generate seurat umaps (cell state, stage, clusters)
  output$dimplot <- renderPlot(DimPlot(dat_list[[input$subset_featureplots]], group.by = input$group),
                               height = function() {session$clientData$output_dimplot_width * 0.8})
  
  
  subset_dataset <- reactive({dat_list[[input$select_lineage_dataset]]})
  ####################################################################
  # Generate gene feature plots
  output$featureplot <- renderPlot(FeaturePlot(dat_list[[input$subset_featureplots]], features = input$gene_id),
                                   height = function() {session$clientData$output_featureplot_width * 0.8})

  
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
  output$coexpression_featureplot_1 <- renderPlot(FeaturePlot(dat_list[[input$coexpression_subset_featureplots]],
                                                              features = input$coexpression_gene_id_1,
                                                              pt.size = 0.75) +
                                                    scale_colour_gradientn(colours = colorRampPalette(c("gray85", input$gene_1_col))(5)),
                                                  height = function() {session$clientData$output_coexpression_featureplot_1_width * 0.75}
  )

  output$coexpression_featureplot_2 <- renderPlot(FeaturePlot(dat_list[[input$coexpression_subset_featureplots]],
                                                              features = input$coexpression_gene_id_2,
                                                              pt.size = 0.75) +
                                                    scale_colour_gradientn(colours = colorRampPalette(c("gray85", input$gene_2_col))(5)),
                                                  height = function() {session$clientData$output_coexpression_featureplot_2_width * 0.75}
  )
  
  output$coexpression_umap <- renderPlot(coexpression_umap(dat_list[[input$coexpression_subset_featureplots]], gene_1 = input$coexpression_gene_id_1,
                                                                   gene_2 = input$coexpression_gene_id_2, col.threshold = input$threshold, two.colors = c(input$gene_1_col, input$gene_2_col),
                                                                   negative.color = 'gray90', limit = 1, highlight_cell_size = 2),
                                         height = function() {session$clientData$output_coexpression_umap_width * 0.7}
                                         )  

}