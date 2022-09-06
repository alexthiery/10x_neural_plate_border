
server <- function(input, output, session){
  
  ####################################################################
  # Generate seurat umaps (cell state, stage, clusters)
  output$dimplot <- renderPlot(DimPlot(dat_list[[input$subset_featureplots]], group.by = input$group) +
                                 my_theme,
                               height = function() {session$clientData$output_dimplot_width * 0.8})
  
  
  subset_dataset <- reactive({dat_list[[input$select_lineage_dataset]]})
  ####################################################################
  # Generate gene feature plots
  updateSelectizeInput(session, "gene_id", choices = gene_ids, selected = "PAX7", options= list(maxOptions = length(gene_ids)), server = TRUE)
  
  output$featureplot <- renderPlot(FeaturePlot(dat_list[[input$subset_featureplots]], features = input$gene_id)+
                                     my_theme,
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
    
    updateSelectizeInput(session, "feature_lineage_dynamics", choices = features, selected = 'latent_time', server = TRUE)
    
  })
  
  output$scvelo_umaps <- renderPlot(
    FeaturePlot(subset_dataset(), features = input$feature_lineage_dynamics, pt.size = 1, order = FALSE) +
      scale_colour_viridis(option=ifelse(input$feature_lineage_dynamics == 'latent_time', 'magma', 'viridis')) +
      my_theme,
    height = function() {session$clientData$output_scvelo_umaps_width * 0.7}
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
    updateCheckboxGroupInput(session, "select_lineage", choices = lineages, selected = 'Neural crest')
  })
  
  # Set lineage dynamic feature options from server side 
  updateSelectizeInput(session, "feature_lineage_dynamics", choices = scvelo_features, server = TRUE)
  
  
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
  
  # Set lineage dynamic gene options from server side 
  updateSelectizeInput(session, "gene_id_lineage_dynamics", choices = gene_ids, selected = "PAX7", options= list(maxOptions = length(gene_ids)), server = TRUE)
  
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
                                                legend.text=element_text(size=20)),
                                        height = function() {session$clientData$output_lineage_dynamics_width * 0.5}
  )
  
  
  ####################################################################
  # Co-expression feature plots
  
  # Set coexpression gene options from server side 
  updateSelectizeInput(session, "coexpression_gene_id_1", choices = gene_ids, selected = "PAX7", options= list(maxOptions = length(gene_ids)), server = TRUE)
  updateSelectizeInput(session, "coexpression_gene_id_2", choices = gene_ids, selected = "SIX1", options= list(maxOptions = length(gene_ids)), server = TRUE)
  
  # Edit co-expression header based on genes selected
  output$coexpression_header <- renderUI({
    h4(paste('Co-expression of', input$coexpression_gene_id_1, 'and', input$coexpression_gene_id_2))
  })
  
  output$coexpression_featureplot_1 <- renderPlot(FeaturePlot(dat_list[[input$coexpression_subset_featureplots]],
                                                              features = input$coexpression_gene_id_1,
                                                              pt.size = 0.75) +
                                                    scale_colour_gradientn(colours = colorRampPalette(c("gray85", input$gene_1_col))(5)) +
                                                    my_theme,
                                                  height = function() {session$clientData$output_coexpression_featureplot_1_width * 0.75}
  )
  
  output$coexpression_featureplot_2 <- renderPlot(FeaturePlot(dat_list[[input$coexpression_subset_featureplots]],
                                                              features = input$coexpression_gene_id_2,
                                                              pt.size = 0.75) +
                                                    scale_colour_gradientn(colours = colorRampPalette(c("gray85", input$gene_2_col))(5)) +
                                                    my_theme,
                                                  height = function() {session$clientData$output_coexpression_featureplot_2_width * 0.75}
  )
  
  output$coexpression_umap <- renderPlot(coexpression_umap(dat_list[[input$coexpression_subset_featureplots]], gene_1 = input$coexpression_gene_id_1,
                                                           gene_2 = input$coexpression_gene_id_2, col.threshold = input$threshold, two.colors = c(input$gene_1_col, input$gene_2_col),
                                                           negative.color = 'gray90', highlight_cell_size = 2) +
                                           my_theme,
                                         height = function() {session$clientData$output_coexpression_umap_width * 0.7}
  )
  
  ####################################################################
  # DEA
  
  # Observe inputs into select a and select b, then remove those selected from the input options for the other
  # (this is to prevent the same cells being selected for both sides of the DEA)
  values <- reactiveValues(group_values = "",
                           a_selected = "",
                           b_selected = "")
  
  # Filter group_by choices based on dataset selected and reset selected cells
  observeEvent(input$subset_dea, {
    if(input$subset_dea %in% c('Full data', 'NPB subset')){
      updateSelectInput(session, "group_dea", choices = group_options, selected = "Cell state")
    }else{
      updateSelectInput(session, "group_dea", choices = group_options[group_options != "Stage"], selected = "Cell state")
    }
    
    # updateSelectInput(session, "select_a", selected = "")
    # updateSelectInput(session, "select_b", selected = "")
    # 
    values$group_values <- dat_list[[input$subset_dea]]@meta.data %>%
      pull(input$group_dea) %>%
      unique() %>%
      sort()
    
    updateSelectInput(session, "select_a", choices = values$group_values, selected = "")
    updateSelectInput(session, "select_b", choices = values$group_values, selected = "")
    
  })
  
  
  
  observeEvent(input$group_dea, {
    values$group_values <- dat_list[[input$subset_dea]]@meta.data %>%
      pull(input$group_dea) %>%
      unique() %>%
      sort()
    
    updateSelectInput(session, "select_a", choices = values$group_values)
    updateSelectInput(session, "select_b", choices = values$group_values)
    
  })
  
  # After selecting a, store selection in memory and remove the selection from select_b
  observeEvent(input$select_a, {
    values$a_selected = input$select_a
    
    b_choices <- isolate(values$group_values)
    b_choices <- b_choices[!b_choices %in% values$a_selected]
    updateSelectInput(session, "select_b", choices = b_choices, selected = values$b_selected)
  })
  
  # After selecting b, store selection in memory and remove the selection from select_a
  observeEvent(input$select_b, {
    values$b_selected = input$select_b
    
    a_choices <- isolate(values$group_values)
    a_choices <- a_choices[!a_choices %in% values$b_selected]
    
    updateSelectInput(session, "select_a", choices = a_choices, selected = values$a_selected)
  })
  
  
  # Render empty data table to hide spinning wheel
  output$dea_table <- DT::renderDT(NULL)
  # isolate execution to prevent updating reactive variables from updating DEA output
  observeEvent(input$run_dea, {
    req(input$select_a, input$select_b, input$padj_threshold, input$fc_threshold)
    output$dea_table <- DT::renderDataTable(
      {isolate(FindMarkers(dat_list[[input$subset_dea]], ident.1 = input$select_a, ident.2 = input$select_b,
                           group.by = input$group_dea, logfc.threshold = input$fc_threshold, base = 2) %>%
                 rownames_to_column('Gene') %>%
                 filter(p_val_adj < input$padj_threshold) %>%
                 rename("p-val" = p_val, "log2FC" = avg_log2FC, "Percent A" = pct.1, "Percent B" = pct.2, "Adj p-val" = p_val_adj) %>%
                 mutate_if(is.numeric, signif, 3) %>%
                 arrange(-log2FC))},
      rownames = FALSE, options = list(scrollX = TRUE), selection = 'none'
    )
  })
  
  # Add download button (hide button if options for DEA are not selected)
  output$download <- renderUI({
    req(input$select_a, input$select_b, input$padj_threshold, input$fc_threshold)
    downloadButton('download_csv', 'Download DEA')
  })
  
  output$download_csv <- downloadHandler("Thiery_et_al_dea.csv", content = function(file) {
    output_data <- isolate(FindMarkers(dat_list[[input$subset_dea]], ident.1 = input$select_a, ident.2 = input$select_b,
                                       group.by = input$group_dea, logfc.threshold = input$fc_threshold, base = 2) %>%
                             rownames_to_column('Gene') %>%
                             filter(p_val_adj < input$padj_threshold) %>%
                             rename("p-val" = p_val, "log2FC" = avg_log2FC, "Percent A" = pct.1, "Percent B" = pct.2, "Adj p-val" = p_val_adj) %>%
                             mutate_if(is.numeric, signif, 3) %>%
                             arrange(-log2FC))
    write.table(output_data, file  ,sep=",",row.names = F)
  })
}
