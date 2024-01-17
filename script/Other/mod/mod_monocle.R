mod_Monocle_ui <- function(id) {
  ns <- NS(id)
  
  # Monocle
  tabItem(tabName = "monocle",
          navbarPage("Analyses: ",
                     tabPanel("Trajectory", 
                              radioGroupButtons(inputId = "color_trajectory", label = "Color by :", choices = singlet@tools$meta_variable, justified = TRUE, checkIcon = list(yes = icon("ok",lib = "glyphicon"))),
                              fluidRow(),
                              plotOutput(ns("cell_trajectory"), width = "100%",  height = "600px")
                     ),
                     tabPanel("Ordering", 
                              radioGroupButtons(inputId = "color_order", label = "Color by :", choices = singlet@tools$meta_variable, justified = TRUE, checkIcon = list(yes = icon("ok",lib = "glyphicon"))),
                              fluidRow(),
                              selectInput('DDvariables', 'Choices', NULL, selectize=TRUE, multiple=TRUE, selected = NULL),
                              uiOutput(ns("DDorder_trajectory"))
                     )
                     
          ),
  )
}

mod_Monocle_server <- function(input, output, session) {
  ns <- session$ns
  
  ## -- Monocle -- ##
  output$cell_trajectory <- renderPlot({plot_cell_trajectory(data, color_by = input$color_trajectory)})
  output$order_trajectory <- renderPlot({
    data_expressed_genes <-  row.names(subset(fData(data), num_cells_expressed >= 10))
    data_filtered <- data[data_expressed_genes,]
    my_genes <- row.names(subset(fData(data_filtered), gene_short_name %in% input$in8))
    cds_subset <- data_filtered[my_genes,]
    p <- plot_genes_in_pseudotime(cds_subset, color_by = input$color_order)
    p + theme(text = element_text(size=20), legend.title = element_text(color = "black", size = 14),legend.text = element_text(color = "black", size = 14),axis.text=element_text(size=14),axis.title=element_text(size=14))
  })
  
  
  updateSelectizeInput(session, 'DDvariables', choices = feature, server = TRUE)
  output$DDorder_trajectory = renderUI({if(length(input$DDvariables)>0){plotOutput("order_trajectory", width = "100%",  height = "1000px")}})
  
  
}


## To be copied in the UI
# menuSubItem("Trajectoires évolutives",     tabName = "monocle", icon = icon("sort-by-attributes-alt", lib = "glyphicon")),
# mod_Monocle_ui("Monocle_ui_1"),

## To be copied in the server
# callModule(mod_Monocle_server, "Monocle_ui_1", r = r)