mod_3D_ui <- function(id) {
  ns <- NS(id)
  # 3D Reduction dimention 
  tabItem(tabName = "3D_RD",
          radioButtons("metadata", "Metadata:",inline=T,singlet@tools$meta_variable),
          
          navbarPage("Reduction dimention",
                     tabPanel("PCA",plotlyOutput(ns("D_PCA"), width = "100%",  height = "800px"),uiOutput(ns("Dfeature_pca"))),
                     tabPanel("UMAP",plotlyOutput(ns("D_UMAP"), width = "100%",  height = "800px"),uiOutput(ns("Dfeature_umap"))),
                     tabPanel("TSNE",plotlyOutput(ns("D_TSNE"), width = "100%",  height = "800px"),uiOutput(ns("Dfeature_tsne")))
          ),
          
          fluidRow(),
          column(1,align="right",checkboxGroupInput("DFeatures", NULL, choices = list("Features" = "Features"), selected = 0)),
          
          fluidRow(),
          conditionalPanel(
            condition = "input.DFeatures == 'Features'",
            fluidRow(column(6,selectInput('Dvariables', 'Choices', NULL, selectize=TRUE, selected = NULL))),
          ),
  )
}

mod_3D_server <- function(input, output, session) {
  ns <- session$ns
  
  D3feature <- reactive({
    gene_feature <- input$Dvariables
    singlet2 <- singlet2()
    all_feature <- rownames(as.matrix(singlet2[["RNA"]]@counts))
    plot.data <- FetchData(object = singlet2, vars = c(all_feature, "PC_1", "PC_2", "PC_3", "tSNE_1", "tSNE_2", "tSNE_3", "UMAP_1", "UMAP_2", "UMAP_3"), slot = 'data')
    plot.data$changed <- ifelse(test = plot.data[[gene_feature]] <1, yes = plot.data[[gene_feature]], no = 1)
    plot.data$label <- paste(rownames(plot.data)," - ", plot.data[[gene_feature]], sep="")
    return(plot.data)
  })  
  ## -- PCA -- ##
  D3plot <- reactive({
    plot.data <- FetchData(object = singlet2(), vars = c(r$dataset@tools$meta_variable, "PC_1", "PC_2", "PC_3", "tSNE_1", "tSNE_2", "tSNE_3", "UMAP_1", "UMAP_2", "UMAP_3"))
    plot.data$label <- paste(rownames(plot.data))
    return(plot.data)
  })
  output$D_PCA <- renderPlotly({
    plot_ly(data = D3plot(), x = ~PC_1, y = ~PC_2, z = ~PC_3,color = singlet@meta.data[[input$metadata]],type = "scatter3d", mode = "markers",marker = list(size = 3, width=2),
            text=~singlet@meta.data[[input$metadata]], hoverinfo="text") %>%layout(title=input$metadata,scene = list(xaxis = list(title = paste0("PCA 1 : ", round(Seurat::Stdev(singlet[["pca"]])[1],2), " %")),yaxis = list(title = paste0("PCA 2 : ", round(Seurat::Stdev(singlet[["pca"]])[2],2), " %")),zaxis = list(title = paste0("PCA 3 : ", round(Seurat::Stdev(singlet[["pca"]])[3],2), " %"))))
  })
  output$D_UMAP <- renderPlotly({plot_ly(data = D3plot(), x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,color = singlet@meta.data[[input$metadata]],type = "scatter3d", mode = "markers",marker = list(size = 3, width=2),text=~singlet@meta.data[[input$metadata]], hoverinfo="text") %>% layout(title=input$metadata)})
  output$D_TSNE <- renderPlotly({plot_ly(data = D3plot(), x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3,color = singlet@meta.data[[input$metadata]],type =  "scatter3d", mode = "markers",marker = list(size = 3, width=2),text=~singlet@meta.data[[input$metadata]], hoverinfo="text") %>% layout(title=input$metadata)})
  
  
  output$Dfeature_pca = renderUI({if(length(input$Dvariables)>0){plotlyOutput(ns("DFeaturePlot_PCA"), width = "100%",  height = "800px")}})
  output$Dfeature_umap = renderUI({if(length(input$Dvariables)>0){plotlyOutput(ns("DFeaturePlot_UMAP"), width = "100%",  height = "800px")}})
  output$Dfeature_tsne = renderUI({if(length(input$Dvariables)>0){plotlyOutput(ns("DFeaturePlot_TSNE"), width = "100%",  height = "800px")}})
  
  
  
  output$DFeaturePlot_PCA <- renderPlotly({
    plot_ly(data = D3feature(), x = ~PC_1, y = ~PC_2, z = ~PC_3,color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
            colors = c('darkgreen', 'red'), opacity = .5,type = "scatter3d", mode = "markers", marker = list(size = 3, width=2),text=~label, hoverinfo="text")# %>% layout(title=gene_feature)
  })
  output$DFeaturePlot_UMAP<- renderPlotly({
    plot_ly(data = D3feature(), x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
            colors = c('darkgreen', 'red'), opacity = .5,type = "scatter3d", mode = "markers", marker = list(size = 3, width=2),text=~label, hoverinfo="text")# %>% layout(title=gene_feature)
  })
  output$DFeaturePlot_TSNE<- renderPlotly({
    plot_ly(data = D3feature(), x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3,color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
            colors = c('darkgreen', 'red'), opacity = .5,type = "scatter3d", mode = "markers", marker = list(size = 3, width=2),text=~label, hoverinfo="text") #%>% layout(title=gene_feature)
  })
  
}


## To be copied in the UI
# menuSubItem("3D Visualisation",            tabName = "3D_RD", icon = icon("indent-right", lib = "glyphicon")),
# mod_3D_ui("3D_ui_1")

## To be copied in the server
# callModule(mod_3D_server, "3D_ui_1")#, r = r)