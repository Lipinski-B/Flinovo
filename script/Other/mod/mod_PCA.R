mod_PCA_ui <- function(id) {
  ns <- NS(id)
  
  # 10 PCA dimensions
  tabItem(tabName = "PCA_facteurs",
          navbarPage("PCA paramètres : ",
                     tabPanel("Facteurs", verbatimTextOutput(ns("PCA10"))),
                     tabPanel("ElbowPlot", plotOutput(ns("ElbowPlot"), width = "100%",  height = "1000px")),
                     tabPanel("DimHeatmap", plotOutput(ns("DimHeatmap"), width = "100%",  height = "1000px")),
                     tabPanel("VizDimLoadings", plotOutput(ns("VizDimLoadings"), width = "100%",  height = "1000px")),
                     tabPanel("JackStrawPlot", plotOutput(ns("JackStrawPlot"), width = "100%",  height = "1000px"))
          ),
  )
  
  
}

mod_PCA_server <- function(input, output, session, r) {
  ns <- session$ns
  
  # -- PCA Figure -- ##
  output$PCA10 <- renderPrint({print(r$dataset[["pca"]], dims = 1:10, nfeatures = 10)})
  output$ElbowPlot <- renderPlot({ElbowPlot(r$dataset, ndims = 50, reduction = "pca")})
  output$DimHeatmap <- renderPlot({DimHeatmap(r$dataset, dims = 1:10, cells = 100, balanced = TRUE)})
  output$VizDimLoadings <- renderPlot({VizDimLoadings(r$dataset, dims = 1:5, reduction = "pca")})
  output$JackStrawPlot <- renderPlot({JackStrawPlot(r$dataset, dims = 1:15)})
  
}


## To be copied in the UI
# menuSubItem("Analyses PCA",                tabName = "PCA_facteurs", icon = icon("lock", lib = "glyphicon")),
# mod_PCA_ui("PCA_ui_1"),

## To be copied in the server
# callModule(mod_PCA_server, "PCA_ui_1")#, r = r)