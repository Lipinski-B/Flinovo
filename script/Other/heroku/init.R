# init.R
#
# Example R code to install packages if not already installed
#

my_packages = c('plotly','shiny','shinybusy','shinydashboard','dashboardthemes','shinyWidgets','BiocManager','SeuratObject','Seurat', 'shinyjs')


install_if_missing = function(p) {
  if (p %in% rownames(installed.packages()) == FALSE) {
    install.packages("document/package/fastmap_1.1.0.tar.gz", repos = NULL, type="source")
    install.packages("document/package/htmltools_0.5.2.tar.gz", repos = NULL, type="source")
    install.packages(p)
    BiocManager::install('dittoSeq', update = F, force = TRUE)
  }
}

invisible(sapply(my_packages, install_if_missing))