library(celldex)
library(SingleR)
library(stringr)
library(dplyr)
library(future)
library(Seurat)
library(umap)
library(tibble)

#library(escape)
#library(DESeq2)
#library(dittoSeq)
#library(monocle3)

# library(enrichR)
# listEnrichrSites() 
# setEnrichrSite("Enrichr") 
# websiteLive <- TRUE 
# dbs <- listEnrichrDbs() 
# if (is.null(dbs)) websiteLive <- FALSE
# rm(dbs)
# rm(websiteLive)

## -- Loading -- ## 
options(future.globals.maxSize = 8000 * 1024^2)
siege <- c("FL140304","FL12C1888","FL09C1164","FL08G0293", "FL02G095","FL05G0330","FL08G0431", "FL06G1206", "FL08G0404") 
