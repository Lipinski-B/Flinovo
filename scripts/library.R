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
# options(future.globals.maxSize = 8000 * 1024^2)
# siege <- c("FL140304","FL12C1888","FL09C1164","FL08G0293", "FL02G095","FL05G0330","FL08G0431", "FL06G1206", "FL08G0404", "FL120212", "FL05G0305", "FL06G1535", "FL180250B", "FL190383B") 
# path <- "/home/boris/Bureau/Flinovo/"
# path.meta <- "/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_meta/"
# path.patient <- "/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_patient"
# source(file = paste0(path,"fonction/workflow.R"))
# setwd(dir = "/home/boris/Documents/lipinskib/Boris_Manon/flinovo/result/")
# condition <- "RCHOP"
# phénotype <- "Other"