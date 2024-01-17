#!/usr/bin/Rscript
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(stringr)
library(SingleR)
library(celldex)
#library(umap)

#score
library(ggplot2)
library(GEOquery)
library(AUCell)
library(data.table)
library(GSEABase)
library(gplots)
library(cowplot) ; theme_set(theme_cowplot())

setwd(dir = "/sps/sallegen/boris/R/")

setwd(dir = "/home/boris/Documents/lipinskib/boris/R/")


# Function ----------------------------------------------------------------
Subset      <- function(singlet, item, sub_item)      {
  Idents(singlet) <- item
  singlet <- subset(singlet, idents = sub_item)
  return(singlet)
}
Integration <- function(singlet, Condition, Phenotype){
  
  singlet <- SplitObject(singlet, split.by = "Sample")
  singlet <- lapply(X = singlet, FUN = function(x) {
    x <- NormalizeData(x, assay = "RNA") ; x <- ScaleData(x, features = rownames(x))
    x <- SCTransform(x, vst.flavor = "v2", conserve.memory = T, verbose = FALSE) %>% RunPCA(npcs = 5, verbose = FALSE)
  })
  
  integ_features <- SelectIntegrationFeatures(object.list = singlet, nfeatures = 3000) 
  singlet        <- PrepSCTIntegration(object.list = singlet, anchor.features = integ_features)
  
  `%!in%` <- Negate(`%in%`)
  if(Condition == "RCHOP" || Condition == "Excipient" || Condition == "Pre-greffe"){
    if("B cells" %!in% Phenotype) {k = 5  ; d = 3  ; s = 5} 
    else                          {k = 20 ; d = 20 ; s = 30}
  } else                          {k = 100; d = 20 ; s = 30}
  
  k = 20 ; d = 20 ; s = 30
  
  integ_anchors <- FindIntegrationAnchors(object.list = singlet, normalization.method = "SCT", anchor.features = integ_features, dim = 1:d, k.score = s) ;
  singlet       <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT", k.weight = k, features.to.integrate = integ_features, dims = 1:d) ; 
  
  rm(integ_anchors, integ_anchors) ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; 
  
  return(singlet)
}
AUCell      <- function(singlet)                      {
  # Score AUCell ------------------------------------------------------------
  Singature_RCHOP <- c("BAX","RPS27L","PHPT1","SRSF3","PSMB4","HIST1H2BK","TRIM22","AEN","PVT1","FDXR","BBC3","SRSF2","MRFAP1","DDB2","EIF2S3","MDM2","HNRNPH1","CHI3L2","CCNG1","LY86","ZMAT3","TXNIP","CD70","SNHG8","P4HA1","CDKN1A","ISG15")
  #Singature_RCHOP <- c("PPIA","RPS19","HIST1H2BK","BAX","RPS27L","SRSF3","PHPT1","PSMB4","TRIM22","VPREB3","BBC3","SNRPF","EIF2S3","CRIP1","FDXR","PSMB4","AEN","HSPA1A") 
  #Singature_RCHOP <- c("PHPT1","HIST1H2BK", "RPS27L","SRSF3","TRIM22","BBC3","FDXR","PSMB4","AEN","HSPA1A","BAX")
  cells_rankings <- AUCell_buildRankings(singlet@assays$RNA@counts, nCores=1, plotStats=FALSE, splitByBlocks=TRUE)
  geneSets <- list(geneSet1=Singature_RCHOP)
  
  score <- AUCell_calcAUC(geneSets, cells_rankings) ; rm(cells_rankings) ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ;
  score <- AUCell::getAUC(score) ;  score <- as.matrix(t(score))
  singlet <- AddMetaData(object=singlet, metadata = score, col.name = "RCHOP_AUC_Score") ; rm(score) ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ;
  return(singlet)
}


# Workflow ----------------------------------------------------------------
Workflow    <- function(singlet, Condition, Phenotype){
  
  singlet <- Subset(singlet, "Condition", Condition)
  singlet <- Subset(singlet, "Phenotype", Phenotype ); gc();gc();gc();
  singlet <- Integration(singlet, Condition, Phenotype)
  singlet <- AUCell(singlet)
  
  return(singlet)
}

# Run ---------------------------------------------------------------------
load(file="All_0.RData")
Idents(singlet) <- "Condition"
singlet2 <- singlet

#All
#singlet <- Workflow(singlet2, c("Pre-greffe", "Excipient","RCHOP"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ; save(singlet, file="All_All_All.RData")
#singlet <- Workflow(singlet2, c("Pre-greffe", "Excipient","RCHOP"), c("CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;            save(singlet, file="All_All_Other.RData")
#singlet <- Workflow(singlet2, c("Pre-greffe", "Excipient","RCHOP"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells")) ;                                                         save(singlet, file="All_All_BT.RData")
#singlet <- Workflow(singlet2, c("Pre-greffe", "Excipient","RCHOP"), c("B cells")) ;                                                                                                    save(singlet, file="All_All_B.RData")

#Pre-greffe-Excipient
#singlet <- Workflow(singlet2, c("Pre-greffe", "Excipient"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;         save(singlet, file="All_Pre-greffe-Excipient_All.RData")
#singlet <- Workflow(singlet2, c("Pre-greffe", "Excipient"), c("CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;                    save(singlet, file="All_Pre-greffe-Excipient_Other.RData")
#singlet <- Workflow(singlet2, c("Pre-greffe", "Excipient"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells")) ;                                                                 save(singlet, file="All_Pre-greffe-Excipient_BT.RData")
#singlet <- Workflow(singlet2, c("Pre-greffe", "Excipient"), c("B cells")) ;                                                                                                            save(singlet, file="All_Pre-greffe-Excipient_B.RData")

#Post-greffe
#singlet <- Workflow(singlet2, c("RCHOP", "Excipient"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;              save(singlet, file="All_Post-greffe_All.RData")
#singlet <- Workflow(singlet2, c("RCHOP", "Excipient"), c("CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;                         save(singlet, file="All_Post-greffe_Other.RData")
#singlet <- Workflow(singlet2, c("RCHOP", "Excipient"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells")) ;                                                                      save(singlet, file="All_Post-greffe_BT.RData")
#singlet <- Workflow(singlet2, c("RCHOP", "Excipient"), c("B cells")) ;                                                                                                                 save(singlet, file="All_Post-greffe_B.RData")

#Pre-greffe
#singlet <- Workflow(singlet2, c("Pre-greffe"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;                      save(singlet, file="All_Pre-greffe_All.RData")
#singlet <- Workflow(singlet2, c("Pre-greffe"), c("CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;                                 save(singlet, file="All_Pre-greffe_Other.RData")
#singlet <- Workflow(singlet2, c("Pre-greffe"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells")) ;                                                                              save(singlet, file="All_Pre-greffe_BT.RData")
#singlet <- Workflow(singlet2, c("Pre-greffe"), c("B cells")) ;                                                                                                                         save(singlet, file="All_Pre-greffe_B.RData")

#Excipient
#singlet <- Workflow(singlet2, c("Excipient"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;                       save(singlet, file="All_Excipient_All.RData")
#singlet <- Workflow(singlet2, c("Excipient"), c("CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;                                  save(singlet, file="All_Excipient_Other.RData")
#singlet <- Workflow(singlet2, c("Excipient"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells")) ;                                                                               save(singlet, file="All_Excipient_BT.RData")
#singlet <- Workflow(singlet2, c("Excipient"), c("B cells")) ;                                                                                                                          save(singlet, file="All_Excipient_B.RData")

#RCHOP
#singlet <- Workflow(singlet2, c("RCHOP"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;                           save(singlet, file="All_RCHOP_All.RData")
#singlet <- Workflow(singlet2, c("RCHOP"), c("CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;                                      save(singlet, file="All_RCHOP_Other.RData")
#singlet <- Workflow(singlet2, c("RCHOP"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells")) ;                                                                                   save(singlet, file="All_RCHOP_BT.RData")
singlet <- Workflow(singlet2, c("RCHOP"), c("B cells")) ;                                                                                                                              save(singlet, file="All_RCHOP_B.RData")
