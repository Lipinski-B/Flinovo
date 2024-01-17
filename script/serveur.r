#!/usr/bin/Rscript
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(stringr)
library(SingleR)
library(celldex)
library(umap)

#score
library(ggplot2)
library(GEOquery)
library(AUCell)
library(data.table)
library(GSEABase)
library(gplots)
library(cowplot) ; theme_set(theme_cowplot())

setwd(dir = "E:/R/")
path = ""
integration <- function(singlet,Condition,Phenotype){
  # Subset ------------------------------------------------------------------
  singlet <- seurat_subset(singlet, "Condition", Condition)
  singlet <- seurat_subset(singlet, "Phenotype", Phenotype ); gc();gc();gc();

  # Integration -------------------------------------------------------------
  singlet <- SplitObject(singlet, split.by = "Sample")
  
  #for (i in c("FL08G0293","FL09C1164", "FL05G0330", "FL140304", "FL02G095", "FL12C1888", "FL08G0431", "FL06G1206", "FL08G0404")) {
    #singlet[[i]] <- NormalizeData(singlet[[i]], verbose = TRUE)
    #singlet[[i]] <- SCTransform(singlet[[i]], assay = "RNA", method = "glmGamPoi", conserve.memory = T, ncells = 10000, variable.features.n = 3000)
  #singlet[[i]] <- NormalizeData(singlet[[i]], assay = "RNA")
  #singlet[[i]] <- ScaleData(singlet[[i]], features = rownames(singlet[[i]]))
  #singlet[[i]] <- SCTransform(singlet[[i]], vst.flavor = "v2", verbose = FALSE) %>% RunPCA(npcs = 30, verbose = FALSE)
  #} 
  
  singlet <- lapply(X = singlet, FUN = function(x) {
    x <- NormalizeData(x, assay = "RNA") ; x <- ScaleData(x, features = rownames(x))
    x <- SCTransform(x, vst.flavor = "v2", ncells = 10000, conserve.memory = T, verbose = FALSE) %>% RunPCA(npcs = 30, verbose = FALSE)
  })

  integ_features <- SelectIntegrationFeatures(object.list = singlet, nfeatures = 3000) 
  singlet <- PrepSCTIntegration(object.list = singlet, anchor.features = integ_features)

  if(Condition == "RCHOP" || Condition == "Excipient"){
    if(Phenotype == c("CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")){k = 5  ; d = 5  ; s = 5 } 
    else                                                                                                                 {k = 20 ; d = 20 ; s = 30}
  } else                                                                                                                 {k = 100; d = 20 ; s = 30}
  
  integ_anchors <- FindIntegrationAnchors(object.list = singlet, normalization.method = "SCT", anchor.features = integ_features, dim = 1:d, k.score = s) ;
  singlet <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT", k.weight = k, features.to.integrate = integ_features, dims = 1:d) ; rm(integ_anchors, integ_anchors) ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; 
  
  # Visualisation -----------------------------------------------------------
  #feature <- VariableFeatures(singlet)[which(!str_detect(VariableFeatures(singlet), "IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]"))];
  #singlet <- RunPCA(singlet, verbose = F, approx=F, feature = feature) %>% RunUMAP(verbose = T, umap.method = "umap-learn", feature = feature) %>% FindNeighbors(feature = feature, verbose = F) %>% FindClusters(resolution = 0.6, verbose = F) %>% RunTSNE(verbose = F, perplexity = 100, feature = feature) #%>% RunVelocity(deltaT = 1, kCells = 50, fit.quantile = 0.02, na.rm = T, ncores = 12)
  
  # Score AUCell ------------------------------------------------------------
  Singature_RCHOP <- c("BAX","RPS27L","PHPT1","SRSF3","PSMB4","HIST1H2BK","TRIM22","AEN","PVT1","FDXR","BBC3","SRSF2","MRFAP1","DDB2","EIF2S3","MDM2","HNRNPH1","CHI3L2","CCNG1","LY86","ZMAT3","TXNIP","CD70","SNHG8","P4HA1","CDKN1A","ISG15")
  #Singature_RCHOP <- c("PPIA","RPS19","HIST1H2BK","BAX","RPS27L","SRSF3","PHPT1","PSMB4","TRIM22","VPREB3","BBC3","SNRPF","EIF2S3","CRIP1","FDXR","PSMB4","AEN","HSPA1A") 
  #Singature_RCHOP <- c("PHPT1","HIST1H2BK", "RPS27L","SRSF3","TRIM22","BBC3","FDXR","PSMB4","AEN","HSPA1A","BAX")
  cells_rankings <- AUCell_buildRankings(singlet@assays$RNA@counts, nCores=1, plotStats=FALSE)
  geneSets <- list(geneSet1=Singature_RCHOP)
  
  score <- AUCell_calcAUC(geneSets, cells_rankings) ; rm(cells_rankings) ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ;
  score <- AUCell::getAUC(score) ;  score <- as.matrix(t(score))
  singlet <- AddMetaData(object=singlet, metadata = score, col.name = "RCHOP_AUC_Score") ; rm(score) ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ;
  
  return(singlet)
}
seurat_subset <- function(singlet, item, sub_item){
  Idents(singlet) <- item
  sub_singlet <- subset(singlet, idents = sub_item)
  return(sub_singlet)
}

# Run ---------------------------------------------------------------------
load(file=paste0(path,"All_0.RData"))
Idents(singlet) <- "Condition"
singlet2 <- singlet

# #All
# singlet <- integration(singlet2, c("Pré-greffe", "Excipient","RCHOP"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ; save(singlet, file=paste0(path,"All_All_All.RData"))
# singlet <- integration(singlet2, c("Pré-greffe", "Excipient","RCHOP"), c("CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;            save(singlet, file=paste0(path,"All_All_Other.RData"))
# singlet <- integration(singlet2, c("Pré-greffe", "Excipient","RCHOP"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells")) ;                                                         save(singlet, file=paste0(path,"All_All_BT.RData"))
# singlet <- integration(singlet2, c("Pré-greffe", "Excipient","RCHOP"), c("B cells")) ;                                                                                                    save(singlet, file=paste0(path,"All_All_B.RData"))
# 
# #Pré-greffe-Excipient
# singlet <- integration(singlet2, c("Pré-greffe", "Excipient"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;         save(singlet, file=paste0(path,"All_Pré-greffe-Excipient_All.RData"))
# singlet <- integration(singlet2, c("Pré-greffe", "Excipient"), c("CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;                    save(singlet, file=paste0(path,"All_Pré-greffe-Excipient_Other.RData"))
# singlet <- integration(singlet2, c("Pré-greffe", "Excipient"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells")) ;                                                                 save(singlet, file=paste0(path,"All_Pré-greffe-Excipient_BT.RData"))
# singlet <- integration(singlet2, c("Pré-greffe", "Excipient"), c("B cells")) ;                                                                                                            save(singlet, file=paste0(path,"All_Pré-greffe-Excipient_B.RData"))

#Post-greffe
singlet <- integration(singlet2, c("RCHOP", "Excipient"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;              save(singlet, file=paste0(path,"All_Post-greffe_All.RData"))
singlet <- integration(singlet2, c("RCHOP", "Excipient"), c("CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;                      save(singlet, file=paste0(path,"All_Post-greffe_Other.RData"))
singlet <- integration(singlet2, c("RCHOP", "Excipient"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells")) ;                                                                      save(singlet, file=paste0(path,"All_Post-greffe_BT.RData"))
singlet <- integration(singlet2, c("RCHOP", "Excipient"), c("B cells")) ;                                                                                                                 save(singlet, file=paste0(path,"All_Post-greffe_B.RData"))

#Pré-greffe
singlet <- integration(singlet2, c("Pré-greffe"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;                      save(singlet, file=paste0(path,"All_Pré-greffe_All.RData"))
singlet <- integration(singlet2, c("Pré-greffe"), c("CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;                                 save(singlet, file=paste0(path,"All_Pré-greffe_Other.RData"))
singlet <- integration(singlet2, c("Pré-greffe"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells")) ;                                                                              save(singlet, file=paste0(path,"All_Pré-greffe_BT.RData"))
singlet <- integration(singlet2, c("Pré-greffe"), c("B cells")) ;                                                                                                                         save(singlet, file=paste0(path,"All_Pré-greffe_B.RData"))

#Excipient
singlet <- integration(singlet2, c("Excipient"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;                       save(singlet, file=paste0(path,"All_Excipient_All.RData"))
singlet <- integration(singlet2, c("Excipient"), c("CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;                                  save(singlet, file=paste0(path,"All_Excipient_Other.RData"))
singlet <- integration(singlet2, c("Excipient"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells")) ;                                                                               save(singlet, file=paste0(path,"All_Excipient_BT.RData"))
singlet <- integration(singlet2, c("Excipient"), c("B cells")) ;                                                                                                                          save(singlet, file=paste0(path,"All_Excipient_B.RData"))

#RCHOP
singlet <- integration(singlet2, c("RCHOP"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;                           save(singlet, file=paste0(path,"All_RCHOP_All.RData"))
singlet <- integration(singlet2, c("RCHOP"), c("CD4+ T cells", "CD8+ T cells", "T cells", "Dendritic cells","Monocytes","NK cells","Progenitors")) ;                                      save(singlet, file=paste0(path,"All_RCHOP_Other.RData"))
singlet <- integration(singlet2, c("RCHOP"), c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells")) ;                                                                                   save(singlet, file=paste0(path,"All_RCHOP_BT.RData"))
singlet <- integration(singlet2, c("RCHOP"), c("B cells")) ;                                                                                                                              save(singlet, file=paste0(path,"All_RCHOP_B.RData"))
