## -- Velocity -- ##
library("velocyto.R")
library("SeuratWrappers")

velocity <- function(patient){#condition
  ## -- Load data -- ##
  load(file = paste0("/home/boris/Bureau/Flinovo/result/analyse_patient/",patient,"/All_All_",patient,".RData")) 
  ldat <- ReadVelocity(file = paste0("/home/boris/Documents/lipinskib/Boris_Manon/flinovo/result/",patient,"/mRNA/CellrangerCount/velocyto/CellrangerCount.loom"))
  
  ## -- Convert velo to seurat -- ##
  bm <- as.Seurat(x = ldat)
  nom <- c() ; for (i in 1:length(colnames(bm))){nom <- c(nom,paste0(str_split(str_split(colnames(bm)[i], ":")[[1]][2],"x")[[1]][1],"-1"))}
  bm <- RenameCells(bm,new.names = nom)
  bm <- subset(bm, cells = colnames(singlet))
  bm <- AddMetaData(bm, singlet@meta.data$Phénotype, col.name="Phénotype")
  bm <- AddMetaData(bm, singlet@meta.data$Phénotype.fine, col.name="Phénotype.fine")
  bm <- AddMetaData(bm, singlet@meta.data$Condition, col.name="Condition")
  bm <- AddMetaData(bm, singlet@meta.data$Sample, col.name="Sample")
  bm <- AddMetaData(bm, singlet@meta.data$Phase, col.name="Phase")
  bm <- AddMetaData(bm, singlet@meta.data$Chaine, col.name="Chaine")
  bm <- AddMetaData(bm, singlet@meta.data$Clonotype, col.name="Clonotype")
  bm <- AddMetaData(bm, singlet@meta.data$CDR3S_AA, col.name="CDR3S_AA")
  bm <- AddMetaData(bm, singlet@meta.data$CDR3, col.name="CDR3")
  bm <- AddMetaData(bm, singlet@meta.data$V, col.name="V")
  bm <- AddMetaData(bm, singlet@meta.data$D, col.name="D")
  bm <- AddMetaData(bm, singlet@meta.data$J, col.name="J")
  bm <- AddMetaData(bm, singlet@meta.data$C, col.name="C")
  #bm <- seurat_subset(bm, "Condition", c(condition))
  
  #bm <- seurat_subset(bm, "Phénotype", 'B cells')
  #bm <- subset(bm, nFeature_spliced > 50)
  
  ## -- Velocity calcul -- ##
  # bm <- SCTransform(bm, assay = "spliced", verbose = F, ncells = 1000, conserve.memory = T) %>% RunPCA(verbose = F) %>% FindNeighbors(dims = 1:20, verbose = F) %>% 
  #   FindClusters(verbose = F) %>% RunUMAP(dims = 1:20, verbose = F) %>% RunTSNE(dims = 1:20, verbose = F) %>% RunVelocity(deltaT = 1, kCells = 50, fit.quantile = 0.02, na.rm = T, ncores = 12)

  return(bm)
}
