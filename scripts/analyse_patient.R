source(file = "/home/boris/Bureau/Flinovo/fonction/library.R")

# Workflow ----------------------------------------------------------------
for (patient in c("FL180250B")){ # "FL08G0293","FL02G095","FL05G0330", "FL06G1206"  "FL05G0305", "FL120212","FL12C1888","FL140304", "FL08G0431", 
  
  # Seurat Object ----------------------------------------------------------
  singlet <- Seurat_object(patient)
  singlet <- SCT_Normalization(singlet)
  singlet <- Visualisation(singlet)
  singlet <- Metadata(singlet)
  
  if(is.null(singlet@meta.data$EZH2_A682G)){singlet@meta.data$EZH2_A682G<- NA}
  if(is.null(singlet@meta.data$EZH2_A692V)){singlet@meta.data$EZH2_A692V<- NA}
  if(is.null(singlet@meta.data$EZH2_Y646C)){singlet@meta.data$EZH2_Y646C<- NA}
  if(is.null(singlet@meta.data$EZH2_Y646F)){singlet@meta.data$EZH2_Y646F<- NA}
  if(is.null(singlet@meta.data$EZH2_Y646H)){singlet@meta.data$EZH2_Y646H<- NA}
  if(is.null(singlet@meta.data$EZH2_Y646N)){singlet@meta.data$EZH2_Y646N<- NA}
  if(is.null(singlet@meta.data$EZH2_Y646S)){singlet@meta.data$EZH2_Y646S<- NA}
  
  singlet@meta.data <- singlet@meta.data[, c("Sample","Condition","Phénotype","Phénotype.fine","Clonotype","Heavy","Light","Chaine_Light","V_Light","D_Light","J_Light","C_Light","CDR3_Light","CDR3_nt_Light","Chaine_Heavy","V_Heavy","D_Heavy","J_Heavy","C_Heavy","CDR3_Heavy","CDR3_nt_Heavy","Phénotype_TRUST","V_Heavy_TRUST","D_Heavy_TRUST","J_Heavy_TRUST","Heavy_TRUST","CDR3_DNA_Heavy_TRUST", "CDR3_AA_Heavy_TRUST","V_Light_TRUST","D_Light_TRUST","J_Light_TRUST","Light_TRUST","CDR3_DNA_Light_TRUST","CDR3_AA_Light_TRUST","Phase","S.Score","G2M.Score","Greffe","State","percent.mt","percent.rb","percent.ig","seurat_clusters","nCount_RNA","nFeature_RNA","nCount_SCT","nFeature_SCT","nCount_HTO","nFeature_HTO")]#,"BCL2_L23L", "BCL2_K22K", "CD79B_Y196H", "EZH2_A682G","EZH2_A692V","EZH2_Y646C","EZH2_Y646F","EZH2_Y646H","EZH2_Y646N","EZH2_Y646S")]
  Idents(singlet) <- "Sample"
  
  metadata <- singlet[[]]
  colnames(metadata)[1] <- "SampleID"
  write.csv(metadata,paste0("/home/boris/Documents/lipinskib/Boris_Manon/flinovo/result/",patient,"/mRNA/",patient,".csv"), row.names = FALSE)
  
  # Sauvegarde All ----------------------------------------------------------
  save(singlet, file = paste0("/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_patient/", patient,"/All_All_", patient,".RData"))
  Subset_Sample(singlet, c("RCHOP","Excipient"), "Post-greffe_All", patient)
  Subset_Sample(singlet, "Pré-greffe", "Pré-greffe_All", patient)
  Subset_Sample(singlet, c("Pré-greffe","Excipient"), "Pré-greffe-Excipient_All", patient)
  
  # Sub sample BT -----------------------------------------------------------
  singlet_BT <- Subset_Object(singlet, "Phénotype", c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells"))
  Subset_Sample(singlet_BT, c("Pré-greffe","RCHOP","Excipient"), "All_BT", patient)
  Subset_Sample(singlet_BT, c("RCHOP","Excipient"), "Post-greffe_BT", patient)
  Subset_Sample(singlet_BT, "Pré-greffe", "Pré-greffe_BT", patient)
  Subset_Sample(singlet_BT, c("Pré-greffe","Excipient"), "Pré-greffe-Excipient_BT", patient)
  
  # Sub sample B ------------------------------------------------------------
  singlet_B <- Subset_Object(singlet, "Phénotype", "B cells")
  Subset_Sample(singlet_B, c("Pré-greffe","RCHOP","Excipient"), "All_B", patient)
  Subset_Sample(singlet_B, c("RCHOP","Excipient"), "Post-greffe_B", patient)
  Subset_Sample(singlet_B, "Pré-greffe", "Pré-greffe_B", patient)
  Subset_Sample(singlet_B, c("Pré-greffe","Excipient"), "Pré-greffe-Excipient_B", patient)
  
  # Diet --------------------------------------------------------------------
  singlet <- Diet(singlet)
  save(singlet, file = paste0("/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/datasets/Patient/All_All_", patient,".RData"))
  
}




# Entropy -----------------------------------------------------------------
source(file = "/home/boris/Bureau/Flinovo/fonction/entropy.R")
entropy <- function(counts, size){
  cell <- sample(colnames(counts), size, replace = TRUE)   # 20 Cellules
  counts.bootstrap <- counts[,cell]       #matrice : 20 cellules x 15000 gènes
  entropy <- apply(counts.bootstrap, 1, estimateur_bub_entropie_discrete) # 15000 entropies
  #entropy.mean <- mean(entropy) #1 moyenne d'entropie par gène 
  return(entropy)
}
entropy_calcul <- function(singlet, condition, times, size){
  Condition <- seurat_subset(singlet, "Condition", condition)
  Condition.counts <- Condition@assays$RNA@counts
  Condition.entropy <- replicate(times, entropy(Condition.counts, size))
  Condition.entropy.mean <- apply(Condition.entropy, 1, mean)
  Condition.entropy.result <- data.frame(
    entropy = Condition.entropy.mean,
    contidion = rep(condition, length(Condition.entropy.mean))
  )
  ggplot(entropy.Excipient, aes(x=entropy)) + geom_histogram(colour="black", fill="white", bins = 5)
  return(Condition.entropy.result)
}

entropy.RCHOP <- entropy_calcul(singlet, "RCHOP", 15, 100)
entropy.Excipient <- entropy_calcul(singlet, "Excipient", 15, 100)
entropy.Pregreffe <- entropy_calcul(singlet, "Pré-greffe", 15, 100)

save(entropy.RCHOP, entropy.Excipient, entropy.Pregreffe, file = "/home/boris/Bureau/entropy_FL14.RData")
load(file = "/home/boris/Bureau/entropy_FL14.RData")

# Update ------------------------------------------------------------------
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
for (patient in c("FL06G1206", "FL08G0404", "FL120212", "FL05G0305")) {#"FL140304", "FL09C1164","FL08G0293", "FL02G0^95","FL05G0330", "FL08G0431", "FL06G1206", "FL12C1888" ("FL08G0431",
  
  ldat <- ReadVelocity(file = paste0("/home/boris/Documents/lipinskib/Boris_Manon/flinovo/result/",patient,"/mRNA/CellrangerCount/velocyto/CellrangerCount.loom"))
  bam <- as.Seurat(x = ldat)
  nom <- c() ; for (i in 1:length(colnames(bam))){nom <- c(nom,paste0(str_split(str_split(colnames(bam)[i], ":")[[1]][2],"x")[[1]][1],"-1"))}
  bam <- RenameCells(bam,new.names = nom)
  
  for(condition in c("All","Post-greffe","Pré-greffe","Pré-greffe-Excipient")){#
    for(phénotype in c("All","BT","B")){#"
      load(file = paste0("/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_patient/",patient,"/",condition,"_",phénotype,"_",patient,".RData"))
      
      bm <- subset(bam, cell = colnames(singlet))
      
      bm <- SCTransform(object = bm, assay = "spliced")
      bm <- RunPCA(object = bm, verbose = FALSE)
      bm <- FindNeighbors(object = bm, dims = 1:20)
      bm <- FindClusters(object = bm)
      bm <- RunUMAP(object = bm, dims = 1:20)
      bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02)
      
      singlet[['spliced']] <- bm[['spliced']]
      singlet[['unspliced']] <- bm[['unspliced']]
      singlet[['ambiguous']] <- bm[['ambiguous']]
      #singlet <- RunVelocity(singlet, deltaT = 1, kCells = 50, fit.quantile = 0.02, na.rm = T)
      
      save(singlet, file = paste0("/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_patient/",patient,"/",condition,"_",phénotype,"_",patient,".RData"))
      singlet <- Diet(singlet)
      save(singlet, file = paste0("/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/datasets/Patient/",condition,"_",phénotype,"_",patient,".RData")) 
    }
  }
}

