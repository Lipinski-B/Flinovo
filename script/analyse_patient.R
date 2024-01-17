path <- "/home/boris/Bureau/Flinovo/"
path.meta <- "/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_meta/"
path.patient <- "/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_patient/"
source(file = paste0(path,"fonction/library.R"))
source(file = paste0(path,"fonction/workflow.R"))
source(file = paste0(path,"fonction/metadata.R"))

setwd(dir = "/home/boris/Documents/lipinskib/Boris_Manon/flinovo/result/")
setwd(dir = "/home/boris/")


# Transcriptomic------------------------------------------------------------
for (patient in c("FL05G0330") ){ #FL08G0431 : Comment GoT & VDJ metadatas  "FL08G0293", "FL140304","FL12C1888","FL09C1164","FL02G095","FL05G0330", "FL08G0404", "FL06G1206" 
  # Seurat Object ----------------------------------------------------------
  singlet <- seurat_object(patient)
  singlet <- normalisation(singlet)
  singlet <- visualisation(singlet)
  
  # Métadata ----------------------------------------------------------------  
  backup(singlet)
  singlet <- metadata(singlet)
  
  if(is.null(singlet@meta.data$EZH2_A682G)){singlet@meta.data$EZH2_A682G<- NA}
  if(is.null(singlet@meta.data$EZH2_A692V)){singlet@meta.data$EZH2_A692V<- NA}
  if(is.null(singlet@meta.data$EZH2_Y646C)){singlet@meta.data$EZH2_Y646C<- NA}
  if(is.null(singlet@meta.data$EZH2_Y646F)){singlet@meta.data$EZH2_Y646F<- NA}
  if(is.null(singlet@meta.data$EZH2_Y646H)){singlet@meta.data$EZH2_Y646H<- NA}
  if(is.null(singlet@meta.data$EZH2_Y646N)){singlet@meta.data$EZH2_Y646N<- NA}
  if(is.null(singlet@meta.data$EZH2_Y646S)){singlet@meta.data$EZH2_Y646S<- NA}
  
  singlet@meta.data <- singlet@meta.data[, c("Sample","Condition","Phénotype","Phénotype.fine","Clonotype","Heavy","Light","Chaine_Light","V_Light","D_Light","J_Light","C_Light","CDR3_Light","CDR3_nt_Light","Chaine_Heavy","V_Heavy","D_Heavy","J_Heavy","C_Heavy","CDR3_Heavy",                                             "CDR3_nt_Heavy",
                                             "Phénotype_TRUST","V_Heavy_TRUST","D_Heavy_TRUST","J_Heavy_TRUST","Heavy_TRUST","CDR3_DNA_Heavy_TRUST", "CDR3_AA_Heavy_TRUST","V_Light_TRUST","D_Light_TRUST","J_Light_TRUST","Light_TRUST","CDR3_DNA_Light_TRUST",                                              "CDR3_AA_Light_TRUST",
                                             "Phase","S.Score","G2M.Score","Greffe","State","percent.mt","percent.rb","percent.ig","seurat_clusters","nCount_RNA","nFeature_RNA","nCount_SCT","nFeature_SCT","nCount_HTO","nFeature_HTO",
                                             "BCL2_L23L", "BCL2_K22K", "CD79B_Y196H", "EZH2_A682G","EZH2_A692V","EZH2_Y646C","EZH2_Y646F","EZH2_Y646H","EZH2_Y646N","EZH2_Y646S")]
  Idents(singlet) <- "Sample"
  
  # Sauvegarde All ----------------------------------------------------------
  save(singlet, file = paste0("/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_patient/", patient,"/All_All_", patient,".RData"))
  sub(singlet, c("RCHOP","Excipient"), "Post-greffe_All", patient)
  sub(singlet, "Pré-greffe", "Pré-greffe_All", patient)
  sub(singlet, c("Pré-greffe","Excipient"), "Pré-greffe-Excipient_All", patient)
  
  # Sub sample BT -----------------------------------------------------------
  singlet_BT <- seurat_subset(singlet, "Phénotype", c("B cells", "CD4+ T cells", "CD8+ T cells", "T cells"))
  sub(singlet_BT, c("Pré-greffe","RCHOP","Excipient"), "All_BT", patient)
  sub(singlet_BT, c("RCHOP","Excipient"), "Post-greffe_BT", patient)
  sub(singlet_BT, "Pré-greffe", "Pré-greffe_BT", patient)
  sub(singlet_BT, c("Pré-greffe","Excipient"), "Pré-greffe-Excipient_BT", patient)
  
  # Sub sample B ------------------------------------------------------------
  singlet_B <- seurat_subset(singlet, "Phénotype", "B cells")
  sub(singlet_B, c("Pré-greffe","RCHOP","Excipient"), "All_B", patient)
  sub(singlet_B, c("RCHOP","Excipient"), "Post-greffe_B", patient)
  sub(singlet_B, "Pré-greffe", "Pré-greffe_B", patient)
  sub(singlet_B, c("Pré-greffe","Excipient"), "Pré-greffe-Excipient_B", patient)
  
  # Diet --------------------------------------------------------------------
  singlet <- diet(singlet)
  save(singlet, file = paste0("/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/datasets/Patient/All_All_", patient,".RData"))
}


# Signature Greffe --------------------------------------------------------
result <- list() ; Condition <- "Condition" ; item1 <- "Pré-greffe" ; item2 <- "Excipient" ; FC <- 1.189207
annotations <- read.csv(paste0(path,"document/annotation.csv"))

for (patient in siege) {
  load(file = paste0("/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_patient/", patient,"/Pré-greffe-Excipient_All_", patient,".RData"))
  Idents(singlet) <- Condition
  
  singlet <- NormalizeData(singlet, assay = "RNA")
  singlet <- ScaleData(singlet, features = rownames(singlet))
  response <- FindMarkers(singlet, ident.1 = item1, ident.2 = item2 , logfc.threshold = log2(FC), assay = "RNA", slot = "data")
  result[[patient]] <- response %>% rownames_to_column(var="gene") %>% left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% as_tibble(response)

  # singlet <- SCTransform(singlet, verbose = F, ncells = 10000, conserve.memory = T, vst.flavor = "v2")
  # singlet <- PrepSCTFindMarkers(singlet)
  # response <- FindMarkers(singlet,ident.1 = item1, ident.2 = item2, logfc.threshold = log2(FC), assay = "SCT", verbose = F)
  # result[[patient]] <- response %>% rownames_to_column(var="gene") %>% left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% as_tibble(response) 
  
}

intersect.gene.greffe <- Reduce(intersect, list(result[["FL140304"]]$gene,result[["FL12C1888"]]$gene,result[["FL09C1164"]]$gene,result[["FL08G0293"]]$gene,result[["FL02G095"]]$gene,result[["FL05G0330"]]$gene,result[["FL08G0431"]]$gene,result[["FL06G1206"]]$gene,result[["FL08G0404"]]$gene)) #492
save(intersect.gene.greffe, file="/home/boris/Bureau/Flinovo/script/Signature/Greffe/DB/intersect.gene.greffe")
load(file="/home/boris/Bureau/Flinovo/script/Signature/Greffe/DB/intersect.gene.greffe")


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
for (patient in siege) {#"FL140304", "FL09C1164","FL08G0293", "FL02G095","FL05G0330", "FL08G0431", "FL06G1206", "FL12C1888"
  for(condition in c("RCHOP", "Excipient","Pré-greffe")){#"All","Post-greffe","Pré-greffe","Pré-greffe-Excipient"
    for(phénotype in c("All","BT","B")){
      load(file = paste0("/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_patient/",patient,"/",condition,"_",phénotype,"_",patient,".RData"))
      
      save(singlet, file = paste0("/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_patient/",patient,"/",condition,"_",phénotype,"_",patient,".RData"))
      singlet <- diet(singlet)
      save(singlet, file = paste0("/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_patient/datasets/Patient/",condition,"_",phénotype,"_",patient,".RData")) 
    }
  }
}
