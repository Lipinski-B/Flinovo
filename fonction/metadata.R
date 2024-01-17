## -- Metadata -- ##
metadata <- function(singlet){
  ## -- Phénotype -- ##
  load(file = paste0("/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_patient/",patient,"/phenotype_", patient,".RDs"))
  singlet <- AddMetaData(singlet, labels, col.name = "Phénotype")
  singlet <- AddMetaData(singlet, labels.fine , col.name = "Phénotype.fine")

  ## -- GoT -- ##
  load(file=paste0("/home/boris/Documents/lipinskib/Boris_Manon/flinovo/result/", patient, "/R/", patient, "_GoT.Rdata"))
  singlet <- AddMetaData(object = singlet, metadata = GOT)
  
  ## -- Other fonction -- ##
  singlet <- cycle(singlet)         ## -- Cycle -- ##
  singlet <- VDJ(singlet, patient)  ## -- VDJ -- ##
  singlet <- other(singlet)         ## -- Other -- ##
  singlet <- sunburst(singlet)      ## -- Sunburst -- ##
  
  ## -- DE -- ##
  load(file = paste0("/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_patient/",patient,"/DE_", patient,".RDs"))
  singlet@tools$RNA_DE_RE <- RNA_DE_RE
  singlet@tools$SCT_DE_RE <- SCT_DE_RE
  singlet@tools$RNA_DE_PE <- RNA_DE_PE
  singlet@tools$SCT_DE_PE <- SCT_DE_PE
  
  Idents(singlet) <- "seurat_clusters"
  
  return(singlet)
}

phenotype <- function(singlet){
  BD <- celldex::MonacoImmuneData()
  
  results.blueprint <- SingleR::SingleR(test = singlet@assays$RNA@data, ref = BD, labels = BD$label.main)
  results.blueprint.fine <- SingleR::SingleR(test = singlet@assays$RNA@data, ref = BD, labels = BD$label.fine)

  labels <- results.blueprint$labels
  labels.fine <- results.blueprint.fine$labels

  save(labels, labels.fine, file = paste0("/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_patient/", patient,"/phenotype_", patient,".RDs"))
}
DE <- function(singlet){
  Idents(singlet)<-"Condition"
  annotations <- read.csv(paste0(path,"document/annotation.csv"))
  singlet@assays$RNA@data <- singlet@assays$RNA@data[rownames(singlet@assays$RNA@data)[which(!str_detect(rownames(singlet@assays$RNA@data), "IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]"))],];
  singlet@assays$RNA@scale.data <- singlet@assays$RNA@scale.data[rownames(singlet@assays$RNA@scale.data)[which(!str_detect(rownames(singlet@assays$RNA@scale.data), "IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]"))],];
  singlet@assays$RNA@counts <- singlet@assays$RNA@counts[rownames(singlet@assays$RNA@counts)[which(!str_detect(rownames(singlet@assays$RNA@counts), "IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]"))],];

  ## RCHOP
  RNA_DE_RE <- FindMarkers(singlet, ident.1 = "RCHOP", ident.2 = "Excipient", test.use = "negbinom", logfc.threshold = 0.15, assay = "RNA", slot = "data")
  RNA_DE_RE <- RNA_DE_RE %>% rownames_to_column(var="gene") %>% left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% as_tibble(RNA_DE_RE)

  singlet <- PrepSCTFindMarkers(singlet)
  SCT_DE_RE <- FindMarkers(singlet, ident.1 = "RCHOP", ident.2 = "Excipient", logfc.threshold = 0.15, assay = "SCT", verbose = F)
  SCT_DE_RE <- SCT_DE_RE %>% rownames_to_column(var="gene") %>% left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% as_tibble(SCT_DE_RE) 
  
  ## Greffe
  RNA_DE_PE <- FindMarkers(singlet, ident.1 = "Pré-greffe", ident.2 = "Excipient", test.use = "negbinom", logfc.threshold = 0.25, assay = "RNA", slot = "data")
  RNA_DE_PE <- RNA_DE_PE %>% rownames_to_column(var="gene") %>% left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% as_tibble(RNA_DE_PE)
  
  SCT_DE_PE <- FindMarkers(singlet, ident.1 = "Pré-greffe", ident.2 = "Excipient", logfc.threshold = 0.15, assay = "SCT", verbose = F)
  SCT_DE_PE <- SCT_DE_PE %>% rownames_to_column(var="gene") %>% left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% as_tibble(SCT_DE_PE) 
  
  save(RNA_DE_RE, SCT_DE_RE, RNA_DE_PE, SCT_DE_PE, file = paste0("/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_patient/", patient,"/DE_", patient,".RDs"))
  
}
backup <- function(singlet){
  phenotype(singlet)
  DE(singlet)
}

cycle  <- function(singlet){
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  singlet <- CellCycleScoring(singlet, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)    # Attribution à chaque cellule d'un état de phase
  #singlet <- ScaleData(singlet, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(singlet)) # Correction de la matrice d'expression dépendemment du cycle cellulaire
  Idents(singlet)<-"seurat_clusters"
  return(singlet)
}
VDJ <- function(singlet, patient){
  ## -- TRUST4 -- ##
  col <- c("Phénotype","V_Heavy","D_Heavy","J_Heavy","Heavy","CDR3_DNA_Heavy","CDR3_AA_Heavy","Read_count_Heavy","Consensus_ID_Heavy","CDR3_germline_similarity_Heavy", "Consensus_full_length_Heavy","V_Light","D_Light","J_Light","Light","CDR3_DNA_Light","CDR3_AA_Light","Read_count_Light","Consensus_ID_Light","CDR3_germline_similarity_Light", "Consensus_full_length_Light")
  VDJ <- as.data.frame(setNames(replicate(length(col),numeric(0), simplify = F),col ))
  path <- paste0('/home/boris/Documents/lipinskib/Boris_Manon/flinovo/script/TRUST4/result/',patient,'/') ; files <- list.files(path = path, pattern = "\\_barcode_report.tsv$")
  result <- read.table(file = paste0(path,files[1]), sep = '\t', header = F, row.names = 1) ; colnames(result) <- c("Phénotype","Heavy","Light","secondary_chain1","secondary_chain2") ; result[,"Phénotype"][which(result[,"Phénotype"] == "abT")] <- "T"
  Heavy <- str_split(result$Heavy, ",", simplify=T) ; colnames(Heavy) <- c("V_Heavy","D_Heavy","J_Heavy","Heavy","CDR3_DNA_Heavy","CDR3_AA_Heavy", "Read_count_Heavy", "Consensus_ID_Heavy", "CDR3_germline_similarity_Heavy", "Consensus_full_length_Heavy")
  Light <- str_split(result$Light, ",", simplify=T) ; colnames(Light) <- c("V_Light","D_Light","J_Light","Light","CDR3_DNA_Light","CDR3_AA_Light", "Read_count_Light", "Consensus_ID_Light", "CDR3_germline_similarity_Light", "Consensus_full_length_Light")
  finale <- cbind(result$barcode, result$Phénotype, Heavy, Light) ; rownames(finale) <- paste0(rownames(result),'-1') ; colnames(finale)[1] <- "Phénotype" #; intercell <- intersect(colnames(singlet),  rownames(finale)) ; cat("Nombre de cellule : ", length(intercell), " / ", length(colnames(singlet)))
  VDJ <- rbind(VDJ, finale)
  
  VDJ <- VDJ[,c("Phénotype","V_Heavy","D_Heavy","J_Heavy","Heavy","CDR3_DNA_Heavy","CDR3_AA_Heavy","V_Light","D_Light","J_Light","Light","CDR3_DNA_Light","CDR3_AA_Light")]
  colnames(VDJ) <- c("Phénotype_TRUST","V_Heavy_TRUST","D_Heavy_TRUST","J_Heavy_TRUST","Heavy_TRUST","CDR3_DNA_Heavy_TRUST","CDR3_AA_Heavy_TRUST","V_Light_TRUST","D_Light_TRUST","J_Light_TRUST","Light_TRUST","CDR3_DNA_Light_TRUST","CDR3_AA_Light_TRUST")
  
  singlet <- AddMetaData(object=singlet, metadata = VDJ)
  
  ## -- Cellranger VDJ -- ##
  bcr <- read.csv(paste0(patient,"/VDJ/CellrangerVDJ/outs/filtered_contig_annotations.csv"))  # bcr$barcode <- gsub("-1", "", bcr$barcode)                                   # Remove the -1 at the end of each barcode: VERIFIER SI IL Y A LES TIRETS DANS OBJECT SINGLET!!!
  bcr <- bcr[,c("barcode", "raw_clonotype_id","chain","v_gene","d_gene","j_gene","c_gene","cdr3", "cdr3_nt")] ; names(bcr)[names(bcr) == "raw_clonotype_id"] <- "clonotype_id"
  light <- rbind(bcr[which(bcr$chain == "IGL"),],bcr[which(bcr$chain == "IGK"),]) ; colnames(light) <- c("barcode", "clonotype_id_light","Chaine_Light","V_Light","D_Light","J_Light","C_Light","CDR3_Light","CDR3_nt_Light")
  heavy <- bcr[which(bcr$chain == "IGH"),] ; colnames(heavy) <- c("barcode", "clonotype_id_heavy","Chaine_Heavy","V_Heavy","D_Heavy","J_Heavy","C_Heavy","CDR3_Heavy","CDR3_nt_Heavy")
  bcr <- merge(light,heavy,by="barcode")
  H <- heavy[which(heavy$barcode %in% setdiff(heavy[,"barcode"], light[,"barcode"])),] ; colnames(H) <- c("barcode", "clonotype_id_heavy","Chaine_Heavy","V_Heavy","D_Heavy","J_Heavy","C_Heavy","CDR3_Heavy","CDR3_nt_Heavy")
  L <- light[which(light$barcode %in% setdiff(light[,"barcode"], heavy[,"barcode"])),] ; colnames(L) <- c("barcode", "clonotype_id_light","Chaine_Light","V_Light","D_Light","J_Light","C_Light","CDR3_Light","CDR3_nt_Light")
  bcr <- rbind(bcr,dplyr::bind_rows(H,L))
  bcr <- bcr[!duplicated(bcr$barcode), ] ;rownames(bcr) <- bcr$barcode
  bcr$barcode <- bcr$clonotype_id_heavy <- NULL
  colnames(bcr)[1] <- "Clonotype"
  bcr$Clonotype <- as.numeric(str_replace(unlist(bcr$Clonotype,"clonotype",use.names=F), "clonotype", ""))
  singlet <- AddMetaData(object=singlet, metadata = bcr)                        #, col.name = colnames(as.data.frame(bcr)))   # Add to the Seurat object's metadata.
  
  
  
  
  ## -- State -- ##
  #levels(singlet@meta.data$Clonotype) <- c(levels(singlet@meta.data$Clonotype),9999999999)
  singlet@meta.data$Clonotype[which(is.na(singlet@meta.data$Clonotype))] <- 9999
  singlet@meta.data$State <- as.numeric(singlet@meta.data$nCount_RNA)
  singlet@meta.data$State[which(as.numeric(singlet@meta.data$Clonotype) < 5)] <- "Tumoral"
  singlet@meta.data$State[which(as.numeric(singlet@meta.data$Clonotype) >= 5)] <- "Non Tumoral"
  
  
  # Update VDJ --------------------------------------------------------------
  # load(file = "/home/boris/Bureau/Flinovo/result/analyse_meta/Trust_VDJ.RData")
  # singlet <- AddMetaData(object=singlet, metadata = Trust_VDJ)                        #, col.name = colnames(as.data.frame(bcr)))   # Add to the Seurat object's metadata.
  # 
  # load(file = "/home/boris/Bureau/Flinovo/result/analyse_meta/CellrangerVDJ.RData")
  # singlet <- AddMetaData(object=singlet, metadata = CellrangerVDJ) 
  # 
  singlet@meta.data$Light <- singlet@meta.data$Heavy <- NA
  
  #CellrangerVDJ
  singlet@meta.data$Chaine_Light <- str_replace_na(singlet@meta.data$Chaine_Light,"Unknown")
  singlet@meta.data[singlet@meta.data$Chaine_Light == "IGK","Light"] <- "IGK"
  singlet@meta.data[singlet@meta.data$Chaine_Light == "IGL","Light"] <- "IGL"
  
  singlet@meta.data$Chaine_Heavy <- str_replace_na(singlet@meta.data$Chaine_Heavy,"Unknown" )
  singlet@meta.data[singlet@meta.data$Chaine_Heavy == "IGH","Heavy"] <- "IGH"
  
  #TRUST4
  singlet@meta.data$Light_TRUST <- str_replace_na(singlet@meta.data$Light_TRUST,"Unknown" )
  singlet@meta.data$Heavy_TRUST <- str_replace_na(singlet@meta.data$Heavy_TRUST,"Unknown" )
  
  singlet@meta.data[which(!(singlet@meta.data$Light_TRUST %in% c("TRAC", "Unknown","","*"))),"Light"] <- substr(singlet@meta.data[which(!(singlet@meta.data$Light_TRUST %in% c("TRAC", "Unknown","","*"))), "Light_TRUST"],1,3)
  singlet@meta.data[which(!(singlet@meta.data$Heavy_TRUST %in% c("TRBC1", "TRBC2", "Unknown","","*"))), "Heavy"] <- substr(singlet@meta.data[which(!(singlet@meta.data$Heavy_TRUST %in% c("TRBC1", "TRBC2", "Unknown","","*"))), "Heavy_TRUST"],1,4)
  
  #Transcriptomic
  VDJ.result <- singlet@assays$SCT@data[c(str_detect(rownames(singlet@assays$SCT@data),"IGHG[1-4]|IGH[MDE]|IGHA[1-2]")),]#which(singlet@meta.data$Phénotype %in% "B cells")] #,"TRBC1","TRBC2",
  VDJ.result.finale <- data.frame(row.names = colnames(VDJ.result),count = rownames(VDJ.result) [ apply(VDJ.result, 2, which.max) ])
  
  singlet@meta.data$Heavy <- str_replace_na(singlet@meta.data$Heavy,"Unknown" )
  singlet@meta.data[which(singlet@meta.data$Heavy %in% c("Unknown")), "Heavy"] <- substr(VDJ.result.finale[which(singlet@meta.data$Heavy %in% c("Unknown")),1],1,4)
  

  ## -- Figure -- ##
  load(file=paste0("/home/boris/Documents/lipinskib/Boris_Manon/flinovo/result/",patient,"/R/VDJ.RData"))
  singlet@tools$Clonotype <- Clonotype
  singlet@tools$Type <- Type
  singlet@tools$Isotype <- Isotype
  singlet@tools$V <- V
  singlet@tools$D <- D
  singlet@tools$J <- J
  singlet@tools$Heavy <- Heavy
  singlet@tools$Light <- Light
  singlet@tools$vloupe <- vloupe
  return(singlet)
}
other <- function(singlet){
  singlet@meta.data$Greffe <- as.character(singlet@meta.data$Condition)   # Dissociation Pré/Post greffe
  singlet@meta.data[singlet@meta.data$Condition == "RCHOP", "Greffe" ] <- "Post-Greffe"
  singlet@meta.data[singlet@meta.data$Condition == "Excipient", "Greffe" ] <- "Post-Greffe"
  singlet@meta.data[singlet@meta.data$Condition == "Pré-greffe", "Greffe" ] <- "Pré-Greffe"

  return(singlet)
}
sunburst <- function(singlet){
  met <- data.frame(
    patient = singlet@meta.data$Sample,
    etat = singlet@meta.data$Greffe,
    condition = singlet@meta.data$Condition,
    phénotype = singlet@meta.data$Phénotype,
    sub = singlet@meta.data$Phénotype.fine,
    phase = singlet@meta.data$Phase,
    value = rep(1, length(singlet@meta.data$Sample)), stringsAsFactors = FALSE
  )

  df <- data.frame(ids=character(), labels=character(), parents=character(), values=integer(),stringsAsFactors=FALSE)

  df0 <- met %>% mutate(path = paste(patient,  sep=";")) %>% dplyr::select(path, value)
  df1 <- met %>% mutate(path = paste(patient, etat, sep=";")) %>% dplyr::select(path, value)
  df2 <- met %>% mutate(path = paste(patient, etat, condition,  sep=";")) %>% dplyr::select(path, value)
  df3 <- met %>% mutate(path = paste(patient, etat, condition, phénotype, sep=";")) %>% dplyr::select(path, value)
  df4 <- met %>% mutate(path = paste(patient, etat, condition, phénotype, sub, sep=";")) %>% dplyr::select(path, value)
  df5 <- met %>% mutate(path = paste(patient, etat, condition, phénotype, sub, phase, sep=";")) %>% dplyr::select(path, value)

  for (dfX in c(df0,df1,df2,df3,df4,df5)) {
    xnames <- row.names(table(dfX))
    max <- 0

    for (rows in xnames) {
      orga = strsplit(rows,";")
      if (max < length(orga[[1]])){max = length(orga[[1]])}}

    for (i in 1:max){
      for (rows in xnames) {
        orga = strsplit(rows,";")
        if ((length(orga[[1]]))==i){
          c = ""
          c2 = ""
          for (j in 1:i){c = paste0(c,orga[[1]][j],"-")}
          for (k in 1:i-1){c2 = paste0(c2,orga[[1]][k],"-")}
          c = str_sub(c,1,-2)
          c2 = str_sub(c2,2,-2)
          df[rows,"labels"]= orga[[1]][i]
          df[rows,"ids"]=c
          df[rows,"parents"]=c2
        }
      }
    }
  }
  df <- df[-2,]
  df$values <- c(table(df0)[,1],table(df1)[,1],table(df2)[,1],table(df3)[,1],table(df4)[,1],table(df5)[,1])

  singlet@tools$sunburst <-  df
  return(singlet)
}

diet <- function(singlet){
  singlet <- Seurat::DietSeurat(singlet, counts = F, data = T, scale.data = F, features = rownames(singlet), assays = c("SCT"), dimreducs = c("pca","umap",'tsne'), graphs = NULL )
  return(singlet)
}
sub <- function(singlet, condition, name, patient){
  ## -- Subset -- ##
  singlet <- seurat_subset(singlet,"Condition", condition)
  
  ## -- RNA normalisation -- ##
  singlet <- NormalizeData(singlet, assay = "RNA")
  singlet <- ScaleData(singlet, features = rownames(singlet))
  
  ## -- Normalisation -- ##
  singlet <- normalisation(singlet)
  singlet <- visualisation(singlet)
  
  ## -- Sauvegarde -- ##
  save(singlet, file = paste0("/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_patient/", patient,"/",name,"_", patient,".RData"))
  singlet <- diet(singlet) ; save(singlet, file = paste0("/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/datasets/Patient/",name,"_", patient,".RData"))
  
}