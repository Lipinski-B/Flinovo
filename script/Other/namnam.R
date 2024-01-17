source(file = "/home/boris/Bureau/scShiny/functions.R")
library(celldex)
library(SingleR)
library(escape)
library(stringr)
library(dplyr)
library(future)
library(DESeq2)
library(cowplot)
library(Seurat)
library(dittoSeq)
library(DT)
library(stringr)
singlet_namanm <- function(path, patient){
  rwa.mRNA <- Read10X(data.dir = path) 
  cso <- CreateSeuratObject(counts = rwa.mRNA, project = patient, min.cells = 3, min.features = 200)
  umis <- GetAssayData(object = cso, slot = "counts")
  singlet <- CreateSeuratObject(counts = umis, assay = "RNA", project = patient)
  return(singlet)
}
metadata_namnam <- function(singlet, path){  
  ##### -- Identification des cellules -- #####
  BD <- celldex::BlueprintEncodeData()
  results.blueprint <- SingleR(test = as.SingleCellExperiment(singlet), ref = BD, labels = BD$label.main)
  results.blueprint.fine <- SingleR(test = as.SingleCellExperiment(singlet), ref = BD, labels = BD$label.fine)
  singlet$Phénotype <- results.blueprint$labels 
  singlet$Phénotype.fine <- results.blueprint.fine$labels 
  
  
  ##### -- Enrichissement des gènes -- ##### 
  GS <- getGeneSets(library = "H")
  ES <- enrichIt(obj = singlet, gene.sets = GS, groups = 1000, cores = 12)
  names(ES) <- str_replace_all(names(ES), "HALLMARK_", "")
  singlet <- AddMetaData(singlet, ES)
  singlet@tools$hallmarks <- names(ES)
  
  
  ##### -- Phase cycle -- #### 
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  singlet <- CellCycleScoring(singlet, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)    # Attribution à chaque cellule d'un état de phase
  #singlet <- ScaleData(singlet, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(singlet)) # Correction de la matrice d'expression dépendemment du cycle cellulaire
  Idents(singlet)<-"seurat_clusters"
  
  
  ##### -- VDJ -- ##### 
  # Cellranger
  bcr <- read.csv(paste0(path, "filtered_contig_annotations.csv"))  # bcr$barcode <- gsub("-1", "", bcr$barcode)                                   # Remove the -1 at the end of each barcode: VERIFIER SI IL Y A LES TIRETS DANS OBJECT SINGLET!!!
  bcr <- bcr[!duplicated(bcr$barcode), ]                                                                  # Subsets so only the first line of each barcode is kept,as each entry for given barcode will have same clonotype.
  bcr <- bcr[,c("barcode", "raw_clonotype_id","chain","v_gene","d_gene","j_gene","c_gene","cdr3")] ; names(bcr)[names(bcr) == "raw_clonotype_id"] <- "clonotype_id"
  clonotype <- read.csv(paste0(path,"clonotypes.csv"))              # Clonotypes
  bcr <- merge(bcr, clonotype[, c("clonotype_id", "cdr3s_aa")])                 # Selection de ceux pour lesquels on a une séquence cdr3
  bcr <- bcr[, c(2,1,3,4,5,6,7,8,9)]                                            # Mettre les codes barres comme nom de la première colonne
  rownames(bcr) <- bcr[,1]
  bcr[,1] <- NULL
  bcr$clonotype_id <- as.numeric(str_replace(unlist(bcr$clonotype_id,"clonotype",use.names=F), "clonotype", ""))
  singlet <- AddMetaData(object=singlet, metadata = bcr)                        #, col.name = colnames(as.data.frame(bcr)))   # Add to the Seurat object's metadata.
  
  # Figure
  #load(file=path)
  singlet@tools$Clonotype <- Clonotype
  singlet@tools$Type <- Type
  singlet@tools$Isotype <- Isotype
  singlet@tools$V <- V
  singlet@tools$D <- D
  singlet@tools$J <- J
  singlet@tools$Heavy <- Heavy
  singlet@tools$Light <- Light
  singlet@tools$vloupe <- vloupe
  
  
  ##### -- Ohter -- #####
  singlet@tools$meta_variable <- c("seurat_clusters", "Phénotype", "orig.ident", "Phase")#, "clonotype_id" "BCL2_K22K", "BCL2_L23L", "CD79B_Y696H")# c("seurat_clusters", "Condition", "Phénotype", "Phase", "K29Q", "L37M", "M11I", "pGln45")
  singlet@meta.data$nFeature_HTO <- NULL
  singlet@meta.data$HTO_secondID <- NULL
  singlet@meta.data$HTO_classification <- NULL
  singlet@meta.data$RNA_snn_res.0.5 <- NULL
  singlet@meta.data$hash.ID <- NULL
  
  #colnames(singlet@meta.data)[which(colnames(singlet@meta.data)=="SingleR.calls")] <- "Phénotype"
  
  
  ##### -- Sunburst -- #### 
  make_sunburst_data <- function(met){
    df <- data.frame(ids=character(), labels=character(), parents=character(), values=integer(),stringsAsFactors=FALSE)
    
    df2 <- met %>% mutate(path = paste(patient,  sep=";")) %>% dplyr::select(path, value)
    df3 <- met %>% mutate(path = paste(patient, phénotype, sep=";")) %>% dplyr::select(path, value)
    df4 <- met %>% mutate(path = paste(patient, phénotype, sub, sep=";")) %>% dplyr::select(path, value)
    df5 <- met %>% mutate(path = paste(patient, phénotype, sub, phase, sep=";")) %>% dplyr::select(path, value)
    
    for (dfX in c(df2,df3,df4,df5)) {
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
    df$values <- c(table(df2)[,1],table(df3)[,1],table(df4)[,1],table(df5)[,1])
    return(df)
  }
  met <- data.frame(
    patient = singlet@meta.data$orig.ident, 
    phénotype = singlet@meta.data$Phénotype, 
    sub = singlet@meta.data$Phénotype.fine, 
    phase = singlet@meta.data$Phase,
    value = rep(1, length(singlet@meta.data$orig.ident)), stringsAsFactors = FALSE
  )
  singlet@tools$sunburst <- make_sunburst_data(met)
  
  return(singlet)
}
count <- function(column){
  names<-c();number<-c()
  for (k in 1:length(column)) {
    for (i in 1:length(names(table(vloupe[column[k]])))) {
      if(names(table(vloupe[column[k]]))[i]!="NA"){
        names <- c(names, names(table(vloupe[column[k]]))[i])
        cell <- subset(vloupe, vloupe[column[k]] == names(table(vloupe[column[k]]))[i])
        number <- c(number, sum(cell$frequency))
      }
    }
  }
  return(c(list(names), list(number)))
}


setwd(dir = "/home/boris/Documents/lipinskib/narimene/")
siege <- c("FL1085",  "FL120316",  "FL1214",  "FL1481")
patient <- siege[4]

for (patient in siege) {
  load(file=paste0("/home/boris/Documents/analyse/namnam/",patient,".RData"))
  
  singlet@assays[["SCT"]]@scale.data <- subset(singlet@assays[["SCT"]]@scale.data, rownames(singlet@assays[["SCT"]]@scale.data) %in% c(rownames(singlet@tools$DE_PE)[1:50], rownames(singlet@tools$DE_RE)[1:50]))
  singlet@assays[["RNA"]]@data <-subset(as.matrix(singlet@assays[["RNA"]]@data), rownames(singlet@assays[["RNA"]]@data) %in% c(rownames(singlet@tools$DE_PE)[1:50], rownames(singlet@tools$DE_RE)[1:50]))
  singlet <- DietSeurat(singlet, counts = FALSE, data = T, scale.data = T,features = NULL, assays = NULL, dimreducs = c("pca","umap",'tsne'), graphs = NULL )
  #singlet@assays[["HTO"]] <- list() ; 
  singlet@assays[["RNA"]]@counts <- as.matrix(0)
  singlet@assays[["RNA"]]@scale.data <- as.matrix(0)
  
  singlet@assays[["SCT"]]@data <- as.matrix(singlet@assays[["SCT"]]@data)
  singlet@assays[["SCT"]]@data <- subset(singlet@assays[["SCT"]]@data, rownames(singlet@assays[["SCT"]]@data) %in% VariableFeatures(singlet)[1:50])
  
  save(singlet, file = paste0("/home/boris/Bureau/scShiny/www/", patient,"/", patient,".RData"))
  
}
for (patient in siege) {
  singlet <- singlet_namanm(paste0(patient,"/Results/",patient,"_CellrangerCount/outs/filtered_feature_bc_matrix/"), patient)
  singlet <- visualisation(singlet)
  singlet@tools$mitochondrie_all <- VlnPlot(singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  vloupe <- as.data.frame(read.table(paste0("/home/boris/Téléchargements/",patient,"_vloupe-clonotypes.csv"), sep = ",", header = T))
  vloupe$clonotype_id <- as.numeric(str_replace(vloupe$clonotype_id,"clonotype",""))
  
  vloupe$lourde <- ifelse(vloupe$igh_c_genes!="", "IGH", "NA")
  vloupe$legere[vloupe$igl_c_genes=="" & vloupe$igk_c_genes==""] <- "NA"
  vloupe$legere[vloupe$igl_c_genes=="" & vloupe$igk_c_genes!=""] <- "IGK"
  vloupe$legere[vloupe$igl_c_genes!="" & vloupe$igk_c_genes==""] <- "IGL"
  vloupe$legere[vloupe$igl_c_genes!="" & vloupe$igk_c_genes!=""] <- "IGK-IGL"
  vloupe$type <- paste0(vloupe$lourde,"-",vloupe$legere)
  
  vloupe$igh_c_genes <- ifelse(vloupe$igh_c_genes!="", vloupe$igh_c_genes, "NA")
  vloupe$V_lourde <- ifelse(vloupe$igh_v_genes!="", vloupe$igh_v_genes, "NA")
  vloupe$D_lourde <- ifelse(vloupe$igh_d_genes!="", vloupe$igh_d_genes, "NA")
  vloupe$J_lourde <- ifelse(vloupe$igh_j_genes!="", vloupe$igh_j_genes, "NA")
  
  vloupe$V_legere <- ifelse(vloupe$igk_v_genes!="", vloupe$igk_v_genes, ifelse(vloupe$igl_v_genes!="", vloupe$igl_v_genes, "NA"))
  vloupe$J_legere <- ifelse(vloupe$igk_j_genes!="", vloupe$igk_j_genes, ifelse(vloupe$igl_j_genes!="", vloupe$igl_j_genes, "NA"))
  vloupe$heavy <-str_sub(vloupe$igh_c_genes, 1, 4)
  
  ### sauvegarde des objets
  Clonotype <- c(list(vloupe$clonotype_id[1:5]), list(vloupe$frequency[1:5]))
  Isotype <- count("igh_c_genes")
  Type <- count("type")
  V <- count(c("V_lourde","V_legere"))
  D <- count("D_lourde")
  J <- count(c("J_lourde","J_legere"))
  Light <- count("legere")
  Heavy <- count("heavy")
  
  singlet <- metadata_namnam(singlet, paste0(patient,"/Results/",patient,"_CellrangerVDJ/outs/"))
  
  save(singlet, file=paste0("/home/boris/Documents/analyse/namnam/",patient,".RData"))
  
  singlet <- DietSeurat(singlet, counts = F, data = T, scale.data = F,features = NULL, assays = NULL, dimreducs = c("pca","umap",'tsne'), graphs = NULL )
  singlet@assays[["RNA"]]@counts <- as.matrix(0)
  singlet@assays[["RNA"]]@scale.data <- as.matrix(0)
  singlet@assays[["RNA"]]<- list()
  
  singlet@assays[["SCT"]]@data <-subset(as.matrix(singlet@assays[["SCT"]]@data), rownames(singlet@assays[["SCT"]]@data) %in% VariableFeatures(singlet)[1:100])
  #singlet@assays$SCT@data<- as.matrix(0)
  
  save(singlet, file = paste0("/home/boris/Bureau/scShiny/www/",patient,"/",patient,".RData"))
  
  rm(singlet) ; gc() ; gc() ; gc()

}


Seurat::DimPlot(object = singlet, label.size = 0.0, pt.size = 2, reduction = 'pca') & 
  theme(title = element_text(size=20),legend.position = "top",legend.title = element_text(size=10),legend.text = element_text(size=10)) & 
  guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 6))) 


load(file = paste0("/home/boris/Documents/analyse/namnam/FL1085.RData"))
load(file = paste0("/home/boris/Bureau/scShiny/www/FL1085/FL1085.RData"))
