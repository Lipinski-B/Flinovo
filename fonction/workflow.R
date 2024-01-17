## -- Worflow -- ##
seurat_object <- function(patient){
  ##### -- Matrice mRNA -- #####
  rwa.mRNA <- Read10X(data.dir = paste0(patient,"/mRNA/CellrangerCount/outs/filtered_feature_bc_matrix/"))
  cso <- CreateSeuratObject(counts = rwa.mRNA, project = patient, min.cells = 3, min.features = 200)
  umis <- GetAssayData(object = cso, slot = "counts")
  
  ##### -- Matrice HTO -- #####
  raw.hto <- Read10X(paste0(patient,"/HTO/result/umi_count/"), gene.column = 1)
  colnames(raw.hto) <- paste0(colnames(raw.hto),"-1")
  hto <- raw.hto[c(1:3),]                                                       # Suppression des séquences unmapped : rownames(raw.hto)
  rownames(hto) <- c("Pré-greffe","Excipient","RCHOP")
  joint.bcs <- intersect(colnames(umis),colnames(hto))                          # Sélection des cellules avec barcode commun HTO / mRNA
  umis <- umis[, joint.bcs]                                                     # Sélection des lignes qui correspondent aux cellules en commun
  hto <- as.matrix(hto[, joint.bcs])
  
  ##### -- Object seurat -- #####
  hashtag <- CreateSeuratObject(counts = umis, assay = "RNA", project = patient)
  hashtag[["HTO"]] <- CreateAssayObject(counts = hto)                           # Ajoute des données HTO comme un nouvel assay indépendant du mRNA
  hashtag <- NormalizeData(hashtag, assay = "HTO", normalization.method = "CLR")
  hashtag <- HTODemux(hashtag, assay = "HTO", positive.quantile = 0.99)         # Association : cellules / échantillons
  Idents(hashtag)<- 'HTO_classification.global' 
  singlet <- subset(hashtag, idents = "Singlet")
  
  ##### -- Clean seurat -- #####
  colnames(singlet@meta.data)[which(colnames(singlet@meta.data)=="HTO_maxID")] <- "Condition"
  colnames(singlet@meta.data)[which(colnames(singlet@meta.data)=="orig.ident")] <- "Sample"
  singlet@meta.data$HTO_secondID <- singlet@meta.data$HTO_margin <- singlet@meta.data$HTO_classification <- singlet@meta.data$HTO_classification.global <- singlet@meta.data$hash.ID <- NULL
  singlet <- subset(singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #& nCount_RNA > 2100 &   & percent.mt < 5
  
  ##### -- RNA normalisation -- #####
  singlet <- NormalizeData(singlet, assay = "RNA")
  singlet <- ScaleData(singlet, features = rownames(singlet))
  
  ##### -- Pourcentage -- #####
  singlet[["percent.mt"]] <- PercentageFeatureSet(singlet, pattern = "^MT-")
  singlet[["percent.rb"]] <- PercentageFeatureSet(singlet, pattern = "^RP[SL]")
  singlet[["percent.ig"]] <- PercentageFeatureSet(singlet, pattern = "^IG")
  
  return(singlet)
}
normalisation <- function(singlet){
  singlet <- SCTransform(singlet, verbose = F, ncells = 10000, conserve.memory = T, vst.flavor = "v2")
  return(singlet)
}
visualisation <- function(singlet){
  singlet <- RunPCA(singlet, verbose = F, approx=F) %>%
    RunUMAP(dims = 1:40, verbose = F, umap.method	= "umap-learn") %>%
    FindNeighbors(dims = 1:40, verbose = F) %>%
    FindClusters(resolution = 0.6, verbose = F) %>% 
    RunTSNE(dims = 1:20, verbose = F, perplexity = 10 )
  return(singlet)
}
seurat_subset <- function(singlet, item, sub_item){
  Idents(singlet) <- item
  sub_singlet <- subset(singlet, idents = sub_item)
  return(sub_singlet)
}
metadata_merge <- function(singlet){
  for (metadata in colnames(singlet@meta.data)) {
    Idents(singlet)<- metadata
    singlet<- StashIdent(singlet, save.name = metadata)
  }
  return(singlet)
}
noig <- function(DE){
  DE@assays$RNA@data <- DE@assays$RNA@data[rownames(DE@assays$RNA@data)[which(!str_detect(rownames(DE@assays$RNA@data), "IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]"))],];
  DE@assays$RNA@scale.data <- DE@assays$RNA@scale.data[rownames(DE@assays$RNA@scale.data)[which(!str_detect(rownames(DE@assays$RNA@scale.data), "IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]"))],];
  DE@assays$RNA@counts <- DE@assays$RNA@counts[rownames(DE@assays$RNA@counts)[which(!str_detect(rownames(DE@assays$RNA@counts), "IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]"))],];
  return(DE)
}
seurat_DE <- function(DE, Condition, ID1, ID2, FC, assay){
  annotations <- read.csv(paste0(path,"document/annotation.csv"))
  Idents(DE) <- Condition
  
  if(assay=="RNA"){
    result <- FindMarkers(DE, ident.1 = ID1, ident.2 = ID2, logfc.threshold = FC, assay = assay, slot = "data", test.use = "negbinom")
  } else {
    DE <- PrepSCTFindMarkers(DE)
    result <- FindMarkers(DE, ident.1 = ID1, ident.2 = ID2, logfc.threshold = FC, assay = assay, verbose = F)
  }
  result <- result %>% rownames_to_column(var="gene") %>% left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% as_tibble(result) 
  return(result)
}
