## -- Other -- ##
correlation <- function(singlet){
  ## -- Phénotype -- ##
  load(file = paste0("/home/boris/Documents/analyse/", patient,"/phenotype/phenotype_", patient,".RDs"))
  singlet <- AddMetaData(singlet, labels, col.name = "Phénotype")
  singlet <- AddMetaData(singlet, labels.fine , col.name = "Phénotype.fine")
  
  linear_correlation <- function(ctrl, stim){
    ctrl$stim <- "CTRL" ; stim$stim <- "STIM"
    immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20, k.filter=80)
    immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20, k.weight=80)
    DefaultAssay(immune.combined) <- "integrated"
    immune.combined <- ScaleData(immune.combined, verbose = FALSE)
    immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
    immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
    immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
    immune.combined <- FindClusters(immune.combined, resolution = 0.5)
    
    Idents(immune.combined)<-"Phénotype" ; DefaultAssay(immune.combined) <- "RNA"
    nk.markers <- FindConservedMarkers(immune.combined, ident.1 = "B-cells", grouping.var = "stim", verbose = FALSE)
    
    b.cells <- subset(immune.combined, idents = "B-cells") ; Idents(b.cells) <- "stim"
    return(as.data.frame(log1p(AverageExpression(b.cells, verbose = FALSE)$RNA)))
  }
  Excipient <- seurat_subset(singlet, "Condition", "Excipient")
  RCHOP <- seurat_subset(singlet, "Condition", "RCHOP")
  Pregreffe <- seurat_subset(singlet, "Condition", "Pré-greffe")
  
  avg.b.cells_RE <- linear_correlation(Excipient, RCHOP)
  avg.b.cells_PE <- linear_correlation(Pregreffe, Excipient)
  
  save(avg.b.cells_RE, avg.b.cells_PE, file = paste0("/home/boris/Documents/analyse/", patient,"/correlation/correlation_", patient,".RDs"))
  
  return(singlet)
}




## -- Correlation -- ##
#load(file = paste0("/home/boris/Bureau/Flinovo/result/analyse_patient/",patient,"/correlation/correlation_", patient,".RDs"))
#singlet@tools$avg.b.cells_RE <- avg.b.cells_RE
#singlet@tools$avg.b.cells_PE <- avg.b.cells_PE

## -- 23-genes -- ##
#singlet <- AddModuleScore(singlet,features = list(c("VPREB1","FOXO1","FCRL2","AFF3","TCF4","RASSF6","GADD45A","E2F5","USP44","CXCR4","SEMA4B","EML6","DCAF12","VCL","RGS10","CXCR4","KIAA0040","TAGAP","ORAI2","METRNL","PRDM15","ABCB1","ALDH2","SHISA8")), name="genes23")



cell_to_keep <- function(emb, velo){
  emb=emb
  vel=velo
  n = 100
  cell.colors = NULL
  corr.sigma = 0.05
  xlab = ""
  ylab = ""
  n.cores = 12
  do.par = T
  show.cell = NULL
  cell.border.alpha = 0.3
  cc = NULL
  
  randomize <- FALSE
  if (do.par) {par(mfrow = c(1, 1), mar = c(3.5, 3.5, 2.5, 1.5), mgp = c(2,0.65, 0), cex = 0.85)}
  celcol <- "white"
  if (is.null(show.cell)) {celcol <- cell.colors[rownames(emb)]}
  #plot(emb, bg = celcol, pch = 21, col = ac(1, alpha = cell.border.alpha), xlab = xlab, ylab = ylab)
  em <- as.matrix(vel$current)
  ccells <- intersect(rownames(emb), colnames(em))
  em <- em[, ccells]
  emb <- emb[ccells, ]
  nd <- as.matrix(vel$deltaE[, ccells])
  cgenes <- intersect(rownames(em), rownames(nd))
  nd <- nd[cgenes, ]
  em <- em[cgenes, ]
  
  colDeltaCorSqrt <- function(e, d, nthreads = 1L) {.Call('_velocyto_R_colDeltaCorSqrt', PACKAGE = 'velocyto.R', e, d, nthreads)}
  
  cc <- colDeltaCorSqrt(em, (sqrt(abs(nd)) * sign(nd)), nthreads = n.cores)
  colnames(cc) <- rownames(cc) <- colnames(em)
  diag(cc) <- 0
  
  if (n > nrow(cc)) {n <- nrow(cc)}
  
  tp <- exp(cc/corr.sigma) 
  tp <- t(t(tp)/Matrix::colSums(tp))
  tp <- as(tp, "dgCMatrix")
  tp <- as.data.frame(tp)
  tp <- t(tp)
  dim(tp)
  
  tp <- na.omit(tp)
  dim(tp)

  to_keep <- rownames(tp)
  return(to_keep)
}
standard_normalisation <- function(singlet){
  singlet <- NormalizeData(singlet, normalization.method = "LogNormalize", scale.factor = 10000) 
  singlet <- FindVariableFeatures(singlet, selection.method = "vst", nfeatures = 2000)
  singlet <- ScaleData(singlet, features = rownames(singlet)) 
  return(singlet)
}

## -- Old -- ##
all_visualisation <- function(singlet, feature = VariableFeatures(singlet)){
  singlet <- SCTransform(singlet, method = "glmGamPoi", verbose = F, ncells = 10000, conserve.memory = T) %>%
    RunPCA(verbose = F, approx=F) %>%
    RunUMAP(dims = 1:40, verbose = F) %>%
    FindNeighbors(dims = 1:40, verbose = F) %>%
    FindClusters(verbose = F) %>%
    RunTSNE(dims = 1:40, verbose = F)
  return(singlet)
}
old_visualisation <- function(singlet){
  singlet <- subset(singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
  singlet <- NormalizeData(singlet, normalization.method = "LogNormalize", scale.factor = 10000)
  singlet <- FindVariableFeatures(singlet, selection.method = "vst", nfeatures = 2000)
  singlet <- ScaleData(singlet, features = rownames(singlet))                                               # Scale data :singlet@assays$RNA@scale.data, singlet[["RNA"]]@scale.data
  singlet <- RunPCA(singlet, features = VariableFeatures(singlet), ndims.print = 1:10, nfeatures.print = 10)# Reduction dimension
  singlet <- FindNeighbors(singlet, reduction = "pca", dims = 1:40, compute.SNN = T)
  singlet <- FindClusters(singlet, resolution = 0.5)                                                        # head(Idents(singlet), 10)
  singlet <- RunUMAP(singlet, reduction = "pca", dims = 1:40)
  singlet <- RunTSNE(singlet, reduction = "pca", dims = 1:40)
  return(singlet)
}


plot_dimension <- function(singlet){
  pca1 <- Seurat::DimPlot(object = singlet, label.size = 5, pt.size = 0.5, reduction = 'pca', label = TRUE, dims = c(1,2)) & theme(title = element_text(size=20),legend.position = "top",legend.title = element_text(size=10),legend.text = element_text(size=10)) & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 6))) & xlab(label = paste0("PCA 1 : ", round(Seurat::Stdev(singlet[["pca"]])[1],2), " %")) & ylab(label = paste0("PCA 2 : ", round(Seurat::Stdev(singlet[["pca"]])[2],2), " %"))
  umap1 <- Seurat::DimPlot(object = singlet, label.size = 5, pt.size = 0.5, reduction = 'umap', label = TRUE) & theme(legend.position = "top") & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4)))
  tsne1 <- Seurat::DimPlot(object = singlet, label.size = 5, pt.size = 0.5, reduction = 'tsne', label = TRUE) & theme(legend.position = "top") & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4)))
  return(cowplot::plot_grid(nrow = 1,ncol = 3, pca1, umap1, tsne1))
}
plot_presentation <- function(singlet){
  pca1 <- Seurat::DimPlot(object = singlet, label.size = 5, pt.size = 0.5, reduction = 'pca', label = TRUE) & Seurat::NoLegend() & Seurat::NoAxes()
  umap1 <- Seurat::DimPlot(object = singlet, label.size = 5, pt.size = 0.5, reduction = 'umap', label = TRUE) & Seurat::NoLegend() & Seurat::NoAxes()
  tsne1 <- Seurat::DimPlot(object = singlet, label.size = 5, pt.size = 0.5, reduction = 'tsne', label = TRUE) & Seurat::NoLegend() & Seurat::NoAxes()
  return(cowplot::plot_grid(nrow = 1,ncol = 3, pca1, umap1, tsne1))
}


gene_subset <- function(singlet, gene, expression){
  expr <- FetchData(object = singlet, vars = gene)
  sub_singlet <- singlet[, which(x = expr > expression)]
  sub_singlet <- visualisation(sub_singlet)
  return(sub_singlet)
}
QC_subset <- function(singlet, maximum_sub, percent_mt_sub){
  if(maximum_sub==""){maximum_sub = 10000}
  if(percent_mt_sub==""){percent_mt_sub = 50}
  
  ## -- QC filtres-- ##
  Idents(singlet)<-"seurat_clusters"
  sub_singlet <- subset(singlet, subset = nFeature_RNA > 200 & nFeature_RNA < maximum_sub & percent.mt < percent_mt_sub)            # QC Filter : tester 15%
  sub_singlet <- visualisation(sub_singlet)
  
  return(sub_singlet)
}

## -- Merge -- ##
numeric_merge <- function(singlet){
  for (metadata in colnames(singlet@meta.data)) {
    if (metadata %in% singlet@tools$hallmarks){
      singlet@meta.data[[metadata]]<- as.numeric(singlet@meta.data[[metadata]])
    }
  }
  return(singlet)
}

## -- Subset -- ##
sub_merge <- function(singlet, condition, name){
  ## -- Load et subset -- ##
  singlet <- seurat_subset(singlet,"Condition", condition); gc(); gc(); gc(); gc();
  
  ## -- Intergration -- ##
  singlet <- integration(singlet)
  
  ## -- Visualisation -- ##
  feature <- VariableFeatures(singlet)[which(!str_detect(VariableFeatures(singlet), "IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]"))];
  singlet <- RunPCA(singlet, verbose = F, approx=F, feature = feature) %>% RunUMAP(verbose = T, umap.method = "umap-learn", feature = feature) %>% FindNeighbors(feature = feature, verbose = F) %>% FindClusters(resolution = 0.6, verbose = F) %>% RunTSNE(verbose = F, perplexity = 100, feature = feature) #%>% RunVelocity(deltaT = 1, kCells = 50, fit.quantile = 0.02, na.rm = T, ncores = 12)
  
  ## -- DE : Réponse RC/RP + RCHOP/Excipient -- ##
  DE <- singlet
  DE@assays$RNA@data <- DE@assays$RNA@data[rownames(DE@assays$RNA@data)[which(!str_detect(rownames(DE@assays$RNA@data), "IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]"))],];
  DE@assays$RNA@scale.data <- DE@assays$RNA@scale.data[rownames(DE@assays$RNA@scale.data)[which(!str_detect(rownames(DE@assays$RNA@scale.data), "IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]"))],];
  DE@assays$RNA@counts <- DE@assays$RNA@counts[rownames(DE@assays$RNA@counts)[which(!str_detect(rownames(DE@assays$RNA@counts), "IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]"))],];
  
  Idents(DE) <- "Reponse"  ; singlet@tools$DE_Reponse <- FindMarkers(DE, assay = "RNA", slot = "data", ident.1 = "RC", ident.2 = "RP", test.use = "negbinom")
  Idents(DE) <- "Condition"; singlet@tools$DE_Condition <- FindMarkers(DE, assay = "RNA", slot = "data", ident.1 = "RCHOP", ident.2 = "Excipient", test.use = "negbinom")
  
  ## -- Sauvegarde -- ##
  save(singlet, file = paste0("/home/boris/Bureau/Flinovo/result/analyse_meta/",name,".RData"))
  
  singlet <- diet_merge(singlet)
  save(singlet, file = paste0("/home/boris/Bureau/scShiny/datasets/All/",name,".RData"))
  
}
diet_merge <- function(singlet){
  return(Seurat::DietSeurat(singlet, counts = F, data = T, scale.data = F, features = rownames(singlet), assays = "integrated", dimreducs = c("pca","umap",'tsne'), graphs = NULL ))
}



enrichissement <- function(singlet){
  GS <- getGeneSets(library = "H")
  ES <- enrichIt(obj = singlet, gene.sets = GS, groups = 1000, cores = 10)
  names(ES) <- str_replace_all(names(ES), "HALLMARK_", "")
  save(ES, file = paste0("/home/boris/Bureau/Flinovo/result/analyse_patient/", patient,"/enrichissement/enrichissement_", patient,".RDs"))
}

## -- Enrichissement -- ##
#load(file = paste0("/home/boris/Bureau/Flinovo/result/analyse_patient/",patient,"/enrichissement/enrichissement_", patient,".RDs"))
#singlet <- AddMetaData(singlet, ES)
#singlet@tools$hallmarks <- names(ES)

## -- Enrichissement -- ##
#GS <- getGeneSets(library = "H")
#ES <- enrichIt(obj = singlet, gene.sets = GS, groups = 1000, cores = 12)
#names(ES) <- str_replace_all(names(ES), "HALLMARK_", "")
#singlet <- AddMetaData(singlet, ES)
#singlet@tools$hallmarks <- names(ES)
#save(ES, file = paste0("/home/boris/Documents/analyse/", patient,"/enrichissement/",name,"_enrichissement_", patient,".RDs"))

# KEGG <- DEenrichRPlot(singlet, assay ="RNA", slot = "data", ident.1 = "RCHOP", ident.2 = "Excipient", test.use = "negbinom" ,max.cells.per.ident = Inf, balanced=T, p.val.cutoff=0.05, return.gene.list=T, num.pathway = 15, enrich.database = "KEGG_2021_Human", max.genes = Inf, logfc.threshold = 0.1)
# GO_Biological <- DEenrichRPlot(singlet, assay ="RNA", slot = "data", ident.1 = "RCHOP", ident.2 = "Excipient", test.use = "negbinom", max.cells.per.ident = Inf, balanced=T, p.val.cutoff=0.05, return.gene.list=T, num.pathway = 15, enrich.database = "GO_Biological_Process_2021", max.genes = Inf, logfc.threshold = 0.1)
# GO_Cellular <- DEenrichRPlot(singlet, assay ="RNA", slot = "data", ident.1 = "RCHOP", ident.2 = "Excipient", test.use = "negbinom", max.cells.per.ident = Inf, balanced=T, p.val.cutoff=0.05, return.gene.list=T, num.pathway = 15, enrich.database = "GO_Cellular_Component_2021", max.genes = Inf, logfc.threshold = 0.1)
# GO_Molecular <- DEenrichRPlot(singlet, assay ="RNA", slot = "data", ident.1 = "RCHOP", ident.2 = "Excipient", test.use = "negbinom", max.cells.per.ident = Inf, balanced=T, p.val.cutoff=0.05, return.gene.list=T, num.pathway = 15, enrich.database = "GO_Molecular_Function_2021", max.genes = Inf, logfc.threshold = 0.1)
#save(DE_RE, DE_PE, KEGG, GO_Biological, GO_Cellular, GO_Molecular, file = paste0("/home/boris/Bureau/Flinovo/result/analyse_patient/", patient,"/DE_", patient,".RDs"))

# singlet@tools$KEGG <- KEGG
# singlet@tools$GO_Biological <- GO_Biological
# singlet@tools$GO_Cellular <- GO_Cellular
# singlet@tools$GO_Molecular <- GO_Molecular



## -- Signature -- ##
#singlet <- AddModuleScore(singlet,features = list(c("VPREB1","FOXO1","FCRL2","AFF3","TCF4","RASSF6","GADD45A","E2F5","USP44","CXCR4","SEMA4B","EML6","DCAF12","VCL","RGS10","CXCR4","KIAA0040","TAGAP","ORAI2","METRNL","PRDM15","ABCB1","ALDH2","SHISA8")), name="genes23")
#Idents(singlet) <- "seurat_clusters"












# aggregateBioVar ---------------------------------------------------------
### Enrichissement fonctionnel à partir des DEGs : GO & KEGG
```{r echo=F}
# go$P.Up <- signif(go$P.Up,3)
# go$P.Down <- signif(go$P.Down,3)
# 
# ko$P.Up <- signif(ko$P.Up,3)
# ko$P.Down <- signif(ko$P.Down,3)
# 
# 
# GO.Up = as.data.frame(topGO(go,sort="up", number = Inf ))
# GO_custom.Up <- data.frame(
#   Term = rep(paste(GO.Up$Ont, ":", GO.Up$Term, "\n", GO.Up$N, "genes - FDR =", GO.Up$P.Up),2),
#   Regulation = c(rep("up",length(GO.Up$Term)),rep("down",length(GO.Up$Term))),
#   Count = c(GO.Up$Up,GO.Up$Down), 
#   pValueUp = c(GO.Up$P.Up, GO.Up$P.Up),
#   label_ypos = c(GO.Up$Up-3,GO.Up$Up+GO.Up$Down+2)
# )
# GO_custom.Up <- GO_custom.Up[1:20,]
# ggplot(data=GO_custom.Up, aes(x=reorder(Term, -pValueUp), y=Count, fill=Regulation), color = Regulation) +
#   geom_bar(stat="identity", width = 0.7)+ xlab("GO term") + ylab("Number of DEGs up and down inside the up-regulated term.") + ggtitle("Up-regulated pathways")  +
#   geom_text(aes(y=label_ypos, label=Count), color="black", size=6) + scale_fill_brewer(palette="Paired") + theme_minimal() + 
#   theme(aspect.ratio = 3/2, axis.text=element_text(size=10), axis.title=element_text(size=15), plot.title = element_text(size=15)) + coord_flip()
# 
# 
# GO.Down = as.data.frame(topGO(go,sort="down"), number = Inf)
# GO_custom.Down <- data.frame(
#   Term = rep(paste(GO.Down$Ont, ":", GO.Down$Term, "\n", GO.Down$N, "genes - FDR =", GO.Down$FDR.Down),2),
#   Regulation = c(rep("up",length(GO.Down$Term)),rep("down",length(GO.Down$Term))),
#   Count = c(GO.Down$Up,GO.Down$Down), 
#   pValueDown = c(GO.Down$P.Down, GO.Down$P.Down),
#   label_ypos = c(GO.Down$Up-3,GO.Down$Up+GO.Down$Down+2)
# )
# GO_custom.Down <- GO_custom.Down[1:20,]
# ggplot(data=GO_custom.Down, aes(x=reorder(Term, -pValueDown), y=Count, fill=Regulation), color = Regulation) +
#   geom_bar(stat="identity", width = 0.7) + xlab("GO term") + ylab("Number of DEGs up and down inside the down-regulated term.") + ggtitle("Down-regulated pathways") +
#   geom_text(aes(y=label_ypos, label=Count),  color="black", size=6) + scale_fill_brewer(palette="Paired") + theme_minimal() + 
#   theme(aspect.ratio = 3/2, axis.text=element_text(size=15), axis.title=element_text(size=15), plot.title = element_text(size=22)) + coord_flip()


```


```{r}
# intersect(rownames(singlet.result), signature)
# 
# Singature_RCHOP=rownames(subj_dds_transf)[1:15]
# Singature_RCHOP1=c("PHPT1","HIST1H2BK", "RPS27L","SRSF3","TRIM22","BBC3","FDXR","PSMB4","PPIA","HSPA1A","BAX","ZMAT3","ZFAS1","CCNG1","COX8A","MDM2","SNHG8","EIF3E")
# Singature_RCHOP2=c("PHPT1","HIST1H2BK", "RPS27L","SRSF3","TRIM22","BBC3","FDXR","PSMB4","PPIA","AEN","HSPA1A","RACK1","MDM2","SNHG8","EIF3E","BAX","ZMAT3","ZFAS1","CCNG1","COX8A")
# Singature_RCHOP3=c("PHPT1","HIST1H2BK", "RPS27L","SRSF3","TRIM22","BBC3","FDXR","PSMB4","AEN","HSPA1A","BAX")#totest auc!!!!
# Singature_RCHOP=c("TRIM22","AEN","PHPT1","FDXR","BAX","MDM2","TRIAP1","BBC3","DDB2","CCNG1","TNFSF8","GADD45A", "CDKN1A" , "PHLDA3" , "ZMAT3","CD70")
# Singature_RCHOP0 <- unique(c(Singature_RCHOP1,Singature_RCHOP2,Singature_RCHOP3,Singature_RCHOP))
# 
# singlet <- AddModuleScore(singlet,features = list(Singature_RCHOP), ctrl=500, name="RCHOP_RNA_Score", assay = "RNA")
# Seurat::VlnPlot(singlet, idents = c("Excipient","RCHOP"), features = c("RCHOP_RNA_Score1"), ncol = 1, assay = 'SCT', slot = "data") & theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
# Seurat::VlnPlot(singlet, idents = c("Excipient","RCHOP"), features = c("RCHOP_AUC_Score"), ncol = 1, assay = 'RNA', slot = "count") & theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
```
