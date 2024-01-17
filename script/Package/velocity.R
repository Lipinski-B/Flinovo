# Velocity ----------------------------------------------------------------
source(file = "/home/boris/Bureau/Flinovo/fonction/velocity.R")
## -- Seurat independent -- ##
for (patient in c("FL06G1206")) { #siege
  bam <- velocity(patient)
  for(condition in c("RCHOP","Excipient","Pré-greffe")){
    bm <- seurat_subset(bam, "Condition", condition)
    bm <- SCTransform(bm, assay = "spliced", verbose = F, ncells = 1000, conserve.memory = T) %>% RunPCA(verbose = F) %>% FindNeighbors(dims = 1:20, verbose = F) %>% FindClusters(verbose = F) %>% RunUMAP(dims = 1:20, verbose = F) %>% RunTSNE(dims = 1:20, verbose = F) %>% RunVelocity(deltaT = 1, kCells = 50, fit.quantile = 0.02, na.rm = T)
    bm <- Seurat::DietSeurat(bm, counts = F, data = T, scale.data = F, features = rownames(bm)[1] , assays = "SCT", dimreducs = c('pca','umap','tsne'), graphs = NULL )
    save(bm, file = paste0("/home/boris/Bureau/scShiny/datasets/Velocity/", patient, "_", condition, "_velocyto.RData"))
  }
}


## -- Seurat dependent -- ##
for (patient in "FL06G1206") {
  ldat <- ReadVelocity(file = paste0("/home/boris/Documents/lipinskib/Boris_Manon/flinovo/result/",patient,"/mRNA/CellrangerCount/velocyto/CellrangerCount.loom"))
  bm <- as.Seurat(x = ldat)
  nom <- c() ; for (i in 1:length(colnames(bm))){nom <- c(nom,paste0(str_split(str_split(colnames(bm)[i], ":")[[1]][2],"x")[[1]][1],"-1"))}
  bm <- RenameCells(bm,new.names = nom)
  #bm <- subset(bm, cells = colnames(singlet))
  
  for(condition in c("Pré-greffe-Excipient", "All","Post-greffe","Pré-greffe")){#
    for(phénotype in c("All","BT","B")){
      load(file = paste0("/home/boris/Bureau/Flinovo/result/analyse_patient/",patient,"/",condition,"_",phénotype,"_",patient,".RData")) 
      
      bam <- subset(bm, cell = colnames(singlet))
      
      singlet[['spliced']] <- bam[['spliced']]
      singlet[['unspliced']] <- bam[['unspliced']]
      singlet[['ambiguous']] <- bam[['ambiguous']]
      singlet <- RunVelocity(singlet, deltaT = 1, kCells = 50, fit.quantile = 0.02, na.rm = T)
      
      save(singlet, file = paste0("/home/boris/Bureau/Flinovo/result/analyse_patient/", patient,"/",condition,"_",phénotype,"_", patient,".RData"))
      
      singlet <- diet(singlet)
      save(singlet, file = paste0("/home/boris/Bureau/scShiny/datasets/Patient/",condition,"_",phénotype,"_", patient,".RData"))
    }
  }
}

## -- Visualisation -- ##
emb = Embeddings(object = singlet, reduction = "umap")
velo <- Tool(object = singlet, slot = "RunVelocity") ;
velo$current <- velo$current[,which(!is.na(colSums(velo$current)))] ;
velo$projected <- velo$projected[,which(!is.na(colSums(velo$projected)))]
#sub <- cell_to_keep(emb, velo)
#singlet <- subset(singlet, cell = sub)

Idents(singlet) <- "Phénotype.fine"
velo <- Tool(object = singlet, slot = "RunVelocity") ; velo$current <- velo$current[,which(!is.na(colSums(velo$current)))] ; velo$projected <- velo$projected[,which(!is.na(colSums(velo$projected)))]
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = singlet))) ; names(x = ident.colors) <- levels(x = singlet)
cell.colors <- ident.colors[Idents(object = singlet)] ; names(x = cell.colors) <- colnames(x = singlet)

show.velocity.on.embedding.cor(emb = Embeddings(object = subset(singlet, nFeature_spliced > 50), reduction = "umap"), vel = velo, n = 200, n.cores = 10, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.2), cex = 1.5, arrow.scale = 0.3,
                               show.grid.flow = TRUE, min.grid.cell.mass = 1.5, grid.n = 100, arrow.lwd = 1.5, do.par = F, cell.border.alpha = 0.1)

DoHeatmap(singlet, features = VariableFeatures(singlet)[1:500], group.by = "Condition") + NoLegend()
Seurat::DimPlot(object = singlet, label.size = 5, pt.size = 1, reduction = "umap", group.by = "Condition") & Seurat::NoLegend() & xlab(label = paste0(input$Reduction, " / PCA 1 : ", round(Seurat::Stdev(r$dataset[["pca"]])[1],2), " %")) & ylab(label = paste0(input$Reduction, " / PCA 2 : ", round(Seurat::Stdev(r$dataset[["pca"]])[2],2), " %"))
