
source(file = "/home/boris/Bureau/scShiny/script/work_tools.R")
load(file = '/home/boris/Bureau/scShiny/script/DE/signature.Rdata')
library(harmony)


seurat <- get(load(filsource(file = "/home/boris/Bureau/scShiny/script/work_tools.R")e = "/home/boris/Documents/analyse/All_Post_greffe_B.Rdata"))
Idents(singlet) <- "Condition" 
singlet <- subset(singlet, idents = c('RCHOP','Excipient'))

singlet <- RunPCA(singlet, verbose = F, feature = signature_RE_bulk, approx=F) 
singlet <- FindNeighbors(singlet, dims = 1:5, verbose = F) 
singlet <- FindClusters(singlet, verbose = F) 
singlet <- RunTSNE(singlet, verbose = F, features = signature_RE_bulk )

singlet <- RunUMAP(singlet, verbose = F, features = rownames(singlet[['SCT']])[which(rownames(singlet[['SCT']]) %in% signature_RE_bulk[1:40])]) #uwot-learn
Seurat::DimPlot(object = singlet, label.size = 5, pt.size = 0.5, reduction = 'umap', label = TRUE) & theme(legend.position = "top") & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4)))

save(singlet, file = "/home/boris/Documents/analyse/FL08G0431/RE_FL08G0431.RData")

singlet <- RunHarmony(singlet, "Condition", assay.use = "SCT")
singlet <- RunUMAP(singlet, reduction = "harmony", features = signature_RE_bulk)
Idents(singlet) <- "Condition" 
Seurat::DimPlot(object = singlet, label.size = 5, pt.size = 0.5, reduction = 'harmony', label = TRUE) & theme(legend.position = "top") & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4)))


seurat <- get(load(file = "/home/boris/Documents/analyse/FL140304/RE_FL140304.RData"))
SaveH5Seurat(seurat, filename = "RE_FL140304.h5Seurat")
Convert("RE_FL140304.h5Seurat", dest = "/home/boris/Documents/analyse/RE_FL140304.h5ad")






pca1 <- Seurat::DimPlot(object = singlet, label.size = 5, pt.size = 0.5, reduction = 'pca', label = TRUE, dims = c(1,2)) & theme(title = element_text(size=20),legend.position = "top",legend.title = element_text(size=10),legend.text = element_text(size=10)) & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 6))) & xlab(label = paste0("PCA 1 : ", round(Seurat::Stdev(singlet[["pca"]])[1],2), " %")) & ylab(label = paste0("PCA 2 : ", round(Seurat::Stdev(singlet[["pca"]])[2],2), " %"))
umap1 <- Seurat::DimPlot(object = singlet, label.size = 5, pt.size = 0.5, reduction = 'umap', label = TRUE) & theme(legend.position = "top") & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4)))
tsne1 <- Seurat::DimPlot(object = singlet, label.size = 5, pt.size = 0.5, reduction = 'tsne', label = TRUE) & theme(legend.position = "top") & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4)))
cowplot::plot_grid(nrow = 1,ncol = 3, pca1, umap1, tsne1)





a <- list()
for (i in 2:50) {
  print(i)
  singlet <- RunUMAP(singlet, dims = 1:i, verbose = F)
  a[[i]] <- Seurat::DimPlot(object = singlet, label.size = 5, pt.size = 0.5, reduction = 'umap', label = TRUE) & NoAxes() & NoLegend()
}

cowplot::plot_grid(nrow = 5,ncol = 10,
                   a[[2]],a[[3]],a[[4]],a[[5]],a[[6]],a[[7]],a[[8]],a[[9]],
                   a[[10]],a[[11]],a[[12]],a[[13]],a[[14]],a[[15]],a[[16]],a[[17]],a[[18]],a[[19]],
                   a[[20]],a[[21]],a[[22]],a[[23]],a[[24]],a[[25]],a[[26]],a[[27]],a[[28]],a[[29]],
                   a[[30]],a[[31]],a[[32]],a[[33]],a[[34]],a[[35]],a[[36]],a[[37]],a[[38]],a[[39]],
                   a[[40]],a[[41]],a[[42]],a[[43]],a[[44]],a[[45]],a[[46]],a[[47]],a[[48]],a[[49]], a[[50]])


b <- list()
for (i in 2:50) {
  print(i)
  singlet <- RunUMAP(singlet, dims = 1:i, verbose = F, n.neighbors =  i, n.epochs = 500) 
  b[[i]] <- Seurat::DimPlot(object = singlet, label.size = 5, pt.size = 0.5, reduction = 'umap', label = TRUE) & NoAxes() & NoLegend()
}


cowplot::plot_grid(nrow = 5,ncol = 10,
                   b[[2]],b[[3]],b[[4]],b[[5]],b[[6]],b[[7]],b[[8]],b[[9]],
                   b[[10]],b[[11]],b[[12]],b[[13]],b[[14]],b[[15]],b[[16]],b[[17]],b[[18]],b[[19]],
                   b[[20]],b[[21]],b[[22]],b[[23]],b[[24]],b[[25]],b[[26]],b[[27]],b[[28]],b[[29]],
                   b[[30]],b[[31]],b[[32]],b[[33]],b[[34]],b[[35]],b[[36]],b[[37]],b[[38]],b[[39]],
                   b[[40]],b[[41]],b[[42]],b[[43]],b[[44]],b[[45]],b[[46]],b[[47]],b[[48]],b[[49]], b[[50]])

