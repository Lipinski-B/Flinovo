source(file = "/home/boris/Bureau/scShiny/document/script/work_tools.R")
patient <- siege[6]
load(file = paste0("/home/boris/Documents/analyse/", patient,"/singlet_", patient,".RData"))

################################################################################################################################################################################################################################################################################################################################################
## -- TEST : DE RCHOP : apop+ / apop- -- ##
patient <- "FL140304"
load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
Idents(singlet) <- "Phénotype" ; sub_singlet <- subset(singlet, idents = "B-cells")
Idents(singlet) <- "Condition" ; singlet <- subset(singlet, idents = c("RCHOP", "Excipient"))
seuil = 0.13
liste = as.matrix(rownames(sub_singlet@meta.data)[which(sub_singlet@meta.data$APOPTOSIS<seuil & sub_singlet@meta.data$Phénotype=="B-cells" & sub_singlet@meta.data$Condition =="RCHOP")])
liste2 = as.matrix(rownames(sub_singlet@meta.data)[which(sub_singlet@meta.data$APOPTOSIS>=seuil & sub_singlet@meta.data$Phénotype=="B-cells" & sub_singlet@meta.data$Condition =="RCHOP")])
sub_singlet@meta.data[sub_singlet@meta.data$Condition == "RCHOP" & colnames(sub_singlet) %in% liste, "Apop" ] <- "+"
sub_singlet@meta.data[sub_singlet@meta.data$Condition == "RCHOP" & colnames(sub_singlet) %in% liste2, "Apop" ] <- "-"


singlet <- visualisation(singlet)
Idents(singlet) <- "Condition" ; DimPlot(object = singlet, dims = c(1, 2), label.size = 7, pt.size = 1.3, reduction = 'pca', label = TRUE) & NoLegend() & theme(title = element_text(size=20), legend.position = "top",legend.title = element_text(size=10),legend.text = element_text(size=10)) & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 6))) & xlab(label = paste0("PCA 1 : ", round(Seurat::Stdev(singlet[["pca"]])[1],2), " %")) & ylab(label = paste0("PCA 2 : ", round(Seurat::Stdev(singlet[["pca"]])[2],2), " %"))
