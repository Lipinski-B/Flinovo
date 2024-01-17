source(file = "/home/boris/Bureau/scShiny/document/script/work_tools.R")
patient <- siege[6]
load(file = paste0("/home/boris/Documents/analyse/", patient,"/singlet_", patient,".RData"))

## -- TEST : DE appareillé : RCHOP/EXP + PG : RC/RP -- ##
## test 1 = Réponse : RP/RC -- ##
load(file = paste0("/home/boris/Documents/analyse/singlet_all.RData"))
all <- seurat_subset(all, "Condition", c("Pré-greffe"))
all@meta.data[all@meta.data$orig.ident == "FL140304", "Reponse"] <- "RP"
all@meta.data[all@meta.data$orig.ident == "FL12C1888", "Reponse"] <- "RP"
all@meta.data[all@meta.data$orig.ident == "FL08G0293", "Reponse"] <- "RP"

all@meta.data[all@meta.data$orig.ident == "FL09C1164", "Reponse"] <- "RC"
all@meta.data[all@meta.data$orig.ident == "FL02G095", "Reponse"] <- "RC"
all@meta.data[all@meta.data$orig.ident == "FL05G0330", "Reponse"] <- "RC"

Idents(all)<-"Reponse" ; all[["RNA"]]@counts <- as.matrix(all[["RNA"]]@counts)+1
all@tools$DE_R <- FindMarkers(all, slot = "counts", ident.1 = "RP", ident.2 = "RC", test.use = "DESeq2", max.cells.per.ident = 2000)
all[["RNA"]]@counts <- as.matrix(all[["RNA"]]@counts)-1 ; Idents(all)<-"seurat_clusters"
save(all, file = "/home/boris/Documents/analyse/singlet_all_PG.RData")

## test 2 = Réponse : RP/RC -- ##


