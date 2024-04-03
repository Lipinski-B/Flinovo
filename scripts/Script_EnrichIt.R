# Load Seurat Object RCHOP B cells
load("/Users/Manon/Documents/Manon/Flinovo/Patients/Object_Seurat/All_RCHOP_B.RData")

# Load and score the RCHOP signature (n= 21 genes)
DefaultAssay(singlet) <- "RNA"
cd_features <- list(c("RPS19","BAX", "PSMB4","PHPT1","RPS27L","SRSF3","HIST1H2BK","TRIM22","FDXR","MRFAP1","LY86","PVT1","DDB2","AEN","PQBP1","CD70","BBC3","HNRNPH1","SRSF2","CDKN1A","MT2A"))
singlet <- AddModuleScore(object = singlet, features = cd_features, assay = "RNA", name = 'RCHOP_signature')

#  Add a new metadata to Seurat object
cluster.8 <- subset(singlet, idents = "8")
cellNames_cl8 <- rownames(cluster.8@meta.data)
singlet$cl8 <- rownames(singlet@meta.data)
singlet@meta.data <- singlet@meta.data %>% mutate(classification = ifelse((singlet$cl8 %in% cellNames_cl8), "cluster 8",  "other cluster"))

# enrichIt function
GS.hallmark <- getGeneSets(species = "Homo sapiens", library = "H")
DefaultAssay(singlet) <- "RNA"
ES <- enrichIt(obj = singlet, gene.sets = GS.hallmark, groups = 1000, cores = 6, min.size = 5)
singlet <- Seurat::AddMetaData(singlet, ES)

ES2 <- data.frame(singlet[[]], Idents(singlet))

output <- getSignificance(ES2, group = "classification", gene.sets = names(ES),fit = "Wilcoxon")
output <- output[order(output$FDR),]
output <- as.matrix(subset(output , FDR < 0.05))
output <- as.matrix(subset(output, output[,5]>output[,4]))
