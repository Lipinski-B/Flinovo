############################################################################################################################################# 
## -- Workflow -- ##
####################
# Partie 1 : Construction d'un objet Seurat "standard"
# 1 - Chargement de la matrice mRNA
# 2 - Création d'un object Seurat et préporcessing des données
# 3 - Ajout de métadata
# 3.1 - Identification des cellules 
# 3.2 - Enrichissement des gènes
# 4 - Filtre sur les gènes mitochondriaux + empty + doublets
# 5 - Familiarisation avec l'object Seurat : View(metadata) + View(matrix)

# Partie 2 : Ajout des outils de réduction de dimension
# 1 - Preprocessing
# 2 - Add phase cycle information 
# 3 - Ajout de PCA
# 4 - Ajout de UMAP & TSNE
# 5 - Visualisation
# 6 - Sauvegarde

############################################################################################################################################# 
## -- Application -- ##
#######################
setwd(dir = "/home/boris/Documents/lipinskib/boris/Cellranger/result/")
getwd()

# Library à charger
library("Seurat")
library("SingleR")
library(scone)
suppressPackageStartupMessages(library(escape))
suppressPackageStartupMessages(library(dittoSeq))

############################################################################################################################################# 
## -- Partie 1 : Construction d'un objet Seurat "standard" -- ##
################################################################
# 1 - Chargement de la matrice mRNA
rwa.mRNA <- Read10X(data.dir ="/home/boris/Documents/lipinskib/narimene/FLG1435/donnees_analysees/FLG1435_CellrangerCount/outs/filtered_feature_bc_matrix")
patient <- "Nam1"
cso <- CreateSeuratObject(counts = rwa.mRNA, project = patient, min.cells = 3, min.features = 200)
umis <- GetAssayData(object = cso, slot = "counts")

# 2 - Création d'un object Seurat et préporcessing des données
hashtag <- CreateSeuratObject(counts = umis, assay = "RNA", project = patient) # Création d'un object Seurat

hashtag <- NormalizeData(hashtag)   # Préprocessing standard de l'objet seurat, explications dans les vignettes seurat : normalisation
hashtag <- FindVariableFeatures(hashtag, selection.method = "vst", nfeatures = 2000, verbose = FALSE) #trouver les gènes les plus variables
hashtag <- ScaleData(hashtag,features = VariableFeatures(hashtag))                                    #scale

singlet <- hashtag

# 3 - Ajout de métadata
# 3.1 - Identification des cellules
# 3.1.1 - DB1: DatabaseImmuneCellExpressionData()
hpca.se <- celldex::DatabaseImmuneCellExpressionData()                                                  #Récupere la base de donnée
results <- SingleR(test = as.SingleCellExperiment(singlet), ref = hpca.se, labels = hpca.se$label.main) #croisement des données
#results.fine <- SingleR(test = as.SingleCellExperiment(singlet), ref = hpca.se, labels = hpca.se$label.fine) #croisement des données

singlet$SingleR.calls <- results$labels                                                                 #ajout des données aux métadata
#singlet$SingleR.pruned.calls <- results$pruned.labels                                                  #ajout des données aux métadata
#singlet$SingleR.pruned.calls.fine <- results.fine$pruned.labels                                        #ajout ...
#singlet$SingleR.calls.fine <- results.fine$labels


# 3.1.2 - DB2: BlueprintEncodeData()
#BD <- celldex::BlueprintEncodeData()                                                                           #Récupere la base de donnée
#results.blueprint <- SingleR(test = as.SingleCellExperiment(singlet), ref = BD, labels = BD$label.main)        #croisement des données
#results.blueprint.fine <- SingleR(test = as.SingleCellExperiment(singlet), ref = BD, labels = BD$label.fine)   #croisement des données

#singlet$SingleR.pruned.calls.blueprint <- results$pruned.labels                                                #ajout des données aux métadata
#singlet$SingleR.calls.blueprint <- results$labels                                                              #
#singlet$SingleR.pruned.calls.blueprint.fine <- results.blueprint.fine$pruned.labels                            #ajout ...
#singlet$SingleR.calls.blueprint.fine <- results.blueprint.fine$labels

# 3.2 - Enrichissement des gènes
GS <- getGeneSets(library = "H")                                         #Récupere la base de donnée
ES <- enrichIt(obj = singlet, gene.sets = GS, groups = 1000, cores = 12) #croisement des données
names(ES) <- str_replace_all(names(ES), "HALLMARK_", "")

singlet <- AddMetaData(singlet, ES)                                      #ajout des données aux métadata
singlet@tools$hallmarks <- names(ES)

# 4 - Filtre sur les gènes mitochondriaux + empty + doublets
singlet[["percent.mt"]] <- PercentageFeatureSet(singlet, pattern = "^MT-")
singlet <- subset(singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)            # QC Filter : empty < 500 < singlet < 2500 < doublet + MT-5%

# 5 - Familiarisation avec l'object Seurat : View(metadata) + View(matrix)
singlet                                     #couvercle de la boite seurat
View(singlet)                               #ouverture de la boite seurat
View(singlet@meta.data)                     #ouverture du tiroir metadata
View(as.matrix(singlet[["RNA"]]@counts))    #ouverture du tiroir matrice de comptage brute


############################################################################################################################################# 
## -- Partie 2 : Ajout des outils de réduction de dimension -- ##
#################################################################
# 1 - Preprocessing standard de l'objet seurat, explications dans les vignettes seurat
singlet <- NormalizeData(singlet, normalization.method = "LogNormalize", scale.factor = 10000) 
singlet <- FindVariableFeatures(singlet, selection.method = "vst", nfeatures = 2000)
singlet <- ScaleData(singlet, features = rownames(singlet))                                               # Scale data :singlet@assays$RNA@scale.data, singlet[["RNA"]]@scale.data

# 2 - Add phase cycle information 
s.genes <- cc.genes$s.genes                                                                               # Récupere la base de donnée des S
g2m.genes <- cc.genes$g2m.genes                                                                           # Récupere la base de donnée des G2M
singlet <- CellCycleScoring(singlet, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)    # Attribution à chaque cellule d'un état de phase + ajout des données aux métadata
#singlet <- ScaleData(singlet, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(singlet)) # Correction de la matrice d'expression dépendemment du cycle cellulaire

# 3 - Ajout de PCA
singlet <- RunPCA(singlet, features = VariableFeatures(singlet), ndims.print = 1:10, nfeatures.print = 30) # PCA
singlet <- FindNeighbors(singlet, reduction = "pca", dims = 1:40)                                          # trouver les cellules similaires et non similaires
singlet <- FindClusters(singlet, resolution = 0.5)                                                         # trouver les clusters

# 4 - Ajout de UMAP & TSNE
singlet <- RunUMAP(singlet, reduction = "pca", dims = 1:40)                                                # UMAP
singlet <- RunTSNE(singlet, reduction = "pca", dims = 1:40)                                                # TSNE
singlet

# 5 - Visualisation
#réduction de dimension
DimPlot(singlet, reduction = "pca")
DimPlot(singlet, reduction = "pca", group.by = "SingleR.calls")
DimPlot(singlet, reduction = "umap", group.by = "SingleR.calls")
DimPlot(singlet, reduction = "tsne")
DimPlot(singlet, reduction = "tsne", group.by = "SingleR.calls")
DimPlot(singlet, reduction = "tsne", group.by = "Phase", split.by = "SingleR.calls")

#feature plot selon des gènes d'intérêts
FeaturePlot(singlet, features = c("MS4A1", "CD19", "CD79A", "CD79B"), reduction='umap')

#analyses mitochondriales
VlnPlot(singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident") # Visualize QC metrics as a violin plot

# 6 - Sauvegarde
setwd(dir = "/home/boris/Bureau/")
save(singlet, file = "Nam1.Rdata")
#load(file = "Nam1.Rdata")