source(file = "/home/boris/Bureau/scShiny/script/work_tools.R")
patient <- siege[7]
load(file = paste0("/home/boris/Documents/analyse/", patient,"/singlet_", patient,".RData"))


################################################################################################################################################################################################################################################################################################################################################
## -- Workflow -- ## cd10+ bcl2
singlet <- processing(patient)
#save(singlet, file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData")); singlet <- diet(singlet)
#save(singlet, file = paste0("/home/boris/Documents/lipinskib/flinovo/result/",patient,"/R/singlet_", patient,".RData"))
rm(singlet) ; gc() ; gc() ; gc()

for (patient in siege) {
  processing(patient)
  gc() ; gc() ; gc()
}


for (patient in siege) {
  load(file = paste0("/home/boris/Documents/analyse/", patient,"/singlet_", patient,".RData"))
  #...
  save(singlet, file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
  singlet <- diet(singlet)
  save(singlet, file = paste0("/home/boris/Bureau/scShiny/datasets/", patient,".RData"))
}



################################################################################################################################################################################################################################################################################################################################################
## -- to test -- ## cd10+ bcl2
source(file = "/home/boris/Bureau/scShiny/script/work_tools.R")
load(file = "/home/boris/Documents/analyse/All_Post_greffe_B.Rdata")
load(file = '/home/boris/Bureau/scShiny/script/DE/signature.Rdata')

library(CoGAPS)
library(rliger)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
backup <- singlet
Idents(singlet) <- "Condition"

singlet <- RunICA(singlet, assay = 'SCT', ica.function = "icaimax")
Seurat::DimPlot(object = singlet, dims = c(1,2), label.size = 5, pt.size = 0.5, reduction = 'ica', label = TRUE) & theme(legend.position = "top") & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4)))


singlet <- RunLDA(singlet, assay = 'SCT', labels = 'Condition')
Seurat::DimPlot(object = singlet, dims = c(1,2), label.size = 5, pt.size = 0.5, reduction = 'lda', label = TRUE) & theme(legend.position = "top") & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4)))


singlet <- RunSPCA(singlet, assay = 'SCT', graph = singlet@graphs$SCT_snn)
Seurat::DimPlot(object = singlet, dims = c(1,2), label.size = 5, pt.size = 0.5, reduction = 'pca', label = TRUE) & theme(legend.position = "top") & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4)))
Seurat::DimPlot(object = singlet, dims = c(1,2), label.size = 5, pt.size = 0.5, reduction = 'umap', label = TRUE) & theme(legend.position = "top") & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4)))
Seurat::DimPlot(object = singlet, dims = c(1,2), label.size = 5, pt.size = 0.5, reduction = 'tsne', label = TRUE) & theme(legend.position = "top") & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4)))


singlet <- RunSPCA(singlet, assay = 'SCT', features = signature_RE_bulk, graph = singlet@graphs$SCT_nn)
Seurat::DimPlot(object = singlet, dims = c(1,2), label.size = 5, pt.size = 0.5, reduction = 'spca', label = TRUE) & theme(legend.position = "top") & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4)))



load(file = "~/Documents/analyse/All_Post_greffe_B.RData")

singlet <- RunOptimizeALS(singlet, k = 20, lambda = 5, split.by = "Condition")
singlet <- RunQuantileNorm(singlet, split.by = "Condition")
# You can optionally perform Louvain clustering (`FindNeighbors` and `FindClusters`) after
# `RunQuantileNorm` according to your needs
singlet <- FindNeighbors(singlet, reduction = "iNMF", dims = 1:20)
singlet <- FindClusters(singlet, resolution = 0.4)
# Dimensional reduction and plotting
singlet <- RunUMAP(singlet, dims = 1:ncol(singlet[["iNMF"]]), reduction = "iNMF")
DimPlot(singlet, group.by = c("Condition", "orig.ident"))



singlet <- RunCoGAPS(object = singlet, nPatterns = 3, nIterations = 5000, outputFrequency = 1000,sparseOptimization = TRUE, nThreads = 1, distributed = "genome-wide", singleCell = TRUE, seed = 891)
DimPlot(singlet, reduction = "CoGAPS", pt.size = 0.5, dims = c(3, 2))




singlet <- RunFastMNN(object.list = SplitObject(singlet, split.by = "orig.ident"))
singlet <- RunUMAP(singlet, reduction = "mnn", dims = 1:30)
singlet <- FindNeighbors(singlet, reduction = "mnn", dims = 1:30)
singlet <- FindClusters(singlet)
DimPlot(singlet, group.by = c("orig.ident", "Condition"), ncol = 3)





# Initial processing to select variable features
m <- GetAssayData(singlet, slot = "counts", assay = "RNA")
devs <- scry::devianceFeatureSelection(m)
dev_ranked_genes <- rownames(singlet)[order(devs, decreasing = TRUE)]
topdev <- head(dev_ranked_genes, 2000)

# run GLM-PCA on Seurat object. 
# Uses Poisson model by default
# Note that data in the counts slot is used
# We choose 10 dimensions for computational efficiency

ndims <- 10
singlet <- RunGLMPCA(singlet, features = topdev, L = ndims)
singlet <- FindNeighbors(singlet, reduction = 'glmpca', dims = 1:ndims, verbose = FALSE)
singlet <- FindClusters(singlet, verbose = FALSE)
singlet <- RunUMAP(singlet, reduction = 'glmpca', dims = 1:ndims, verbose = FALSE)
Seurat::DimPlot(object = singlet, dims = c(1,2), label.size = 5, pt.size = 0.5, reduction = 'glmpca', label = TRUE) & theme(legend.position = "top") & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4)))



SCTResults(singlet, assay = "SCT", slot = "feature.attributes")
TopCells(singlet)



library(harmony)
singlet <- RunHarmony(singlet, "Condition", assay.use = "SCT")
singlet <- RunUMAP(singlet, reduction = "harmony", features = signature_RE_bulk)



################################################################################################################################################################################################################################################################################################################################################
## -- Fonction Seurat to test -- ##
Idents(singlet) <- "Condition"
head(AverageExpression(singlet))

library(ape)
Idents(singlet) <- colnames(singlet)
singlet <- BuildClusterTree(singlet)
Tool(singlet, slot= "BuildClusterTree")
PlotClusterTree(object = singlet)

GroupCorrelationPlot(singlet,assay = NULL,feature.group = "Condition",cor = "Condition")

plot <- DimPlot(object = singlet)
HoverLocator(plot)

IFeaturePlot(singlet, rownames(singlet@tools$DE_RE), dims = c(1, 2), reduction = NULL, slot = "data")
GetTissueCoordinates(singlet)
ISpatialDimPlot(singlet, image = NULL, group.by = NULL, alpha = c(0.3, 1))
LinkedDimPlot(singlet)
MixscapeHeatmap(singlet, ident.1 = "RCHOP")
PCASigGenes(singlet, pcs.use = 1:2)
PlotPerturbScore(singlet)
PolyDimPlot(singlet)
RelativeCounts
singlet <- RunMixscape(singlet)
VariableFeaturePlot(object = singlet)