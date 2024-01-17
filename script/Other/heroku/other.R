#global.R
#FL09 : réponse modéré

#all@meta.data$SingleR.calls.fine<- all@meta.data$d_gene <- all@meta.data[70] <- all@meta.data$v_gene <- all@meta.data$j_gene <- all@meta.data$c_gene <- all@meta.data$old.ident <- all@meta.data$BCL2_L23L <- all@meta.data$BCL2_K22K <- all@meta.data$CD79B_Y696H <- all@meta.data$EZH2_A682G_1 <- all@meta.data$EZH2_A682G_2  <- all@meta.data$EZH2_A692V_1  <- all@meta.data$EZH2_A692V_2  <- all@meta.data$EZH2_Y646S  <- all@meta.data$EZH2_Y646N  <- all@meta.data$EZH2_Y646H  <- all@meta.data$EZH2_Y646F  <- all@meta.data$EZH2_Y646C <- NULL

#annotations <- read.csv("/home/boris/Bureau/scShiny/annotation_FindAllMarkers.csv")
#Idents(singlet) <- "Phénotype"
#singlet@commands[["FindAllMarkers"]] <- FindAllMarkers(singlet, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
#invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)), detach, character.only = TRUE, unload = TRUE))

##### -- identification des cellules -- ##### 
# DB: DatabaseImmuneCellExpressionData()
#hpca.se <- celldex::DatabaseImmuneCellExpressionData()
#results <- SingleR(test = as.SingleCellExperiment(singlet), ref = hpca.se, labels = hpca.se$label.main)
#results.fine <- SingleR(test = as.SingleCellExperiment(singlet), ref = hpca.se, labels = hpca.se$label.fine)

#singlet$SingleR.pruned.calls <- results$pruned.labels
#singlet$SingleR.calls <- results$labels
#singlet$SingleR.pruned.calls.fine <- results.fine$pruned.labels
#singlet$SingleR.calls.fine <- results.fine$labels


# DB: BlueprintEncodeData()
#singlet$SingleR.pruned.calls.blueprint <- results$pruned.labels
#singlet$SingleR.pruned.calls.blueprint.fine <- results.blueprint.fine$pruned.labels


#singlet <- RunSPCA(singlet, features = VariableFeatures(singlet), ndims.print = 1:10, nfeatures.print = 30, graph = singlet@graphs[["RNA_snn"]])# Reduction dimension
## -- Ohter PCA analysis -- ## 
#singlet <- JackStraw(singlet, num.replicate = 100, reduction = "pca")
#singlet <- ScoreJackStraw(singlet, dims = 1:20)




## -- FindAllMarkers -- ## 
singlet@commands[["FindAllMarkers"]] <- FindAllMarkers(singlet, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
#singlet@commands[["FindAllMarkers"]] <- merge(singlet@commands[["FindAllMarkers"]], annotations, by.x="gene", by.y="gene_name")
#singlet@commands[["FindAllMarkers"]] %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)


#seurat.R
save(all, file = paste0(patient,"/R/singlet_",patient,".RData")) #all[["RNA"]]@counts


## -- Summary metadata -- ## 
result=list() ; resume=list()
summary <- list(all@meta.data)
names(summary) <- patient
for (elm in 1:length(names(summary[[patient]]))) {
  if (names(summary[[patient]])[elm] %in% meta_variable) {
    resume[[names(summary[[1]])[elm]]] = as.list(as.vector(table(summary[[1]][[elm]])))
    names(resume[[names(summary[[1]])[elm]]]) = rownames(as.matrix(table(summary[[1]][[elm]])))
  }
}
result[[patient]] <- resume
save(result, file = "/home/boris/Documents/analyse/shiny/result.RData")
load(file = "/home/boris/Documents/analyse/shiny/result.RData")



singlet <- seurat_object(patient)
singlet <- visualitation(singlet, maximum=2500, percent_mt=30)
singlet <- metadata(singlet)
all<-singlet
save(all, file = paste0("/home/boris/Documents/analyse/singlet_FL12_30.RData"))

######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## To test                                                                                                                       ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
#annotations_finale <- all@meta.data[which(all@meta.data$gene_name %in% top10),] 
#all <- BuildClusterTree(all)
#Tool(all,'BuildClusterTree')
#PlotClusterTree(all)

#all <- SCTransform(all)
#Stdev(all[["pca"]])

sce.all <- as.SingleCellExperiment(all)
colData(sce.all)

sce.all@metadata <- all@meta.data
metadata(sce.all)$which_qc
fluidigm <- as.SingleCellExperiment(all)

#Add QC for Mitochondrial genes
fluidigm <- addPerCellQC(fluidigm)
fluidigm <- addPerFeatureQC(fluidigm)

## -- Monocle -- ##
monocle <- function(singlet){
  importCDS <- function (otherCDS, seurat_scale=F, import_all = FALSE) {
    if (class(otherCDS)[1] == "Seurat") {
      requireNamespace("Seurat")
      if (!seurat_scale) {
        data <- otherCDS@assays$RNA@counts
      } else {
        data <- otherCDS@assays$RNA@scale.data
      }
      if (class(data) == "data.frame") {
        data <- as(as.matrix(data), "sparseMatrix")
      }
      pd <- tryCatch({
        pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
        pd
      }, error = function(e) {
        pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
        pd <- new("AnnotatedDataFrame", data = pData)
        message("This Seurat object doesn't provide any meta data")
        pd
      })
      if (length(setdiff(colnames(data), rownames(pd))) > 0) {
        data <- data[, rownames(pd)]
      }
      fData <- data.frame(gene_short_name = row.names(data), 
                          row.names = row.names(data))
      fd <- new("AnnotatedDataFrame", data = fData)
      #lowerDetectionLimit <- otherCDS@is.expr
      if (all(data == floor(data))) {
        expressionFamily <- negbinomial.size()
        expr <- "negbinomial.size"
      }
      else if (any(data < 0)) {
        expressionFamily <- uninormal()
        expr <- "unimormal"
      }
      else {
        expressionFamily <- tobit()
        expr <- "tobit"
      }
      print(paste0("expressionFamily ",expr))
      # valid_data <- data[, row.names(pd)]
      monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
                                    #lowerDetectionLimit = lowerDetectionLimit,
                                    expressionFamily = expressionFamily)
      if (import_all) {
        if ("Monocle" %in% names(otherCDS@misc)) {
          otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
          otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
          monocle_cds <- otherCDS@misc$Monocle
          mist_list <- otherCDS
        }
        else {
          mist_list <- otherCDS
        }
      }
      else {
        mist_list <- list()
      }
      if ("var.genes" %in% slotNames(otherCDS)) {
        var.genes <- setOrderingFilter(monocle_cds, otherCDS@var.genes)
      }
      monocle_cds@auxClusteringData$seurat <- mist_list
    }
    else if (class(otherCDS)[1] == "SCESet") {
      requireNamespace("scater")
      message("Converting the exprs data in log scale back to original scale ...")
      data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset
      fd <- otherCDS@featureData
      pd <- otherCDS@phenoData
      experimentData = otherCDS@experimentData
      if ("is.expr" %in% slotNames(otherCDS)) 
        lowerDetectionLimit <- otherCDS@is.expr
      else lowerDetectionLimit <- 1
      if (all(data == floor(data))) {
        expressionFamily <- negbinomial.size()
      }
      else if (any(data < 0)) {
        expressionFamily <- uninormal()
      }
      else {
        expressionFamily <- tobit()
      }
      if (import_all) {
        mist_list <- otherCDS
      }
      else {
        mist_list <- list()
      }
      monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
                                    lowerDetectionLimit = lowerDetectionLimit, expressionFamily = expressionFamily)
      monocle_cds@auxOrderingData$scran <- mist_list
    }
    else {
      stop("the object type you want to export to is not supported yet")
    }
    return(monocle_cds)
  } # Importer object seurat (sous forme de matrix) dans monocle: CellDataSet (CDS)
  ## 1 -- Store Data in a CellDataSet Object & Filtering low-quality cells
  data <- importCDS(singlet, import_all = TRUE)
  data <- estimateSizeFactors(data)         # Estimate size factors 
  data <- estimateDispersions(data)         # and dispersions 
  
  #valid_cells <- row.names(subset(pData(data), Cells.in.Well == 1, Control == FALSE & Clump == FALSE & Debris == FALSE & Mapped.Fragments > 1000000))
  #data <- data[,valid_cells]
  pData(data)$Total_mRNAs <- Matrix::colSums(exprs(data))
  data <- data[,pData(data)$Total_mRNAs < 1e6]
  upper_bound <- 10^(mean(log10(pData(data)$Total_mRNAs)) + 2*sd(log10(pData(data)$Total_mRNAs)))
  lower_bound <- 10^(mean(log10(pData(data)$Total_mRNAs)) - 2*sd(log10(pData(data)$Total_mRNAs)))
  #qplot(Total_mRNAs, data = pData(data),  geom = "density") + geom_vline(xintercept = lower_bound) + geom_vline(xintercept = upper_bound)
  
  data <- data[,pData(data)$Total_mRNAs > lower_bound & pData(data)$Total_mRNAs < upper_bound]
  data <- detectGenes(data, min_expr = 0.1)
  
  ## 2 -- Classify and Counting Cells cells with known marker genes (Done by seurat)
  ## 3 -- Cluster your cells (Done by seurat)
  
  ## 4 -- Workflow to order cells in pseudotime along Single Cell Trajectories
  # a - choosing genes that define progress
  expressed_genes <- row.names(subset(fData(data), num_cells_expressed >= 10))
  diff_test_res <- differentialGeneTest(data[expressed_genes,], fullModelFormulaStr = "~HTO_maxID") #find all genes that are differentially expressed in response to the switch from the model
  ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
  data <- setOrderingFilter(data, ordering_genes) # methods to select genes that require no knowledge of the design of the experiment at all
  
  # b - reducing the dimensionality of the data 
  data <- reduceDimension(data, max_components = 2, method = 'DDRTree')
  
  # c - ordering the cells in pseudotime 
  data <- orderCells(data)
  
  return(data)
}
data <- monocle(all)




######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## ANALYSES                                                                                                                      ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
###############################################################################################################################################
## -- Main visualisation -- ##
##############################
#  TSNE single and multi label
DimPlot(singlet, reduction = "tsne", label = TRUE)
BD <- "SingleR.calls" #"SingleR.calls.fine", "SingleR.calls.blueprint", "SingleR.calls.blueprint.fine"

Idents(all)<-"HTO_maxID"
tokeep <- levels(Idents(all))
tokeep <- tokeep[tokeep %in% "pregreffe"]
all2 <- subset(all, idents = tokeep)
Idents(all2)<-"seurat_clusters"

plots <- TSNEPlot(object = all2, group.by = c("HTO_maxID", BD), split.by = "Phase", label.size = 0.0, pt.size = 1)
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 1)))

# Heatmap de correlation
top10 <- all@commands[["FindAllMarkers"]] %>% group_by("HTO_maxID") %>% top_n(n = 100, wt = avg_log2FC)
DoHeatmap(all, features = top10$gene, group.by = "HTO_maxID") + NoLegend()

# PCA
Idents(all2)<-"seurat_clusters"
DimPlot(singlet, reduction = "pca")
print(all[["pca"]], dims = 1:10, nfeatures = 10)
ElbowPlot(all, ndims = 25, reduction = "pca")
DimHeatmap(all, dims = 1:10, cells = 100, balanced = TRUE)
VizDimLoadings(all, dims = 1:5, reduction = "pca")
JackStrawPlot(all, dims = 1:15)

# UMAP
DimPlot(all, reduction = "umap", label = TRUE)
DimPlot(all, reduction = "umap", group.by = "HTO_maxID")

# TSNE
DimPlot(all, reduction = "tsne", group.by = "HTO_maxID")
DimPlot(all, reduction = "tsne", group.by = BD)
DimPlot(all, reduction = "tsne", group.by = "Phase")
TSNEPlot(object = all, group.by = "clonotype_id", label.size = 0.0, pt.size = 1)

# Mitochondrie 
VlnPlot(all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "HTO_classification")
VlnPlot(all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident") # Visualize QC metrics as a violin plot
FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt") + FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "HTO_classification")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "HTO_classification")
plot1 + plot2

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(all), 10)
VariableFeaturePlot(all)
LabelPoints(plot = VariableFeaturePlot(all), points = top10, repel = TRUE) # plot variable features with and without labels

annotations_finale <- annotations[which(annotations$gene_name %in% top10),] 

# Visualize the distribution of cell cycle markers across
RidgePlot(all, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

# Verification Expression marqueurs B_T
feature <-c("IGHM", "IGHG1")
FeaturePlot(all, features = feature, cols = RColorBrewer::brewer.pal(6,"YlOrBr"), combine = TRUE, blend = TRUE) & NoAxes() & NoLegend()

all1 <- subset(all, subset=IGHM>1)

FeaturePlot(all1, features = feature , reduction='pca')#, split.by = "Phase")
FeaturePlot(all, features = c("MS4A1", "CD19", "CD79A", "CD79B"), reduction='umap')
FeaturePlot(all, features = c("MS4A1", "CD19", "CD79A", "CD79B"), reduction='tsne')
VlnPlot(all, features = c("MS4A1", "CD19", "CD79A", "CD79B"), group.by = "HTO_classification")
FeaturePlot(all, features = c("CD3E", "CD3D", "CD56", "CD4", "CD2", "CD25", "CD62L", "CD197")) # Expression marqueurs T cells
VlnPlot(all, features = c("CD3E", "CD3D", "CD56", "CD4", "CD2", "CD25", "CD62L", "CD197"), group.by = "HTO_classification")



###############################################################################################################################################
## -- Monocle -- ##
###################
# Trajectories visualisations
plot_ordering_genes(data)
plot_cell_trajectory(data, color_by = "HTO_maxID")
plot_cell_trajectory(data, color_by = "SingleR.calls")
plot_cell_trajectory(data, color_by = "State") + facet_wrap(~State, nrow = 1)
plot_cell_trajectory(data, color_by = "Pseudotime")
plot_cell_trajectory(data, color_by = "Cluster")

#jitter plot to pick figure out which state corresponds to rapid proliferation
blast_genes <- row.names(subset(fData(data), gene_short_name %in% c("CD19", "MS4A1")))
plot_genes_jitter(data[blast_genes,], grouping = "State", min_expr = 0.1)

#confirm that the ordering is correct 
data_expressed_genes <-  row.names(subset(fData(data), num_cells_expressed >= 10))
data_filtered <- data[data_expressed_genes,]
my_genes <- row.names(subset(fData(data_filtered), gene_short_name %in% c("CD3E","MS4A1", "CD19")))
cds_subset <- data_filtered[my_genes,]
p <- plot_genes_in_pseudotime(cds_subset, color_by = "SingleR.calls")


###############################################################################################################################################
## -- 3D Visualisation -- ##
############################
# Main visualitation
# cell plotting : Embeddings(object = all, reduction = "umap") --> Visualize what headings are called so that you can extract them to form a dataframe
all2 <- all
all2 <- RunUMAP(all2, reduction = "pca", dims = 1:40, n.components = 3L)
all2 <- RunTSNE(all2, reduction = "pca", dims = 1:40, dim.embed = 3)
meta_variable = c("orig.ident", "HTO_maxID", "SingleR.calls", "clonotype_id","chain","v_gene", "d_gene", "j_gene","c_gene", "cdr3", "Phase", "old.ident", "seurat_clusters")
plot.data <- FetchData(object = all2, vars = c(meta_variable, "PC_1", "PC_2", "PC_3", "tSNE_1", "tSNE_2", "tSNE_3", "UMAP_1", "UMAP_2", "UMAP_3"))
plot.data$label <- paste(rownames(plot.data))

plot_ly(data = plot.data, x = ~PC_1, y = ~PC_2, z = ~PC_3, # ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, or ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
        color = ~seurat_clusters, colors = RColorBrewer::brewer.pal(length(levels(all@meta.data$seurat_clusters)),"Spectral"),
        type = "scatter3d", mode = "markers", 
        marker = list(size = 3, width=2),
        text=~label, hoverinfo="text")

# Feature plot
gene_feature <- "CD19"
all_feature <- rownames(as.matrix(all2[["RNA"]]@counts))
plot.data <- FetchData(object = all2, vars = c(all_feature, "PC_1", "PC_2", "PC_3", "tSNE_1", "tSNE_2", "tSNE_3", "UMAP_1", "UMAP_2", "UMAP_3"), slot = 'data')
plot.data$changed <- ifelse(test = plot.data[[gene_feature]] <1, yes = plot.data[[gene_feature]], no = 1)
plot.data$label <- paste(rownames(plot.data)," - ", plot.data[[gene_feature]], sep="")


plot_ly(data = plot.data, x = ~PC_1, y = ~PC_2, z = ~PC_3, # ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, or ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
        color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        colors = c('darkgreen', 'red'), opacity = .5,
        type = "scatter3d", mode = "markers",
        marker = list(size = 3, width=2), 
        text=~label, hoverinfo="text"
) %>% layout(title=gene_feature)

# Feature plot aleternatif
Cutoff <- quantile(plot.data[,gene_feature], probs = .95)
plot.data$"ExprCutoff" <- ifelse(test = plot.data[,gene_feature] < Cutoff, yes = plot.data[,gene_feature], no = Cutoff)
plot.data$label <- paste(rownames(plot.data)," - ", plot.data[,gene_feature], sep="")

plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, name = gene_feature, # Plot your data, in this example my Seurat object had 21 clusters (0-20), and cells express a gene called ACTB
        color = ~ExprCutoff, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        colors = c('darkgrey', 'red'), opacity = .5,
        type = "scatter3d", mode = "markers", marker = list(size = 1), 
        text=~label,hoverinfo="text"
) %>% layout(title=gene_feature)



###############################################################################################################################################
## -- ESCAPE : Enrichissement de gène -- ##
###########################################
#The Heatmap
all@meta.data$active.idents <- all@active.ident
dittoHeatmap(all, genes = NULL, metas = all@tools$hallmarks, heatmap.colors = rev(colorblind_vector(50)),
             annot.by = "HTO_maxID", cluster_cols = F, fontsize = 7, order.by = "HTO_maxID")

dittoHeatmap(all, genes = NULL, metas = names(ES), heatmap.colors = rev(colorblind_vector(50)),
             annot.by = meta_variable, cluster_cols = TRUE, fontsize = 7)

colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
dittoHeatmap(all, genes = NULL, metas = c("HALLMARK_APOPTOSIS", "HALLMARK_DNA_REPAIR", "HALLMARK_P53_PATHWAY"), 
             heatmap.colors = rev(colorblind_vector(50)), annot.by = meta_variable, cluster_cols = TRUE, fontsize = 7)

#The Violin Plot
dittoPlot(all, "HALLMARK_DNA_REPAIR", group.by = "SingleR.calls") #+ scale_fill_manual(values = colorblind_vector(5))

#Hex Density Enrichment Plots
dittoScatterHex(all,x.var = "HALLMARK_DNA_REPAIR", y.var = "HALLMARK_MTORC1_SIGNALING", do.contour = TRUE) + theme_classic() + 
  scale_fill_gradientn(colors = rev(colorblind_vector(11))) + geom_vline(xintercept = 0, lty=2) + geom_hline(yintercept = 0, lty=2)  

dittoScatterHex(all, x.var = "HALLMARK_DNA_REPAIR", y.var = "HALLMARK_MTORC1_SIGNALING", do.contour = TRUE, split.by = "SingleR.calls") + 
  theme_classic() + scale_fill_gradientn(colors = rev(colorblind_vector(11))) + geom_vline(xintercept = 0, lty=2) + geom_hline(yintercept = 0, lty=2) 

# Enrichment along a Ridge Plot
ES2 <- data.frame(all[[]], Idents(all))
colnames(ES2)[ncol(ES2)] <- "cluster"
ridgeEnrichment(ES2, gene.set = "HALLMARK_DNA_REPAIR", group = "SingleR.calls", add.rug = TRUE)
ridgeEnrichment(ES2, gene.set = "HALLMARK_DNA_REPAIR", group = "cluster", facet = "SingleR.calls", add.rug = TRUE)

# The Split Violin Plot
splitEnrichment(ES2, split = "SingleR.calls", gene.set = "HALLMARK_DNA_REPAIR")
splitEnrichment(ES2, x.axis = "cluster", split = "SingleR.calls", gene.set = "HALLMARK_DNA_REPAIR")

# Expanded Analysis
ES2 <- data.frame(all[[]], Idents(all))
PCA <- performPCA(enriched = ES2, groups = c("cluster", "SingleR.calls"))
pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = TRUE)
pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = FALSE, facet = "cluster") 
masterPCAPlot(ES2, PCx = "PC1", PCy = "PC2", top.contribution = 10)

#Signficance
output <- getSignificance(ES2, group = "cluster", fit = "linear.model")





######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## VDJ - Métadata                                                                                                                ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
## ATTENTION pour que installation scRepertoire ok il faut d'abord désintaller ggplot2 et Seurat, installer scRepertoire puis re installer ggplot2 et Seurat
suppressMessages(library(ggplot2))
suppressMessages(library(Seurat))
library(scRepertoire)

# Chargement des données VDJ
patient1 <- "hFL_180008B"
csv1 <- list(read.csv(paste0(patient1,"/VDJ/",patient1,"_CellrangerVDJ/outs/filtered_contig_annotations.csv"), stringsAsFactors = FALSE))
patient2 <- "hFL_130337"
csv2 <- list(read.csv(paste0(patient2,"/VDJ/",patient2,"_CellrangerVDJ/outs/filtered_contig_annotations.csv"), stringsAsFactors = FALSE))

combined <- combineBCR(c(csv1,csv2), samples = c(patient1, patient2), ID = c("BCR","BCR"))

quantContig(combined, cloneCall="gene+nt", scale = TRUE)
quantContig_output <- quantContig(combined, cloneCall="gene+nt", scale = TRUE, exportTable = TRUE)
quantContig_output




######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## CellSIUS - Cell Subtype Identification from Upregulated gene Sets                                                             ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
all_matrix <- as.matrix(all[["RNA"]]@counts)
CellSIUS.out <- CellSIUS(all_matrix, as.character(Idents(all)))
Result_List = CellSIUS_GetResults(CellSIUS.out=CellSIUS.out)

my_tsne <- Embeddings(all[["tsne"]])
d = CellSIUS_plot(coord = my_tsne, CellSIUS.out = CellSIUS.out)
d + scale_color_brewer(palette = "Blues")

Final_Clusters = CellSIUS_final_cluster_assignment(CellSIUS.out=CellSIUS.out, group_id=Idents(all), min_n_genes = 3)
table(Final_Clusters)
p <- qplot(my_tsne[,1],my_tsne[,2],xlab="tSNE1",ylab="tSNE2",color=as.factor(Final_Clusters))
p + theme(legend.position = "none")




######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## GENOTYPAGE - Long Ranger                                                                                                      ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
#vawk '{print $1,$2}' variants.vcf > SNV.loci.txt
#sed -i 's/\s/:/g' SNV.loci.txt 
setwd(dir = "/home/boris/Documents/lipinskib/")

# Construct the data.frame
barcodes <- read.table("/home/boris/Documents/lipinskib/flinovo/datas/datas_analysees/hFL_180008B/hFL_180008B_CellrangerCount/outs/filtered_feature_bc_matrix/barcodes.tsv", header = F) #read in the cell barcodes output by Cell Ranger
snps <- read.table("/home/boris/Documents/lipinskib/boris/longranger/result/wgs/hFL_180008B/outs/SNV.loci.txt", header = F) # Construct the final table to add to the Seurat object

# charge raw data
snv_matrix <- readMM("/home/boris/Documents/lipinskib/boris/longranger/result/vartex/hFL_180008B/matrix.mtx") # Read in the sparse genotype matrix
colnames(snv_matrix) <- barcodes$V1
row.names(snv_matrix) <- snps$V1 #View(head(snv_matrix[1:10]))

snv_matrix <- as.data.frame(as.matrix(t(snv_matrix))) # convert the matrix to a dataframe #View(head(snv_matrix[1:10]))
#sub_snv_matrix <- snv_matrix[98640:98825] #BCL2 target

colnames(snv_matrix)
taget_name <-"chr18:63319316" #or colnames(snv_matrix)[1]
target <- data.frame(snv_matrix[taget_name])

row.names(target) <- barcodes$V1 ; colnames(target) <- taget_name

target[taget_name] <- str_replace(as.character(target[taget_name]), "0", "No Call") # No reads detected
target[taget_name] <- str_replace(as.character(target[taget_name]), "1", "ref/ref") # Only ref detected
target[taget_name] <- str_replace(as.character(target[taget_name]), "2", "alt/alt") # Only alt detected
target[taget_name] <- str_replace(as.character(target[taget_name]), "3", "alt/ref") # Both alleles detected

target <- t(target)
target <- target[ , which(colnames(target) %in% colnames(all))] 

all <- AddMetaData(object = all, metadata = as.data.frame(target), col.name = colnames(as.data.frame(target)))
TSNEPlot(object = all, group.by = "target", label.size = 0.0, pt.size = 1)
#do.label = T, colors.use = c("azure2","black", "yellow","red"), plot.title = "chr18_63319316", do.return = T, plot.order = c("alt/alt" ,"alt/ref", "ref/ref","No Call"))
#tableau <- data.frame(c(1:638),c(1:638),c(1:638),c(1:638)); for (i in 1:length(snv_matrix)){tableau[i,] <- table(factor(as.vector(snv_matrix[i]), levels = c(0:3)))} ;colnames(tableau) <- c('0','1','2','3')



#library(monocle)
#library("tidyverse") # for data manipulation & plots
#library("lubridate") # right color
#library("BUSpaRse")
#library("TENxBUSData")
#library("simpleSingleCell") # for scRNA storage & manipulation
#library("scater") # for QC control
#library("scran") # analysis pipeline
#library("uwot") # UMAP dim-red
#library("DropletUtils") #utility functions for handling single-cell (RNA-seq)
#library("AnnotationHub") # ensbl query
#library("AnnotationDbi") # ensbl query
#library("sctransform") # sc normalization,
#library(Matrix)
#library(scRNAseq)
#require(RColorBrewer)
#library(ggplot2)
#library(rowr)
#library(grid)
#library(gridExtra)
#library(Signac)
#library(ape)
#library(plotly)
#library(lemon)
#library(cowplot)
#library(clues)
#library(cellrangerRkit)
#library(treemap)
#suppressPackageStartupMessages(library(escape))
#suppressPackageStartupMessages(library(dittoSeq))
#library(sunburstR)
#library(CellSIUS)
















## -- fonction.R -- ## 
## -- Namanm -- ## 
seurat_object_namanm <- function(i){
  ##### -- Préprocessing -- #####
  # Matrice mRNA + object seurat
  rwa.mRNA <- Read10X(data.dir = paste0("/home/boris/Documents/lipinskib/out_nam/",i,"/outs/filtered_feature_bc_matrix/")) 
  cso <- CreateSeuratObject(counts = rwa.mRNA, project = as.character(i), min.cells = 3, min.features = 200)
  umis <- GetAssayData(object = cso, slot = "counts")
  
  # Matrice HTO
  #raw.hto <- Read10X(paste0(patient,"/HTO/umi_count/"), gene.column = 1)        
  #colnames(raw.hto) <- paste0(colnames(raw.hto),"-1")                          
  #hto <- raw.hto[c(1:3),]                                                       # Suppression des séquences unmapped : rownames(raw.hto)
  #rownames(hto) <- c("Pré-greffe","Excipient","RCHOP")    
  
  #joint.bcs <- intersect(colnames(umis),colnames(hto))                          # Sélection des cellules avec barcode commun HTO / mRNA
  
  #umis <- umis[, joint.bcs]                                                     # Sélection des lignes qui correspondent aux cellules en commun
  #hto <- as.matrix(hto[, joint.bcs])
  
  # Object seurat : HTOs + mRNA
  hashtag <- CreateSeuratObject(counts = umis, assay = "RNA", project = as.character(i))
  hashtag <- NormalizeData(hashtag)
  hashtag <- FindVariableFeatures(hashtag, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  hashtag <- ScaleData(hashtag,features = VariableFeatures(hashtag))
  
  #hashtag[["HTO"]] <- CreateAssayObject(counts = hto)                           # Ajoute des données HTO comme un nouvel assay indépendant du mRNA
  #hashtag <- NormalizeData(hashtag, assay = "HTO", normalization.method = "CLR")
  
  # Association : cellules / échantillons
  #hashtag <- HTODemux(hashtag, assay = "HTO", positive.quantile = 0.99)
  #table(hashtag$HTO_maxID)  
  #table(hashtag$HTO_classification.global)                                      # Result
  
  Idents(hashtag)<- 'HTO_classification.global'
  
  singlet <- hashtag
  #colnames(singlet@meta.data)[which(colnames(singlet@meta.data)=="HTO_maxID")] <- "Condition"
  
  singlet[["percent.mt"]] <- PercentageFeatureSet(singlet, pattern = "^MT-")
  singlet@tools$mitochondrie_all <- VlnPlot(singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  return(singlet)
}
metadata_namnam <- function(singlet){  
  singlet <- all
  # DB: BlueprintEncodeData()
  BD <- celldex::BlueprintEncodeData()
  results.blueprint <- SingleR(test = as.SingleCellExperiment(singlet), ref = BD, labels = BD$label.main)
  results.blueprint.fine <- SingleR(test = as.SingleCellExperiment(singlet), ref = BD, labels = BD$label.fine)
  
  #singlet$SingleR.pruned.calls.blueprint <- results$pruned.labels
  singlet$SingleR.calls <- results.blueprint$labels
  #singlet$SingleR.pruned.calls.blueprint.fine <- results.blueprint.fine$pruned.labels
  singlet$SingleR.calls.fine <- results.blueprint.fine$labels
  
  
  ##### -- Enrichissement des gènes -- ##### 
  #GS <- getGeneSets(library = "H")
  #ES <- enrichIt(obj = singlet, gene.sets = GS, groups = 1000, cores = 12)
  #names(ES) <- str_replace_all(names(ES), "HALLMARK_", "")
  #singlet <- AddMetaData(singlet, ES)
  #singlet@tools$hallmarks <- names(ES)
  #singlet@tools$meta_variable <- c("seurat_clusters", "Condition", "Greffe", "Phénotype", "clonotype_id", "Phase", "BCL2_K22K", "BCL2_L23L", "CD79B_Y696H")# c("seurat_clusters", "Condition", "Phénotype", "Phase", "K29Q", "L37M", "M11I", "pGln45")
  
  
  ##### -- Phase cycle -- #### 
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  singlet <- CellCycleScoring(singlet, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)    # Attribution à chaque cellule d'un état de phase
  #singlet <- ScaleData(singlet, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(singlet)) # Correction de la matrice d'expression dépendemment du cycle cellulaire
  Idents(singlet)<-"seurat_clusters"
  
  ## -- Ohter -- ## 
  singlet@meta.data$nFeature_HTO <- NULL
  singlet@meta.data$HTO_secondID <- NULL
  singlet@meta.data$HTO_classification <- NULL
  singlet@meta.data$RNA_snn_res.0.5 <- NULL
  singlet@meta.data$hash.ID <- NULL
  colnames(singlet@meta.data)[which(colnames(singlet@meta.data)=="SingleR.calls")] <- "Phénotype"
  
  all<- singlet
  
  ##### -- Sunburst -- #### 
  make_sunburst_data <- function(met){
    df <- data.frame(ids=character(), labels=character(), parents=character(), values=integer(),stringsAsFactors=FALSE)
    
    df1 <- met %>% mutate(path = paste(patient,  sep=";")) %>% dplyr::select(path, value)
    df3 <- met %>% mutate(path = paste(patient, phénotype, sep=";")) %>% dplyr::select(path, value)
    df4 <- met %>% mutate(path = paste(patient, phénotype, sub, sep=";")) %>% dplyr::select(path, value)
    df5 <- met %>% mutate(path = paste(patient, phénotype, sub, phase, sep=";")) %>% dplyr::select(path, value)
    
    for (dfX in c(df1,df3,df4,df5)) {
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
    df$values <- c(table(df1)[,1],table(df3)[,1],table(df4)[,1],table(df5)[,1])
    return(df)
  }
  met <- data.frame(
    patient = singlet@meta.data$orig.ident, 
    phénotype = singlet@meta.data$Phénotype, 
    sub = singlet@meta.data$SingleR.calls.fine, 
    phase = singlet@meta.data$Phase,
    value = rep(1, length(singlet@meta.data$orig.ident)), stringsAsFactors = FALSE
  )
  singlet@tools$sunburst <- make_sunburst_data(met)
  
  
  return(singlet)
}


## -- Other -- ##
#.libPaths( c( .libPaths(), '/home/boris/R/x86_64-pc-linux-gnu-library/4.1/') )
#library(dplyr)
#library(SingleR)
#library(stringr)
#library(celldex)
#library(monocle)
#library(escape)

#library(dashboardthemes)
#library(shinyjs)
#library(ggplot2)
#library(shinybusy)
#library(Seurat)
#library(dittoSeq)

visualisation2 <- function(singlet){
  ## -- Pre-processing  -- ##
  singlet <- NormalizeData(singlet, normalization.method = "LogNormalize", scale.factor = 10000) 
  singlet <- FindVariableFeatures(singlet, selection.method = "vst", nfeatures = 2000)
  singlet <- ScaleData(singlet, features = rownames(singlet))                                               # Scale data :singlet@assays$RNA@scale.data, singlet[["RNA"]]@scale.data
  singlet <- RunPCA(singlet, features = VariableFeatures(singlet), ndims.print = 1:10, nfeatures.print = 30)# Reduction dimension
  singlet <- FindNeighbors(singlet, reduction = "pca", dims = 1:40, compute.SNN = T)
  singlet <- FindClusters(singlet, resolution = 0.5)                                                        # head(Idents(singlet), 10) 
  singlet <- RunUMAP(singlet, reduction = "pca", dims = 1:40)
  singlet <- RunTSNE(singlet, reduction = "pca", dims = 1:40)
  
  ## -- FindAllMarkers -- ## 
  #singlet@commands[["FindAllMarkers"]] <- FindAllMarkers(singlet, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
  
  return(singlet)
}