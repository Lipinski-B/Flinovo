source(file = "/home/boris/Bureau/Flinovo/script/work_tools.R")
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
set.seed(1234)


BD <- celldex::MonacoImmuneData()
for (patient in siege){
  for (condition in c("singlet","Post_greffe","Pré_greffe")){
    load(file = paste0("/home/boris/Documents/analyse/",patient,"/",condition,"_",patient,".RData"))
    results.blueprint <- SingleR::SingleR(test = singlet@assays$RNA@data, ref = BD, labels = BD$label.main) ; labels <- results.blueprint$pruned.labels
    results.blueprint.fine <- SingleR::SingleR(test = singlet@assays$RNA@data, ref = BD, labels = BD$label.fine) ; labels.fine <- results.blueprint.fine$pruned.labels
    singlet <- AddMetaData(singlet, labels, col.name = "Monaco")
    singlet <- AddMetaData(singlet, labels.fine , col.name = "Monaco.fine")
    percent.rb <- singlet[["percent.rb"]] <- PercentageFeatureSet(singlet, pattern = "^RP[SL]")
    percent.ig <- singlet[["percent.ig"]] <- PercentageFeatureSet(singlet, pattern = "^IG")
    save(singlet, file = paste0("/home/boris/Documents/analyse/",patient,"/",condition,"_",patient,".RData"))
    
    for (sub_population in c("All","BT","B")){
      if (condition == "singlet"){sub_condition = "All"}
      if (condition == "Post_greffe"){sub_condition = "Post-greffe"}
      if (condition == "Pré_greffe"){sub_condition = "Pré-greffe"}
      
      load(file = paste0("/home/boris/Bureau/scShiny/datasets/Patient/",sub_condition,"_",sub_population,"_",patient,".RData"))
      one <- results.blueprint[intersect(colnames(singlet), rownames(results.blueprint)),]
      one <- one$pruned.labels
      
      two <- results.blueprint.fine[intersect(colnames(singlet), rownames(results.blueprint.fine)),]
      two <- two$pruned.labels
      
      tree <- percent.rb[which(colnames(singlet) %in% rownames(percent.rb)),]
      
      four <- percent.ig[which(colnames(singlet) %in% rownames(percent.ig)),]
      
      
      singlet <- AddMetaData(singlet, one, col.name = "Monaco")
      singlet <- AddMetaData(singlet, two , col.name = "Monaco.fine")
      singlet <- AddMetaData(singlet, tree , col.name = "percent.rb")
      singlet <- AddMetaData(singlet, four , col.name = "percent.ig")
      
      save(singlet, file = paste0("/home/boris/Bureau/scShiny/datasets/Patient/",sub_condition,"_",sub_population,"_",patient,".RData"))
    }
  }
}


singlet <- seurat_subset(singlet, "Monaco.fine", names(table(results.blueprint.fine$pruned.labels)[which(table(results.blueprint.fine$pruned.labels) > 10)]))
singlet <- SCTransform(singlet, method = "glmGamPoi", vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F, ncells = 10000, conserve.memory = T)
singlet <- RunPCA(singlet, verbose = F, approx=F) %>% RunUMAP(dims = 1:20, verbose = F) %>% FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(resolution = 0.4, verbose = F, method = "igraph", algorithm=1) %>% RunTSNE(dims = 1:20, verbose = F, perplexity = 100 )




Seurat::DimPlot(object = singlet, group.by = "Chain", label.size = 5, pt.size = 1, reduction = "umap", label = T)
Seurat::DimPlot(object = singlet, group.by = c("V","D","J","C"), label.size = 5, pt.size = 1, reduction = "umap", label = T)
Seurat::DimPlot(object = singlet, group.by = "Monaco.fine", label.size = 5, pt.size = 1, reduction = "umap", label = T)
table(singlet@meta.data$Monaco.fine)


# singlet@assays$RNA@data <- singlet@assays$RNA@data[rownames(singlet@assays$RNA@data)[which(!str_detect(rownames(singlet@assays$RNA@data), "^IG"))],]
# singlet@assays$RNA@scale.data <- singlet@assays$RNA@scale.data[rownames(singlet@assays$RNA@scale.data)[which(!str_detect(rownames(singlet@assays$RNA@scale.data), "^IG"))],]
# singlet@assays$RNA@counts <- singlet@assays$RNA@counts[rownames(singlet@assays$RNA@counts)[which(!str_detect(rownames(singlet@assays$RNA@counts), "^IG"))],]
Seurat::FeaturePlot(singlet, features = "IGHD", reduction = "umap", pt.size = 1, combine = T ) & Seurat::NoAxes()

pbmc.markers <- FindAllMarkers(singlet, min.pct = 0.25, test.use = "negbinom", logfc.threshold = 0.25, assay = "RNA", slot = "data")
pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(singlet, features = top10$gene) + NoLegend()




erythroid.cds <- as.cell_data_set(singlet)
erythroid.cds <- cluster_cells(cds = erythroid.cds, reduction_method = "UMAP")
erythroid.cds <- learn_graph(erythroid.cds, use_partition = TRUE)
erythroid.cds <- order_cells(erythroid.cds, reduction_method = "UMAP", root_cells = colnames(singlet)[which(singlet@meta.data$Monaco.fine=="Progenitor cells")])

plot_cells(cds = erythroid.cds,color_cells_by = "pseudotime",cell_size = 1,show_trajectory_graph = TRUE, label_groups_by_cluster=FALSE,reduction_method = "UMAP",label_leaves=TRUE,label_branch_points=TRUE,graph_label_size=5, group_label_size = 4)

singlet <- AddMetaData(object = singlet,metadata = erythroid.cds@principal_graph_aux@listData$UMAP$pseudotime,col.name = "Erythroid")
FeaturePlot(singlet, c("Erythroid"), pt.size = 0.1) & scale_color_viridis_c()



immdata <- repLoad("/home/boris/Documents/lipinskib/Boris_Manon/flinovo/result/FL140304/VDJ/CellrangerVDJ")
