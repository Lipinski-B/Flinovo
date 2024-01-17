source(file = "/home/boris/Bureau/Flinovo/script/work_tools.R")
load("/home/boris/Documents/analyse/FL08G0431/singlet_FL08G0431.RData")
singlet <- seurat_subset(singlet,"Condition","Pré-greffe")
singlet <- seurat_subset(singlet,"Phénotype","B-cells")
singlet <- seurat_subset(singlet,"Phénotype.fine",c("Class-switched memory B-cells","Memory B-cells","naive B-cells","Plasma cells"))
singlet <- subset(singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 2500  )

singlet@assays$RNA@data <- singlet@assays$RNA@data[rownames(singlet@assays$RNA@data)[which(!str_detect(rownames(singlet@assays$RNA@data), "^IG"))],]
singlet@assays$RNA@scale.data <- singlet@assays$RNA@scale.data[rownames(singlet@assays$RNA@scale.data)[which(!str_detect(rownames(singlet@assays$RNA@scale.data), "^IG"))],]
singlet@assays$RNA@counts <- singlet@assays$RNA@counts[rownames(singlet@assays$RNA@counts)[which(!str_detect(rownames(singlet@assays$RNA@counts), "^IG"))],]

#singlet <- normalisation(singlet) %>% visualisation()
singlet <- NormalizeData(singlet)
singlet <- FindVariableFeatures(singlet, selection.method = "vst", nfeatures = 2000)
singlet <- ScaleData(singlet, features = rownames(singlet))
singlet <- visualisation(singlet)


Seurat::DimPlot(singlet, group.by =c("Sample","Condition","Phénotype","Phénotype.fine","Phase","Greffe","seurat_clusters","SCT_snn_res.0.6"), label.size = 5, pt.size = 1, reduction = "umap", label = T) & Seurat::NoLegend()
Seurat::FeaturePlot(singlet, features = c("S.Score","G2M.Score","percent.mt","nCount_RNA","nFeature_RNA","nCount_SCT","nFeature_SCT","nCount_HTO","nFeature_HTO"), reduction = "umap", pt.size = 1, combine = T ) & Seurat::NoAxes()