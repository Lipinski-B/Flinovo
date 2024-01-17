path <- "/home/boris/Bureau/Flinovo/"
path.meta <- "/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_meta/"
path.patient <- "/media/boris/bae7e14e-21e5-48b8-80d6-f94583367b83/Flinovo/analyse_patient/"
source(file = paste0(path,"fonction/library.R"))
source(file = paste0(path,"fonction/workflow.R"))

setwd(dir = "/home/boris/Documents/lipinskib/Boris_Manon/flinovo/result/")
setwd(dir = "/home/boris/")


# Load patient ------------------------------------------------------------
# Get metadata
P1 <- get(load(file = paste0(path.patient, siege[1],"/All_All_", siege[1],".RData"))); rownames(P1@meta.data) <- paste0("FL14_",rownames(P1@meta.data)) ; P1 <- metadata_merge(P1) ; gc() ; gc() ; gc() 
P2 <- get(load(file = paste0(path.patient, siege[2],"/All_All_", siege[2],".RData"))); rownames(P2@meta.data) <- paste0("FL12_",rownames(P2@meta.data)) ; P2 <- metadata_merge(P2) ; gc() ; gc() ; gc() 
P3 <- get(load(file = paste0(path.patient, siege[3],"/All_All_", siege[3],".RData"))); rownames(P3@meta.data) <- paste0("FL09_",rownames(P3@meta.data)) ; P3 <- metadata_merge(P3) ; gc() ; gc() ; gc() 
P4 <- get(load(file = paste0(path.patient, siege[4],"/All_All_", siege[4],".RData"))); rownames(P4@meta.data) <- paste0("FL08_",rownames(P4@meta.data)) ; P4 <- metadata_merge(P4) ; gc() ; gc() ; gc() 
P5 <- get(load(file = paste0(path.patient, siege[5],"/All_All_", siege[5],".RData"))); rownames(P5@meta.data) <- paste0("FL02_",rownames(P5@meta.data)) ; P5 <- metadata_merge(P5) ; gc() ; gc() ; gc() 
P6 <- get(load(file = paste0(path.patient, siege[6],"/All_All_", siege[6],".RData"))); rownames(P6@meta.data) <- paste0("FL05_",rownames(P6@meta.data)) ; P6 <- metadata_merge(P6) ; gc() ; gc() ; gc() 
P7 <- get(load(file = paste0(path.patient, siege[7],"/All_All_", siege[7],".RData"))); rownames(P7@meta.data) <- paste0("FL082_",rownames(P7@meta.data)) ; P7 <- metadata_merge(P7) ; gc() ; gc() ; gc() 
P8 <- get(load(file = paste0(path.patient, siege[8],"/All_All_", siege[8],".RData"))); rownames(P8@meta.data) <- paste0("FL06_",rownames(P8@meta.data)) ; P8 <- metadata_merge(P8) ; gc() ; gc() ; gc() 
P9 <- get(load(file = paste0(path.patient, siege[9],"/All_All_", siege[9],".RData"))); rownames(P9@meta.data) <- paste0("FL083_",rownames(P9@meta.data)) ; P9 <- metadata_merge(P9) ; gc() ; gc() ; gc() 
# P10 <- get(load(file = paste0(path.patient, siege[9],"/All_All_", siege[10],".RData"))); rownames(P10@meta.data) <- paste0("FL083_",rownames(P10@meta.data)) ; P10 <- metadata_merge(P10) ; gc() ; gc() ; gc() 
# P11 <- get(load(file = paste0(path.patient, siege[9],"/All_All_", siege[11],".RData"))); rownames(P11@meta.data) <- paste0("FL083_",rownames(P11@meta.data)) ; P11 <- metadata_merge(P11) ; gc() ; gc() ; gc() 
metadata <- rbind(P1[[]],P2[[]],P3[[]],P4[[]],P5[[]],P6[[]],P7[[]],P8[[]],P9[[]])#,P10[[]],P11[[]]) ; 
rm(P1,P2,P3,P4,P5,P6,P7,P8,P9)#,P10,P11)
metadata$Condition <- metadata$Sample <- metadata$nCount_HTO <- metadata$nCount_RNA <- metadata$nCount_SCT <- metadata$nFeature_HTO <- metadata$nFeature_RNA <- metadata$nFeature_SCT <- metadata$Phénotype<-  metadata$Phénotype.fine<-NULL
save(metadata, file = paste0(path.meta,"metadata.RData"))

# Get patient
P1 <- seurat_object(siege[1]); rownames(P1@meta.data) <- paste0("FL14_",rownames(P1@meta.data)) ; P1 <- metadata_merge(P1) ; gc() ; gc() ; gc() 
P2 <- seurat_object(siege[2]); rownames(P2@meta.data) <- paste0("FL12_",rownames(P2@meta.data)) ; P2 <- metadata_merge(P2) ; gc() ; gc() ; gc() 
P3 <- seurat_object(siege[3]); rownames(P3@meta.data) <- paste0("FL09_",rownames(P3@meta.data)) ; P3 <- metadata_merge(P3) ; gc() ; gc() ; gc() 
P4 <- seurat_object(siege[4]); rownames(P4@meta.data) <- paste0("FL08_",rownames(P4@meta.data)) ; P4 <- metadata_merge(P4) ; gc() ; gc() ; gc() 
P5 <- seurat_object(siege[5]); rownames(P5@meta.data) <- paste0("FL02_",rownames(P5@meta.data)) ; P5 <- metadata_merge(P5) ; gc() ; gc() ; gc() 
P6 <- seurat_object(siege[6]); rownames(P6@meta.data) <- paste0("FL05_",rownames(P6@meta.data)) ; P6 <- metadata_merge(P6) ; gc() ; gc() ; gc() 
P7 <- seurat_object(siege[7]); rownames(P7@meta.data) <- paste0("FL082_",rownames(P7@meta.data)) ; P7 <- metadata_merge(P7) ; gc() ; gc() ; gc() 
P8 <- seurat_object(siege[8]); rownames(P8@meta.data) <- paste0("FL06_",rownames(P8@meta.data)) ; P8 <- metadata_merge(P8) ; gc() ; gc() ; gc() 
P9 <- seurat_object(siege[9]); rownames(P9@meta.data) <- paste0("FL083_",rownames(P9@meta.data)) ; P9 <- metadata_merge(P9) ; gc() ; gc() ; gc() 
# P10 <- seurat_object(siege[10]); rownames(P10@meta.data) <- paste0("FL06_",rownames(P10@meta.data)) ; P10 <- metadata_merge(P10) ; gc() ; gc() ; gc() 
# P11 <- seurat_object(siege[11]); rownames(P11@meta.data) <- paste0("FL083_",rownames(P11@meta.data)) ; P11 <- metadata_merge(P11) ; gc() ; gc() ; gc() 
singlet <- merge(P1, y = c(P2,P3,P4,P5,P6,P7,P8,P9), add.cell.ids = c("FL14","FL12","FL09","FL08","FL02","FL05", "FL082","FL06","FL083"), project = "FL")
rm(P1,P2,P3,P4,P5,P6,P7,P8,P9)#,P10,P11)


# Phenotype
BD <- celldex::MonacoImmuneData()
results.blueprint <- SingleR::SingleR(test = singlet@assays$RNA@data, ref = BD, labels = BD$label.main)
results.blueprint.fine <- SingleR::SingleR(test = singlet@assays$RNA@data, ref = BD, labels = BD$label.fine)
labels <- results.blueprint$labels
labels.fine <- results.blueprint.fine$labels
singlet <- AddMetaData(singlet, labels, col.name = "Phenotype")
singlet <- AddMetaData(singlet, labels.fine , col.name = "Phenotype.fine")

# Reponse
singlet@meta.data[singlet@meta.data$Sample == "FL140304", "Reponse"] <- "RP"
singlet@meta.data[singlet@meta.data$Sample == "FL12C1888","Reponse"] <- "RP"
singlet@meta.data[singlet@meta.data$Sample == "FL08G0293","Reponse"] <- "RP"
singlet@meta.data[singlet@meta.data$Sample == "FL09C1164","Reponse"] <- "RC"
singlet@meta.data[singlet@meta.data$Sample == "FL02G095", "Reponse"] <- "RC"
singlet@meta.data[singlet@meta.data$Sample == "FL05G0330","Reponse"] <- "RC"
singlet@meta.data[singlet@meta.data$Sample == "FL08G0431","Reponse"] <- "RC"
singlet@meta.data[singlet@meta.data$Sample == "FL06G1206","Reponse"] <- "RP"
singlet@meta.data[singlet@meta.data$Sample == "FL08G0404","Reponse"] <- "RC"
#+2RC

# Save
save(singlet, file = paste0(path.meta,"All_0.RData"))

#Integration on server
# /home/boris/Documents/lipinskib/boris/Flinovo/integration/serveur.r
# /home/boris/Bureau/Flinovo/fonction/integration.r

# Metadata ----------------------------------------------------------------
for(condition in c("All", "Pre-greffe", "RCHOP", "Excipient","All","Pre-greffe-Excipient","Post-greffe")){ #
  for(phénotype in c("Other","B","All","BT")){ #
    # Visualisation -----------------------------------------------------------
    feature <- VariableFeatures(singlet)[which(!str_detect(VariableFeatures(singlet), "IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]"))];
    singlet <- RunPCA(singlet, verbose = F, approx=F, feature = feature) %>% RunUMAP(verbose = T, umap.method = "umap-learn", feature = feature) %>% 
      FindNeighbors(feature = feature, verbose = F) %>% FindClusters(resolution = 0.6, verbose = F) %>% RunTSNE(verbose = F, perplexity = 100, feature = feature) 
      #%>% RunVelocity(deltaT = 1, kCells = 50, fit.quantile = 0.02, na.rm = T, ncores = 12)
    
    # Add metadata ------------------------------------------------------------
    load(file = paste0(path,"result/R/All_",condition,"_",phénotype,".RData"))
    load(file = paste0(path.meta,"metadata.RData"))
    
    ## -- All -- ##
    colnames(singlet@meta.data)[which(colnames(singlet@meta.data) == "Phenotype")] <- "Phénotype"
    colnames(singlet@meta.data)[which(colnames(singlet@meta.data) == "Phenotype.fine")] <- "Phénotype.fine"
    metadata <- metadata[which(rownames(metadata) %in% rownames(singlet@meta.data)),]
    singlet@meta.data <- cbind(singlet@meta.data,metadata)
    
    ## -- Percent -- ##
    DefaultAssay(singlet) <- "RNA"
    singlet[["percent.rb.meta"]] <- PercentageFeatureSet(singlet, pattern = "^RP[SL]")
    singlet[["percent.ig.meta"]] <- PercentageFeatureSet(singlet, pattern = "^IG")
    singlet[["percent.mt.meta"]] <- PercentageFeatureSet(singlet, pattern = "^MT-")
    DefaultAssay(singlet) <- "integrated"
    
    ## -- DE : Réponse RC/RP + RCHOP/Excipient -- ##
    DE <- noig(singlet)
    Idents(DE) <- "Reponse" ; 
    RNA_DE_Reponse <- FindMarkers(DE, ident.1 = "RC", ident.2 = "RP", logfc.threshold = 0.25, assay = "RNA", slot = "data", test.use = "negbinom")
    singlet@tools$RNA_DE_Reponse <- RNA_DE_Reponse %>% rownames_to_column(var="gene") %>% left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% as_tibble(RNA_DE_Reponse)
    
    DE <- PrepSCTFindMarkers(DE)
    SCT_DE_Reponse <- FindMarkers(DE, ident.1 = "RC", ident.2 = "RP", logfc.threshold = 0.25, assay = "SCT", verbose = F)
    singlet@tools$SCT_DE_Reponse <- SCT_DE_Reponse %>% rownames_to_column(var="gene") %>% left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% as_tibble(SCT_DE_Reponse) 
    
    if(condition=="Pre-greffe-Excipient"){
      singlet@tools$RNA_DE_PE <- seurat_DE(DE, "Condition", "Pre-greffe", "Excipient", 0.25, "RNA")
      singlet@tools$SCT_DE_PE <- seurat_DE(DE, "Condition", "Pre-greffe", "Excipient", 0.25, "SCT")
    }
    
    if(condition=="Post-greffe"){
      singlet@tools$RNA_DE_RE <- seurat_DE(DE, "Condition", "RCHOP", "Excipient", 0.15, "RNA")
      singlet@tools$SCT_DE_RE <- seurat_DE(DE, "Condition", "RCHOP", "Excipient", 0.15, "SCT")
    }
    
    if(condition=="All"){
      singlet@tools$RNA_DE_PE <- seurat_DE(DE, "Condition", "Pre-greffe", "Excipient", 0.25, "RNA")
      singlet@tools$SCT_DE_PE <- seurat_DE(DE, "Condition", "Pre-greffe", "Excipient", 0.25, "SCT")
      
      singlet@tools$RNA_DE_RE <- seurat_DE(DE, "Condition", "RCHOP", "Excipient", 0.15, "RNA")
      singlet@tools$SCT_DE_RE <- seurat_DE(DE, "Condition", "RCHOP", "Excipient", 0.15, "SCT")
    } ; rm(DE); gc();gc();gc();gc();gc();gc()
    
  
    
    ## -- Signature RCHOP -- ##
    Singature_RCHOP <- c("BAX","RPS27L","PHPT1","SRSF3","PSMB4","HIST1H2BK","TRIM22","AEN","PVT1","FDXR","BBC3","SRSF2","MRFAP1","DDB2","EIF2S3","MDM2","HNRNPH1","CHI3L2","CCNG1","LY86","ZMAT3","TXNIP","CD70","SNHG8","P4HA1","CDKN1A","ISG15")
    #Singature_RCHOP <- c("PPIA","RPS19","HIST1H2BK","BAX","RPS27L","SRSF3","PHPT1","PSMB4","TRIM22","VPREB3","BBC3","SNRPF","EIF2S3","CRIP1","FDXR","PSMB4","AEN","HSPA1A") 
    #Singature_RCHOP <- c("PHPT1","HIST1H2BK", "RPS27L","SRSF3","TRIM22","BBC3","FDXR","PSMB4","AEN","HSPA1A","BAX")
    if(phénotype!="Other"){
      singlet <- AddModuleScore(singlet,features = list(Singature_RCHOP), ctrl=100, replace =T, name="RCHOP_SCT_Score", assay = "SCT")
      singlet <- AddModuleScore(singlet,features = list(Singature_RCHOP), ctrl=100, replace =T, name="RCHOP_RNA_Score", assay = "RNA")
      singlet <- AddModuleScore(singlet,features = list(Singature_RCHOP), ctrl=100, replace =T, name="RCHOP_i_Score", assay = "integrated")
    }
    
    # Clean -------------------------------------------------------------------
    if(phénotype=="Other"){
      singlet@meta.data <- singlet@meta.data[,c("Sample","State","Condition","Phénotype","Phénotype.fine","Heavy","Light",
                                                "Clonotype","Chaine_Light","V_Light","D_Light","J_Light","C_Light","CDR3_Light","CDR3_nt_Light","Chaine_Heavy", "V_Heavy","D_Heavy","J_Heavy","C_Heavy","CDR3_Heavy","CDR3_nt_Heavy",
                                                "Phénotype_TRUST","V_Heavy_TRUST","D_Heavy_TRUST","J_Heavy_TRUST","Heavy_TRUST","CDR3_DNA_Heavy_TRUST","CDR3_AA_Heavy_TRUST","V_Light_TRUST","D_Light_TRUST","J_Light_TRUST","Light_TRUST","CDR3_DNA_Light_TRUST","CDR3_AA_Light_TRUST",
                                                "BCL2_L23L","BCL2_K22K","CD79B_Y196H","EZH2_Y646C","EZH2_Y646F","EZH2_Y646H","EZH2_Y646N","EZH2_Y646S","EZH2_A692V","EZH2_A682G","Phase","Reponse","Greffe","S.Score","G2M.Score", "RCHOP_AUC_Score",
                                                "percent.mt","percent.rb","percent.ig","percent.mt.meta","percent.rb.meta","percent.ig.meta","nCount_RNA","nFeature_RNA","nCount_SCT","nFeature_SCT", "nCount_HTO","nFeature_HTO", "integrated_snn_res.0.6")]
    }else{
      singlet@meta.data <- singlet@meta.data[,c("Sample","State","Condition","Phénotype","Phénotype.fine","Heavy","Light",
                                                "Clonotype","Chaine_Light","V_Light","D_Light","J_Light","C_Light","CDR3_Light","CDR3_nt_Light","Chaine_Heavy", "V_Heavy","D_Heavy","J_Heavy","C_Heavy","CDR3_Heavy","CDR3_nt_Heavy",
                                                "Phénotype_TRUST","V_Heavy_TRUST","D_Heavy_TRUST","J_Heavy_TRUST","Heavy_TRUST","CDR3_DNA_Heavy_TRUST","CDR3_AA_Heavy_TRUST","V_Light_TRUST","D_Light_TRUST","J_Light_TRUST","Light_TRUST","CDR3_DNA_Light_TRUST","CDR3_AA_Light_TRUST",
                                                "BCL2_L23L","BCL2_K22K","CD79B_Y196H","EZH2_Y646C","EZH2_Y646F","EZH2_Y646H","EZH2_Y646N","EZH2_Y646S","EZH2_A692V","EZH2_A682G","Phase","Reponse","Greffe","S.Score","G2M.Score", "RCHOP_AUC_Score", "RCHOP_RNA_Score1", "RCHOP_SCT_Score1", "RCHOP_i_Score1", 
                                                "percent.mt","percent.rb","percent.ig","percent.mt.meta","percent.rb.meta","percent.ig.meta","nCount_RNA","nFeature_RNA","nCount_SCT","nFeature_SCT", "nCount_HTO","nFeature_HTO", "integrated_snn_res.0.6")]
    }
    
    save(singlet, file = paste0(path.meta,"All_",condition,"_",phénotype,".RData")) ; rm(singlet); gc(); gc(); gc(); gc(); gc(); gc(); gc(); gc();
  }
}


# Correlation -------------------------------------------------------------
load(file = paste0(path.meta,"All_Pré-greffe-Excipient_All.RData"))
load(file = paste0(path.meta,"All_RCHOP_sign.RData"))
load(file = paste0(path.meta,"All_RCHOP_B.RData"))

Condition <- "Condition"
item1 <- "RCHOP"
item2 <- "Excipient"
fc_seuil <- 0.15

DE <- noig(singlet) ; Idents(DE) <- Condition ; 

#FindConservedMarkers
result <- FindConservedMarkers(DE, ident.1 = item1, ident.2 = item2 , grouping.var = "Sample", logfc.threshold = fc_seuil , assay = "RNA", slot = "data") ; annotations <- read.csv(paste0(path,"document/annotation.csv"))
result <- result %>% rownames_to_column(var="gene") %>% left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% as_tibble(result) %>% dplyr::select("gene", "max_pval", "minimump_p_val", "description")

#FindMarkers
response <- FindMarkers(DE, ident.1 = item1, ident.2 = item2 , logfc.threshold = fc_seuil, assay = "RNA", slot = "data")

gene <- intersect(rownames(response),result$gene)
length(intersect(gene, rownames(response)[which(response$avg_log2FC<0)]))

avg <- as.data.frame(log1p(AverageExpression(DE, verbose = FALSE)$RNA)) ; avg$gene <- rownames(avg)
p <- ggplot(avg, aes(RCHOP, Excipient)) + geom_point() + geom_point(data=avg[gene,], color='red',size=2) + stat_smooth(formula =y ~ x) + ggtitle("B Cells : log(1+ Average Normalised Counts)") 
LabelPoints(plot = p, points = gene, repel = T)





# Figure ------------------------------------------------------------------
# QC
load(file = paste0(path,"result/analyse_meta/All_All_All.RData"))

# Nombre de cellule
mpg <- as.data.frame(table(singlet$Sample, singlet$Condition)) ; colnames(mpg) <- c("Patient","Condition","Frequence")
mpg %>% mutate(Patient = fct_reorder(mpg$Patient, mpg$Frequence, .fun='min')) 
ggplot(mpg, aes(reorder(Patient, Frequence), Frequence, color=Condition)) + geom_point(aes(fill=Condition), size = 5) +
  labs(title="\nCell number per condition per patient", x="Patient", y="Cell number") + geom_hline(yintercept=400, color = "red") +
  theme_classic() +
  theme(plot.title = element_text(size = rel(1.2), hjust = 0.5), axis.title = element_text(size = rel(1))) 

# Nombre de reads
mpg <- as.data.frame(singlet[[c('Sample', 'Condition', 'nCount_RNA', 'nFeature_RNA')]]) ; 
colnames(mpg) <- c("Patient","Condition", "NbReads", "NbGenes")
mpg$Patient <- as.factor(mpg$Patient)
mpg$NbReads <- as.numeric(mpg$NbReads)
mpg$NbGenes <- as.numeric(mpg$NbGenes)

ggplot(mpg, aes(reorder(Patient,NbReads), NbReads, fill=Condition)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=T) + 
  ggtitle("\nRead number per condition per patient") + xlab("Patient") + ylab("Read number") +
  theme_classic() +
  theme(plot.title = element_text(size = rel(1.2), hjust = 0.5), axis.title = element_text(size = rel(1))) 

ggplot(mpg, aes(reorder(Patient,NbGenes), NbGenes, fill=Condition)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=T) + 
  ggtitle("\nExpressed genes per condition per patient") + xlab("Patient") + ylab("Expressed genes number") +
  theme_classic() +
  theme(plot.title = element_text(size = rel(1.2), hjust = 0.5), axis.title = element_text(size = rel(1))) 






mpg %>% mutate(Patient = fct_reorder(mpg$Patient, mpg$Frequence, .fun='min')) 
ggplot(mpg, aes(reorder(Patient, Frequence), Frequence, color=Condition)) + geom_point(aes(fill=Condition), size = 5) + theme_classic() +
  labs(title="Cell number by patient by condition", x="Patient", y="Cell number") + geom_hline(yintercept=400, color = "red")

# Nombre de gene
mpg <- as.data.frame(table(singlet$Sample, singlet$Condition, singlet$nCount_RNA, singlet$nFeature_RNA)) ; colnames(mpg) <- c("Patient","Condition","Frequence", "NbReads", "NbGenes")
mpg %>% mutate(Patient = fct_reorder(mpg$Patient, mpg$Frequence, .fun='min')) 
ggplot(mpg, aes(reorder(Patient, Frequence), Frequence, color=Condition)) + geom_point(aes(fill=Condition), size = 5) + theme_classic() +
  labs(title="Cell number by patient by condition", x="Patient", y="Cell number") + geom_hline(yintercept=400, color = "red")


# Transcriptomic
Idents(singlet) <- "Condition" 
DefaultAssay(singlet) <- "RNA"

Seurat::DimPlot(singlet, group.by = "Condition", label.size = 5, pt.size = 0.5, reduction = "umap", label = F, dims = c(1, 2)) #& Seurat::NoLegend()
Seurat::DimPlot(singlet, group.by =c("Sample","Condition","Phénotype","Phénotype.fine","Phase","Greffe","seurat_clusters","SCT_snn_res.0.6"), label.size = 5, pt.size = 1, reduction = "umap", label = T) & Seurat::NoLegend()
Seurat::DimPlot(object = singlet, group.by = "Phénotype", label.size = 5, pt.size = 0.7, reduction = 'umap', label = T) & Seurat::NoLegend() #& xlab(label = paste0(input$Reduction, " / PCA 1 : ", round(Seurat::Stdev(r$dataset[["pca"]])[1],2), " %")) & ylab(label = paste0("PCA 2 : ", round(Seurat::Stdev(r$dataset[["pca"]])[2],2), " %"))

Seurat::VlnPlot(singlet, idents = c("RC","RP"), features = c("BAX", "LDHA",  "RPS19"), ncol = 3, assay = 'RNA', slot = "data") & theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
Seurat::VlnPlot(singlet, idents = c("Pré-greffe","Excipient"), features = c("Greffe_AUC_Score"), ncol = 1, assay = 'RNA', slot = "count") & theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
Seurat::VlnPlot(singlet, group.by = "Sample", split.by = "Condition", features = c("RCHOP_AUC_Score"),ncol = 1,pt.size = 0.1) & theme(plot.title = element_text(size=10)) 
Seurat::VlnPlot(singlet, idents = c(item1,item2), features = c("signature_SCT1", "signature_RNA1", "signature_i1", "Reponse_AUC_Score"), ncol = 3, assay = 'RNA', slot = "count") & theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")

Seurat::FeaturePlot(singlet, features = c("percent.rb.meta","percent.mt.meta"), reduction = "umap", pt.size = 0.7, combine = T , ncol = 2) & Seurat::NoAxes()
Seurat::FeaturePlot(singlet, features = "percent.mt.meta", reduction = "umap", pt.size = 0.7, combine = T , ncol = 1) & Seurat::NoAxes()
Seurat::FeaturePlot(singlet, features = c("LDHA"), reduction = "umap", pt.size = 0.7, combine = T , ncol = 1) & Seurat::NoAxes()
Seurat::FeaturePlot(singlet, features = "nFeature_SCT", reduction = "umap", pt.size = 0.7, combine = T , ncol = 1) & Seurat::NoAxes()
Seurat::FeaturePlot(singlet, features = c("MS4A1", "CD3E"), reduction = "umap", pt.size = 0.7, combine = T , ncol = 2) & Seurat::NoAxes()
Seurat::FeaturePlot(singlet, features = "percent.ig", reduction = "umap", pt.size = 0.7, combine = T , ncol = 2) & Seurat::NoAxes()
Seurat::FeaturePlot(singlet, features = c("percent.rb", "percent.ig", "nCount_RNA","nFeature_RNA","nCount_SCT","nFeature_SCT","nCount_HTO","nFeature_HTO"), reduction = "umap", pt.size = 1, combine = T ) & Seurat::NoAxes()



# Trajectory
erythroid.cds <- as.cell_data_set(singlet)
erythroid.cds <- cluster_cells(cds = erythroid.cds, reduction_method = "UMAP")
erythroid.cds <- learn_graph(erythroid.cds, use_partition = TRUE)
erythroid.cds <- order_cells(erythroid.cds, reduction_method = "UMAP", root_cells = colnames(singlet)[which(singlet@meta.data$Phénotype.fine=="Progenitor cells")])
plot_cells(cds = erythroid.cds,color_cells_by = "Condition",cell_size = 1,show_trajectory_graph = TRUE, label_groups_by_cluster=FALSE,reduction_method = "UMAP",label_leaves=TRUE,label_branch_points=TRUE,graph_label_size=5, group_label_size = 4)


# first figure presentation app
load(file = paste0(path.meta,"All_0.RData"))
sunburst <- function(singlet){
  met <- data.frame(
    projet = rep("FLINOVO", length(singlet@meta.data$Sample)),
    patient = singlet@meta.data$Sample,
    condition = singlet@meta.data$Condition,
    phénotype = singlet@meta.data$Phénotype,
    sub = singlet@meta.data$Phénotype.fine,
    value = rep(1, length(singlet@meta.data$Sample)), stringsAsFactors = FALSE
  )
  
  df <- data.frame(ids=character(), labels=character(), parents=character(), values=integer(),stringsAsFactors=FALSE)
  
  df0 <- met %>% mutate(path = paste(projet, sep=";")) %>% dplyr::select(path, value)
  df1 <- met %>% mutate(path = paste(projet, patient, sep=";")) %>% dplyr::select(path, value)
  df2 <- met %>% mutate(path = paste(projet, patient, condition,  sep=";")) %>% dplyr::select(path, value)
  df3 <- met %>% mutate(path = paste(projet, patient, condition, phénotype, sep=";")) %>% dplyr::select(path, value)
  df4 <- met %>% mutate(path = paste(projet, patient, condition, phénotype, sub, sep=";")) %>% dplyr::select(path, value)
  
  
  for (dfX in c(df0,df1,df2,df3,df4)) {
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
  df$values <- c(table(df0)[,1],table(df1)[,1],table(df2)[,1],table(df3)[,1],table(df4)[,1])
  
  singlet@tools$sunburst <-  df
  return(singlet)
}
data <- sunburst(singlet)
sunburst <- data@tools$sunburst

result <- tidyr::crossing(singlet[["Sample"]], singlet[["Condition"]])
tables <- table(data.frame(Sample = singlet[["Sample"]],Condition = singlet[["Condition"]]))
for(Condition in names(table(singlet[["Condition"]]))){for(Sample in names(table(singlet[["Sample"]]))){result$Number[which(result$Sample == Sample & result$Condition == Condition)] <- tables[Sample,Condition]}
}

save(sunburst, result, file="/home/boris/Bureau/scShiny/inst/app/www/presentation.RData")
# Update ------------------------------------------------------------------
for(condition in c("All", "Post-greffe", "RCHOP", "Excipient","Pré-greffe","Pré-greffe-Excipient")){
  for(phénotype in c("All", "Other","B","BT")){
    load(file = paste0(path.meta,"All_",condition,"_",phénotype,".RData"))
    save(singlet, file = paste0(path.meta,"All_",condition,"_",phénotype,".RData"))
  }
}