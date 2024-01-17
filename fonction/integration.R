## -- Intégration -- ##
load(file = paste0(path.meta,"All_0.RData"))

singlet <- seurat_subset(singlet, "Condition", c("RCHOP", "Excipient"))
singlet <- seurat_subset(singlet, "Phenotype", "B cells"); 

singlet <- SplitObject(singlet, split.by = "Sample")
singlet <- singlet[c("FL09C1164", "FL140304", "FL02G095", "FL12C1888", "FL08G0431", "FL08G0404")] #FL08G0293 FL06G1206 FL05G0330

for (i in c("FL09C1164", "FL140304", "FL02G095", "FL12C1888", "FL08G0431", "FL08G0404")) {
  singlet[[i]] <- NormalizeData(singlet[[i]], verbose = TRUE)
  singlet[[i]] <- SCTransform(singlet[[i]], assay = "RNA", method = "glmGamPoi", conserve.memory = T, ncells = 10000, variable.features.n = 3000)
} 
integ_features <- SelectIntegrationFeatures(object.list = singlet, nfeatures = 3000) 
singlet <- PrepSCTIntegration(object.list = singlet, anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = singlet, normalization.method = "SCT", anchor.features = integ_features)
singlet <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT", k.weight = 100 ,features.to.integrate = integ_features, dims = 1:15)

save(singlet, file = paste0(path.meta,"All_RCHOP_sign.RData"))