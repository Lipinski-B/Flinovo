library("muscat")

## -- Chargement des données -- ##
seurat <- get(load(file = "/home/boris/Documents/analyse/All_Post_greffe_B.Rdata")) 

## -- Formatage -- ##
sample_id <- paste0(seurat@meta.data$orig.ident, seurat@meta.data$Condition) #Condition
patient_id <- seurat@meta.data$orig.ident
group_id <- seurat@meta.data$Condition #Reponse
cluster_id <- seurat@meta.data$Phénotype

seurat@meta.data <- seurat@meta.data[1]
seurat@meta.data$patient_id <- factor(patient_id)
seurat@meta.data$sample_id <- factor(sample_id)
seurat@meta.data$cluster_id <- factor(cluster_id)
seurat@meta.data$group_id <- factor(group_id)
seurat@meta.data <- seurat@meta.data[-1]

sce <- SingleCellExperiment(assays = list(counts = seurat@assays$RNA@counts), colData = seurat@meta.data)  # Create single cell experiment object


# remove undetected genes
sce <- sce[rowSums(counts(sce) > 0) > 0, ]

# Calculate quality control (QC) metrics + filters outlier
sce <- addPerCellQC(sce)
sce$is_outlier <- isOutlier(metric = sce$total,nmads = 2, type = "both", log = TRUE) # Get cells w/ few/many detected genes
sce <- sce[, !sce$is_outlier] ; #dim(sce)                  # Remove outlier cells  : 20476 31556
sce <- sce[rowSums(counts(sce) > 1) >= 10, ] ; #dim(sce)   # Remove lowly expressed genes which have less than 10 cells with any counts : 10966 31556

# compute sum-factors & normalize
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)

# or
#library(sctransform)
#assays(sce)$vstresiduals <- vst(counts(sce), show_progress = FALSE)$y

sce <- prepSCE(sce, 
               kid = "cluster_id", # subpopulation assignments
               gid = "group_id",  # group IDs (ctrl/stim)
               sid = "sample_id",   # sample IDs (ctrl/stim.1234)
               drop = TRUE) #


# compute UMAP using 1st 20 PCs
sce <- runUMAP(sce, pca = 20)


# compute pseudobulks (sum of counts)
pb <- aggregateData(sce, 
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))

# one sheet per subpopulation
assayNames(pb)
# pseudobulks for 1st subpopulation
t(head(assay(pb)))




(pb_mds <- pbMDS(pb))
# use very distinctive shaping of groups & change cluster colors
pb_mds <- pb_mds + 
  scale_shape_manual(values = c(17, 4)) +
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Set2"))
# change point size & alpha
pb_mds$layers[[1]]$aes_params$size <- 5
pb_mds$layers[[1]]$aes_params$alpha <- 0.6
pb_mds






# run pseudobulk (aggregation-based) DS analysis
res <- pbDS(pb, method = "edgeR")

# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)
# view results for 1st cluster
k1 <- tbl[[1]]
head(format(k1[, -ncol(k1)], digits = 2))


# construct design & contrast matrix
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts("RCHOP-Excipient", levels = mm)

# run DS analysis
pbDS(pb, design = mm, contrast = contrast)


# filter FDR < 5%, abs(logFC) > 1 & sort by adj. p-value
tbl_fil <- lapply(tbl, function(u) {
  u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
  dplyr::arrange(u, p_adj.loc)
})

# nb. of DS genes & % of total by cluster
n_de <- vapply(tbl_fil, nrow, numeric(1))
p_de <- format(n_de / nrow(sce) * 100, digits = 3)
data.frame("#DS" = n_de, "%DS" = p_de, check.names = FALSE)

# view top 2 hits in each cluster
top2 <- bind_rows(lapply(tbl_fil, top_n, 2, p_adj.loc))
format(top2[, -ncol(top2)], digits = 2)
