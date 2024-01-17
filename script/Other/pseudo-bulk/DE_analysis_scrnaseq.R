source(file = "/home/boris/Bureau/scShiny/script/work_tools.R")
setwd(dir = "/home/boris/Bureau/scShiny/script/DE_analysis_scrnaseq/")
load(file="/home/boris/Bureau/scShiny/www/DE.RData")

#############################################################################################################################################################################################################################
## -- Load libraries -- ##
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(SingleCellExperiment)
library(DEGreport)
library(scPipe)
library(edgeR)
library(ggVennDiagram)
x <- 4000
DE_bulk <- function(seurat){
  ## -- Construction de l'object SCE -- ##
  sce <- SingleCellExperiment(assays = list(counts = seurat@assays$RNA@counts), colData = seurat@meta.data)  # Create single cell experiment object
  #assays(sce) # Explore the raw counts for the dataset : Check the assays present
  #dim(counts(sce)) ; counts(sce)[1:6, 1:6] # Count :    33382  20476 
  #dim(colData(sce)) ; head(colData(sce))   # Metadata : 33382     4
  
  
  ## -- Acquiring necessary metrics for aggregation across cells in a sample -- ##
  # Shape : clusters / samples
  kids <- purrr::set_names(levels(sce$cluster_id)) ; #kids # Named vector of cluster names
  nk <- length(kids) ; #nk                                 # Total number of clusters
  sids <- purrr::set_names(levels(sce$sample_id)) ; #sids  # Named vector of sample names
  ns <- length(sids) ; #ns                                 # Total number of samples 
  
  # Generate sample level metadata
  n_cells <- as.numeric(table(sce$sample_id)) ; #table(sce$sample_id)  # Determine the number of cells per sample
  m <- match(sids, sce$sample_id)                                     # Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
  ei <- data.frame(colData(sce)[m, ], n_cells, row.names = NULL) %>% select(-"cluster_id") ; # View(ei) # Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
  
  # Calculate quality control (QC) metrics + filters outlier
  sce <- addPerCellQC(sce)
  sce$is_outlier <- isOutlier(metric = sce$total,nmads = 2, type = "both", log = TRUE) # Get cells w/ few/many detected genes
  sce <- sce[, !sce$is_outlier] ; #dim(sce)                  # Remove outlier cells  : 20476 31556
  sce <- sce[rowSums(counts(sce) > 1) >= 10, ] ; #dim(sce)   # Remove lowly expressed genes which have less than 10 cells with any counts : 10966 31556
  
  
  #############################################################################################################################################################################################################################
  ## -- Count aggregation to sample level -- ##
  # Aggregate the counts per sample_id and cluster_id
  groups <- colData(sce)[, c("cluster_id", "sample_id")]                  # Identify groups for aggregation of counts : Subset metadata to only include the cluster and sample IDs to aggregate across
  pb <- Matrix.utils::aggregate.Matrix(t(counts(sce)),groupings = groups, fun = "sum") # Aggregate across cluster-sample groups
  splitf <- sapply(stringr::str_split(rownames(pb),pattern = "_",n = 2),`[`, 1)
  #class(pb); dim(pb); pb[1:6, 1:6]
  
  # Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
  pb <- split.data.frame(pb,factor(splitf)) %>% lapply(function(u)set_colnames(t(u),stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
  #class(pb) ; str(pb) ; options(width = 100) ; table(sce$cluster_id, sce$sample_id)   # Print out the table of cells in each cluster-sample group
  
  # Get sample names for each of the cell type clusters
  get_sample_ids <- function(x){ pb[[x]] %>% colnames() } # prep. data.frame for plotting
  de_samples <- map(1:length(kids), get_sample_ids) %>% unlist()
  
  # Get cluster IDs for each of the samples
  samples_list <- map(1:length(kids), get_sample_ids)
  get_cluster_ids <- function(x){rep(names(pb)[x],each = length(samples_list[[x]]))}
  de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>% unlist()
  
  # Create a data frame with the sample IDs, cluster IDs and condition
  gg_df <- data.frame(cluster_id = de_cluster_ids,sample_id = de_samples)
  gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id")]) 
  metadata <- gg_df %>% dplyr::select(cluster_id, sample_id, group_id) 
  clusters <- levels(as.factor(metadata$cluster_id))
  
  # Subsetting dataset to cluster(s) of interest
  cluster_metadata <- metadata[which(metadata$cluster_id == clusters[1]), ]                          # Subset the metadata to only the B cells
  rownames(cluster_metadata) <- cluster_metadata$sample_id                                           # Assign the rownames of the metadata to be the sample IDs
  counts <- pb[[clusters[1]]]                                                                        # Subset the counts to only the B cells
  cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
  all(rownames(cluster_metadata) == colnames(cluster_counts))                                        # Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
  
  return(list(counts=cluster_counts, metadata=cluster_metadata))
}
DES <- function(cluster){
  #############################################################################################################################################################################################################################
  cluster_counts <- cluster$counts
  cluster_counts <- cluster_counts[rownames(cluster_counts)[which(!str_detect(rownames(cluster_counts), "^IG"))],]
  
  cluster_metadata <- cluster$metadata
  cluster_metadata$patient <- c("A","A","B","B","C","C","D","D","E","E")
  
  
  ## --  Differential gene expression with DESeq2 : RC/RP -- ##
  # DESeq2 Object
  dds <- DESeqDataSetFromMatrix(cluster_counts, colData = cluster_metadata, design = ~ patient+group_id)
  rld <- rlog(dds, blind=TRUE)                                                                       # Transform counts for data visualization
  DESeq2::plotPCA(rld, intgroup = "group_id")
  
  # Extract the rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id"), drop=F])     # Hierarchical clustering : Plot heatmap
  
  # Run DESeq2 + Résult
  dds <- DESeq(dds)                                                             # Run DESeq2 differential expression analysis
  plotDispEsts(dds)                                                             # Plot dispersion estimates
  levels(cluster_metadata$group_id)[2] ; levels(cluster_metadata$group_id)[1]   # Output results of Wald test for contrast for stim vs ctrl
  
  contrast <- c("group_id", levels(cluster_metadata$group_id)[2], levels(cluster_metadata$group_id)[1]) ; resultsNames(dds)
  res <- results(dds,contrast = contrast,alpha = 0.05)
  res <- lfcShrink(dds,contrast =  contrast, res=res, type='normal')
  res_tbl <- res %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble() ; res_tbl
  #write.csv(res_tbl, paste0("results/", clusters[1], "_", levels(cluster_metadata$sample)[2], "_vs_", levels(cluster_metadata$sample)[1], "_all_genes.csv"),quote = FALSE, row.names = FALSE)
  
  padj_cutoff <- 0.05                                                                       
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>% dplyr::arrange(padj) ; sig_res  # Subset the significant results
  normalized_counts <- counts(dds, normalized = TRUE)
  
  top20_sig_genes <- sig_res %>% dplyr::arrange(padj) %>% dplyr::pull(gene) %>% head(n=20)  # Order results by padj values
  top20_sig_norm <- data.frame(normalized_counts) %>%  rownames_to_column(var = "gene") %>% dplyr::filter(gene %in% top20_sig_genes)
  gathered_top20_sig <- top20_sig_norm %>% gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")
  gathered_top20_sig <- inner_join(ei[, c("sample_id", "group_id" )], gathered_top20_sig, by = c("sample_id" = "samplename"))
  
  ## plot using ggplot2
  ggplot(gathered_top20_sig) +
    geom_point(aes(x = gene, y = normalized_counts,color = group_id), position=position_jitter(w=0.1,h=0)) +  
    scale_y_log10() + xlab("Genes") + ylab("log10 Normalized Counts") + ggtitle("Top 20 Significant DE Genes") + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5))
  
  sig_norm <- data.frame(normalized_counts) %>% rownames_to_column(var = "gene") #%>%dplyr::filter(gene %in% sig_res$gene) # Extract normalized counts for only the significant genes
  
  # Run pheatmap using the metadata data frame for the annotation
  pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], color = brewer.pal(6, "YlOrRd") , cluster_rows = T, show_rownames = F, annotation = cluster_metadata[, c("group_id", "cluster_id")], 
           border_color = NA, fontsize = 10, scale = "row", fontsize_row = 10, height = 20)
  
  # Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
  res_table_thres <- res_tbl %>% mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58)
  
  # Volcano plot
  ggplot(res_table_thres) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) + ggtitle("Volcano plot of stimulated B cells excipient to pregraft") +
    xlab("log2 fold change") + ylab("-log10 adjusted p-value") + scale_x_continuous(limits = c(-2.5,2.5)) + scale_y_continuous(limits = c(0,6)) +
    theme(legend.position = "none", plot.title = element_text(size = rel(1.5), hjust = 0.5), axis.title = element_text(size = rel(1.25)))                    
  
  return(sig_res)
}

#############################################################################################################################################################################################################################
## -- Chargement des données -- ##
seurat <- get(load(file = "/home/boris/Bureau/Flinovo/result/analyse_meta/All_RCHOP_sign.RData")) 
seurat <- load(file = "/home/boris/Documents/lipinskib/boris/Flinovo/integration/analyse_meta (copie)/All_RCHOP_sign.RData")

seurat <- singlet

## -- Formatage -- ##
sample_id <- paste0(seurat@meta.data$Sample, seurat@meta.data$Condition) #Condition
patient_id <- seurat@meta.data$Sample
group_id <- seurat@meta.data$Condition #Reponse
cluster_id <- seurat@meta.data$Phenotype <- "B cells"

seurat@meta.data <- seurat@meta.data[1]
seurat@meta.data$patient_id <- factor(patient_id)
seurat@meta.data$sample_id <- factor(sample_id)
seurat@meta.data$cluster_id <- factor(cluster_id)
seurat@meta.data$group_id <- factor(group_id)
seurat@meta.data <- seurat@meta.data[-1]


#############################################################################################################################################################################################################################
## -- Signature Excipient / RCHOP -- ##
## -- Differential gene expression with EdgeR : RCHOP/Excipient -- ##
cluster <- DE_bulk(seurat)
cluster_counts <- cluster$counts
cluster_counts <- cluster_counts[rownames(cluster_counts)[which(!str_detect(rownames(cluster_counts), "^IG"))],]
cluster_metadata <- cluster$metadata



y <- DGEList(counts=as.matrix(cluster_counts), group = c(1,1,2,2,3,3,4,4,5,5,6,6,7,7))

Subject <- factor(cluster_metadata$patient_id)
Treat <- factor(cluster_metadata$group_id, levels=c("Excipient","RCHOP"))
design <- model.matrix(~Subject+Treat)
y <- estimateDisp(y,design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
res<-topTags(qlf, n = 242)$table
View(res$table)

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "http://www.ensembl.org")
genes <- biomaRt::getBM(filters = "hgnc_symbol", attributes = c("hgnc_symbol","entrezgene_id"), values = rownames(qlf), mart = mart)
for (i in 1:length(rownames(qlf))) {
  if(rownames(qlf)[i] %in% genes$hgnc_symbol){
    if(!is.na(genes$entrezgene_id[which(genes$hgnc_symbol == rownames(qlf)[i] )])){rownames(qlf)[i] <-  genes$entrezgene_id[which(genes$hgnc_symbol == rownames(qlf)[i] )]} 
  }
}

go <- goana(qlf, species="Hs") ; topGO(go, n=15)
keg <- kegga(qlf, species="Hs") ; topKEGG(keg, n=15)


## Bulk
signature_RE_bulk <- rownames(topTags(qlf, n = 242)$table)
count <- cluster_counts[which(rownames(cluster_counts) %in% signature_RE_bulk),]
count <- count[, c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)]
pheatmap(count[ , 1:length(colnames(count))], color = brewer.pal(6, "YlOrRd") , cluster_rows = T, show_rownames = F, annotation = cluster_metadata[, c("patient_id", "group_id", "cluster_id")], cluster_cols = F, border_color = NA, fontsize = 10, scale = "row", fontsize_row = 10, height = 20)

count <- cluster_counts[which(rownames(cluster_counts) %in% signature_RE_bulk),]
pheatmap(count[ , 1:length(colnames(count))], color = brewer.pal(6, "YlOrRd") , cluster_rows = T, show_rownames = F, annotation = cluster_metadata[, c("patient_id", "group_id", "cluster_id")], cluster_cols = F, border_color = NA, fontsize = 10, scale = "row", fontsize_row = 10, height = 20)

load(file="/home/boris/Bureau/scShiny/script/DE/signature.Rdata")
save(cluster_counts, cluster_metadata ,res,signature_RE_bulk,signature_RE_sc,file="/home/boris/Bureau/scShiny/www/DE.RData")

## PCA
#seurat <- get(load(file = "/home/boris/Documents/analyse/ALL_B_Post_greffe.Rdata")) 
singlet <- SCTransform(singlet, method = "glmGamPoi", verbose = F, ncells = 10000, conserve.memory = T) 
singlet <- RunPCA(singlet, verbose = F, feature = signature_RE_bulk, approx=F) %>% RunUMAP(dims = 1:5, verbose = F) %>% FindNeighbors(dims = 1:5, verbose = F) %>% FindClusters(verbose = F) %>% RunTSNE(dims = 1:5, verbose = F)
Idents(singlet) <- "Reponse" ; plot_dimension(singlet)



## Single-cell
#signature_RE_sc <- rownames(singlet@tools$DE_RE)
#Idents(singlet)<-"Condition"
#cell <- c(sample(colnames(singlet)[which(singlet$Condition %in% "Excipient")], x), sample(colnames(singlet)[which(singlet$Condition %in% "RCHOP")], x))
#DoHeatmap(singlet, features = signature_RE_sc, cell = cell, size = 3, assay = 'SCT', slot = "scale.data")

## Intersect
signature_RCHOP <- Reduce(intersect, list(signature_RE_sc, signature_RE_bulk))
sign <- list(Bulk = signature_RE_bulk,Single_cell = signature_RE_sc)
ggVennDiagram(sign, label_alpha = 1, label = "count") + scale_fill_gradient(low="white",high = "white") + theme(legend.position = "none") + ggtitle("Expression différentielle :\nRCHOP vs Excipient") + theme(plot.title = element_text(hjust = 0.5))


#############################################################################################################################################################################################################################
## -- Signature RC / RP -- ##
## -- Differential gene expression with DESeq2 : RC/RP -- ##
cluster <- DE_bulk(seurat)
result <- DES(cluster)

## Bulk
signature_RPRC_bulk <- sig_res$gene
sig_norm <- data.frame(normalized_counts) %>% rownames_to_column(var = "gene") %>%dplyr::filter(gene %in% signature_RPRC_bulk) # Extract normalized counts for only the significant genes
pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], color = brewer.pal(6, "YlOrRd") , cluster_rows = T, show_rownames = F, annotation = cluster_metadata[, c("group_id", "cluster_id")], border_color = NA, fontsize = 10, scale = "row", fontsize_row = 10, height = 20)

## PCA
#signature_Reponse <- Reduce(intersect, list(signature_RPRC_bulk, signature_RPRC_sc))
seurat <- get(load(file = "/home/boris/Documents/analyse/ALL_B_PG.Rdata")) 
singlet <- SCTransform(singlet, method = "glmGamPoi", verbose = F, ncells = 10000, conserve.memory = T) 
singlet <- RunPCA(singlet, verbose = F, feature = signature_RPRC_bulk, approx=F) %>% RunUMAP(dims = 1:5, verbose = F) %>% FindNeighbors(dims = 1:5, verbose = F) %>% FindClusters(verbose = F) %>% RunTSNE(dims = 1:5, verbose = F)
Idents(singlet) <- "Reponse" ; plot_dimension(singlet)

## Single-cell
#signature_RPRC_sc <- rownames(singlet@tools$DE_RC_RP)[1:150]
#Idents(singlet)<-"Reponse"
#cell <- c(sample(colnames(singlet)[which(singlet$Condition %in% "RC")], x), sample(colnames(singlet)[which(singlet$Condition %in% "RP")], x))
#DoHeatmap(singlet, features = signature_RPRC_sc, cell = cell, size = 3, assay = 'SCT', slot = "scale.data")

## Intersect
#signature_Reponse <- Reduce(intersect, list(signature_RPRC_bulk, signature_RPRC_sc))
#sign <- list(Bulk = signature_RPRC_bulk, Single_cell = signature_RPRC_sc)
#ggVennDiagram(sign, label_alpha = 1, label = "count") + scale_fill_gradient(low="white",high = "white") + theme(legend.position = "none") + ggtitle("Expression différentielle :\nRéponse complète vs Réponse partielle") + theme(plot.title = element_text(hjust = 0.5))



#############################################################################################################################################################################################################################
#AddModuleScore(singlet,features = list(signature_RE_bulk), name="RE")
#FeaturePlot(singlet, features = "RE1", reduction = "tsne", pt.size = 1, combine = T, label = TRUE, repel = TRUE) & Seurat::NoAxes() & scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))
#DimPlot(object = singlet, label.size = 5, pt.size = 1, reduction = 'tsne', label = TRUE) & theme(legend.position = "top") & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4)))
