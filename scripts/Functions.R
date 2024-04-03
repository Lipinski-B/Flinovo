library(ggplot2)
library(limma)
library(edgeR)
library(RColorBrewer)
library(gplots)
library(genefilter)
library(LPS)
library(viridis)
library(fastICA)
library(umap)
library(Rtsne)
library(magrittr)
library(stringr)
library(Seurat)
library(SingleCellExperiment)
library(scuttle)
library(purrr)
library(dplyr)
library(pheatmap)
library(scater)
library(tibble)
library(tidyr)
library(ggVennDiagram)
library("readxl")
library(tidyverse)
library(ggridges)
library(grid)
library(fgsea)
library(org.Hs.eg.db)
library(biomaRt)
library(NMF)
library(DT)


suppressPackageStartupMessages(library(escape))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dittoSeq))

options(digits = 3)

GS <- getGeneSets(library = "H")
#C5 <- getGeneSets(library = "C5")

# BP <- getGeneSets(library = "C5", subcategory = "BP")
# CC <- getGeneSets(library = "C5", subcategory = "CC")
# MF <- getGeneSets(library = "C5", subcategory = "MF")

# Log2(FC) = 0.25 ; FC = 2^0.25 = 1.189207
# log2(FC) = 0.15 ; FC = 2^0.15 = 1.109569 
# means that the expression of that gene is increased in WT relative to KO by a multiplicative factor of 2^log(FC) ≈ FC

# Function ----------------------------------------------------------------
ATACseq <- function(dge_list = omni, TSSfilter = 2,substracting_background = FALSE, noise_prob = 0.1,filter_by_expr = FALSE,min.count = 10, min.sample = 3,
                     filter_top_peaks_by_sample = FALSE, n_top_peaks = 10000, lib.size = "final reads", # colSums or "final reads"
                     norm.package = 'edgeR', norm_on_high_average_peaks = FALSE, norm = "TMM", # or TMMwsp, upperquartile, RLE
                     norm.p = 0.75, norm.logratiotrim = 0.3, norm.sumtrim = 0.05 # default 0,75, 0.3, 0.05
                    ){
  
  #dge_list <- dge_list[,dge_list$samples$TSS]
  dge_list$samples <- dge_list$samples[dge_list$samples$TSS>= TSSfilter,] # removing low TSS 
  dge_list$samples$group <- as.character(dge_list$samples$group) # Renaming groups

  # Subtract background
  if (substracting_background == TRUE ) {     message('substracting background')
    for(i in 1:ncol(dge_list)) {
      size= round(dge_list$samples$`Final`[i])# * (1 - dge_list$samples$`Fraction of reads in called peak regions [2]`[i]))
      dge_list$samples$noiseLevel[i] <- qbinom(p = noise_prob, size= size, prob=500/2176355877, lower.tail=FALSE)
      dge_list$counts[,i] <- pmax(dge_list$counts[,i] - dge_list$samples$noiseLevel[i], 0L)
    }
  }
  
  # Filter genes with min.count read counts in min.prop patients
  if (filter_by_expr == TRUE) {               message ('filtering peaks with min count')
    filtered <- filterByExpr(dge_list$counts, min.count=min.count,  group = dge_list$samples$group)
    dge_list <- dge_list[filtered, ]
  }
  
  # Filter top peaks by sample
  if (filter_top_peaks_by_sample == TRUE) {   message ('filtering top peaks by sample')
    merged_top_peaks <- vector()
    for ( i in 1:ncol(dge_list$counts)) {
      ordered_peaks <- order(dge_list$counts[,i])
      top_peaks_by_sample <- ordered_peaks[1:n_top_peaks]
      merged_top_peaks <- append(merged_top_peaks, top_peaks_by_sample)
      merged_top_peaks <- unique(merged_top_peaks)
    }
    dge_list <- dge_list[merged_top_peaks, ]
  }

  # Reassign lib size
  message('Reassign lib size')
  if (lib.size == "final reads") {dge_list$samples$lib.size <- dge_list$samples$`Final`}
  if (lib.size == "colSums") {    dge_list$samples$lib.size <- colSums(dge_list$counts)}
  
  # Normalize
  if (norm.package == "edgeR"){ # Recompute normalization factors TMM par edgeR
    if (norm_on_high_average_peaks ==TRUE) {  message('Recompute normalization factors TMM par edgeR excluding low average genes')
      dge_list$genes$rowmeans <- aveLogCPM(cpm(dge_list))
      filter <- which(dge_list$genes$rowmeans > -2)
      dge_list2 <- calcNormFactors(dge_list[filter,], method = norm , logratioTrim = norm.logratiotrim, sumTrim = norm.sumtrim, p = norm.p)
      dge_list$samples$norm.factors <- dge_list2$samples$norm.factors
    }else{                                    message('Recompute normalization factors TMM par edgeR')
      dge_list <- calcNormFactors(dge_list, method = norm , logratioTrim = norm.logratiotrim, sumTrim = norm.sumtrim, p = norm.p)
    }
  }}

PB     <- function(singlet, target="Condition"){
  seurat <- singlet
  sample_id <- paste0(seurat@meta.data$Sample, seurat@meta.data[, target]) #Condition
  Sample <- seurat@meta.data$Sample
  Condition <- seurat@meta.data[, target] #Reponse
  cluster_id <- seurat@meta.data$Phenotype <- "B cells"
  
  seurat@meta.data <- seurat@meta.data[1]
  seurat@meta.data$Sample <- factor(Sample)
  seurat@meta.data$sample_id <- factor(sample_id)
  seurat@meta.data$cluster_id <- factor(cluster_id)
  seurat@meta.data$Condition <- factor(Condition)
  seurat@meta.data <- seurat@meta.data[-1]
  
  # Construction de l'object SCE
  sce <- SingleCellExperiment(assays = list(counts = seurat@assays$RNA@counts), colData = seurat@meta.data)  # Create single cell experiment object
  kids <- purrr::set_names(levels(sce$cluster_id)) ; nk <- length(kids)
  sids <- purrr::set_names(levels(sce$sample_id)) ; ns <- length(sids)
  
  # Generate sample level metadata
  n_cells <- as.numeric(table(sce$sample_id)) ; #table(sce$sample_id)  # Determine the number of cells per sample
  m <- match(sids, sce$sample_id)                                     # Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
  ei <- data.frame(colData(sce)[m, ], n_cells, row.names = NULL) %>% dplyr::select(-"cluster_id") ; # View(ei) # Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
  
  # Calculate quality control (QC) metrics + filters outlier
  sce <- addPerCellQC(sce)
  sce$is_outlier <- isOutlier(metric = sce$total,nmads = 2, type = "both", log = TRUE) # Get cells w/ few/many detected genes
  sce <- sce[rowSums(counts(sce) > 1) >= 10, !sce$is_outlier] ; # Remove outlier cells  : 20476 31556  #dim(sce)   # Remove lowly expressed genes which have less than 10 cells with any counts : 10966 31556
  
  # Aggregate the counts per sample_id and cluster_id
  groups <- colData(sce)[, c("cluster_id", "sample_id")]                  # Identify groups for aggregation of counts : Subset metadata to only include the cluster and sample IDs to aggregate across
  pb <- Matrix.utils::aggregate.Matrix(t(counts(sce)),groupings = groups, fun = "sum") # Aggregate across cluster-sample groups
  splitf <- sapply(stringr::str_split(rownames(pb),pattern = "_",n = 2),`[`, 1)
  pb <- split.data.frame(pb,factor(splitf)) %>% lapply(function(u)set_colnames(t(u),stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+"))) # Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
  
  # Get sample names for each of the cell type clusters
  get_sample_ids <- function(x){ pb[[x]] %>% colnames() } # prep. data.frame for plotting
  de_samples <- map(1:length(kids), get_sample_ids) %>% unlist()
  
  # Get cluster IDs for each of the samples
  samples_list <- map(1:length(kids), get_sample_ids)
  get_cluster_ids <- function(x){rep(names(pb)[x],each = length(samples_list[[x]]))}
  de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>% unlist()
  
  # Create a data frame with the sample IDs, cluster IDs and condition
  metadata <- data.frame(cluster_id = de_cluster_ids,sample_id = de_samples)
  metadata <- left_join(metadata, ei[, c("sample_id", "Condition")])
  metadata <- metadata %>% dplyr::select(cluster_id, sample_id, Condition)
  clusters <- levels(as.factor(metadata$cluster_id))
  
  # Subsetting dataset to cluster(s) of interest
  cluster_metadata <- metadata[which(metadata$cluster_id == clusters[1]), ]                          # Subset the metadata to only the B cells
  rownames(cluster_metadata) <- cluster_metadata$sample_id                                           # Assign the rownames of the metadata to be the sample IDs
  counts <- pb[[clusters[1]]]                                                                        # Subset the counts to only the B cells
  cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
  all(rownames(cluster_metadata) == colnames(cluster_counts))                                        # Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
  cluster_metadata$Sample <- sub("[^0-9]+$", "\\1", cluster_metadata$sample_id)
  
  # Calculate quality control (QC) metrics + filters outlier
  cutoff <- edgeR::filterByExpr(cluster_counts, group=metadata$Condition)
  cluster_counts <- cluster_counts[cutoff,]
  
  return(list("counts"=cluster_counts, "metadata"=cluster_metadata))
}

EdgeR.Object <- function(counts, metadata, group, gene = NULL){
  y <- DGEList(counts = counts, samples = metadata, group = group, gene = NULL)
  y <- calcNormFactors(y, method = "TMM")
  return(y)
}
EdgeR.glmTreat  <- function(y, lfc = 1.2, FDR = 0.05, protocol){
  design <- model.matrix(protocol, data=y$samples)
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTrendedDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  
  # Test for Differential Expression Relative to a Threshold : much better than glmQLFTest
  fit <- glmQLFit(y, design) 
  tr  <- glmTreat(fit, lfc = log2(lfc)) #tr <- glmQLFTest(fit)
  result <- EdgeR.Resuls(y, tr, lfc/2, FDR = FDR)
  
  return(result)
}
EdgeR.glmTest   <- function(y, lfc = 1.2, FDR = 0.05, protocol){
  design <- model.matrix(protocol, data=y$samples)
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTrendedDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  
  # Test for Differential Expression
  fit <- glmQLFit(y, design) 
  tr  <- glmQLFTest(fit) #tr <- glmQLFTest(fit)z
  result <- EdgeR.Resuls(y, tr, lfc/2, FDR = FDR)
  
  return(result)
}
EdgeR.exactTest <- function(y, lfc = 1.2, FDR = 0.05 ){
  y  <- estimateDisp(y)
  tr <- exactTest(y)
  result <- EdgeR.Resuls(y, tr, lfc/2, FDR = FDR)
  
  return(result)
}
EdgeR.Resuls    <- function(y, tr, FC = 0.6, FDR = 0.05){
  tr <- tr[!is.na(tr$table$PValue),]
  tr$table$FDR <- p.adjust(tr$table$PValue, method="BH") ; 
  tr$table$threshold <- tr$table$FDR < 0.05 ; 
  tr$table[tr$table$logFC > 0 & tr$table$FDR < 0.05, "Expression"] <- "Up" #tr$table[tr$table$logFC > FC & tr$table$FDR < 0.05, "Expression"] <- "Up"
  tr$table[tr$table$logFC < 0 & tr$table$FDR < 0.05, "Expression"] <- "Down"#tr$table[tr$table$logFC < -FC & tr$table$FDR < 0.05, "Expression"] <- "Down"
  tr$table[tr$table$FDR >= 0.05,  "Expression"] <- "Not Sig"
  #tr$table[tr$table$logFC > -FC & tr$table$logFC < FC, "Expression"] <- "Not Sig"
  tr$table$genes <- tr$genes$symbol
  if(!is.null(y$genes$symbol)){tr$table$genes <- tr$genes$symbol} else {tr$table$genes <- rownames(tr$table)}
  
  result <- c()
  result <- topTags(tr, n = sum(tr$table$FDR < FDR))
  num <- c("logFC", "logCPM","PValue","FDR")   #num <- c("logFC","unshrunk.logFC", "logCPM","PValue","FDR")
  result$table[,num] <- signif(result$table[,num],3)

  return(list("y"= y,"tr" = tr, "result" = result)) #, "go" = go, "ko" = ko))
}

EdgeR.Workflow  <- function(counts, metadata, group, gene = NULL, protocol, FC = 0.6, custom = F, FDR = 0.05, treat = T, glm = T, exact = T ){
  
  object  <- EdgeR.Object(counts, metadata, group, gene = NULL); if(custom==T){object <- object[,which(!object$samples$sample_id == "FL05G0330RCHOP")]}
  result  <- list()
  
  if (treat == T) {result$treat <- EdgeR.glmTreat (object, lfc = FC, FDR = FDR, protocol = protocol)}
  if (glm == T)   {result$glm   <- EdgeR.glmTest  (object, lfc = FC, FDR = FDR, protocol = protocol)}
  if (exact == T) {result$exact <- EdgeR.exactTest(object, lfc = FC, FDR = FDR)}
  
  result$y  <- result$exact$y
  result$tr <- result$exact$tr
  
  return(result)
}

noig      <- function(DE){
  DE@assays$RNA@data <- DE@assays$RNA@data[rownames(DE@assays$RNA@data)[which(!str_detect(rownames(DE@assays$RNA@data), "IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]"))],];
  DE@assays$RNA@scale.data <- DE@assays$RNA@scale.data[rownames(DE@assays$RNA@scale.data)[which(!str_detect(rownames(DE@assays$RNA@scale.data), "IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]"))],];
  DE@assays$RNA@counts <- DE@assays$RNA@counts[rownames(DE@assays$RNA@counts)[which(!str_detect(rownames(DE@assays$RNA@counts), "IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]"))],];
  return(DE)
}
result    <- function(data){
  response <- as.data.frame(data)
  response <- response %>% rownames_to_column(var="gene") %>% left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% as_tibble(response) 
  response1 <- signif(as.data.frame(response[,c("PValue", "logFC", "FDR")]),16)
  response <- cbind(response1, response[,c("gene", "description")])
  response[response$logFC > log2(FC) & response$FDR < 0.05, "Expression"] <- "Up"
  response[response$logFC < -log2(FC) & response$FDR < 0.05, "Expression"] <- "Down"
  response[response$FDR >= 0.05,  "Expression"] <- "Not Sig"
  response[response$logFC > -log2(FC) & response$logFC < log2(FC), "Expression"] <- "Not Sig"
  response <- response[,c("gene","PValue","logFC","FDR", "Expression", "description")]
  return(response)
}

# Figure ------------------------------------------------------------------
Clustering   <- function(data, meta = c("Condition","Sample")) {
  data <- data[['y']]
  rld.cor <- cor(edgeR::cpm(data, log = T, prior.count = 1), method = "pearson")
  pheatmap(rld.cor, #color = RColorBrewer::brewer.pal(data$samples[,c("Condition","Sample")], "Blues"), 
           main = "Clustering using log(TMM+1) pics counts. \nDistance: Pearson correlation.",
           annotation = data$samples[,meta], 
           border_color=NA, fontsize = 10, fontsize_row = 10, height=20)
}
Volcano_plot <- function(data, FC=2, xscale = 5, yscale=4, seuil=0.6, label=NULL) {
  data <- data[['tr']]$table
  
  if(length(unique(data$Expression))==2){color = c("grey","#F8766D")}else{color = c("#619CFF","grey","#F8766D")}
  
  if(is.null(label)){filter <- data %>% filter(-log10(data$FDR)>1.2 & abs(data$logFC)>FC)}else{filter <- data %>% filter(data$genes %in% label & -log10(data$FDR)>1.2 & abs(data$logFC)>FC) %>% arrange(FDR)}
  
  filter <- filter %>% arrange(FDR)
  
  ggplot(data) +
    geom_point(aes(x = logFC, y = -log10(FDR), colour = Expression)) + 
    ggrepel::geom_label_repel(data=filter[1:20,], aes(label=genes, x = logFC, y = -log10(FDR))) +
    ggtitle("\nVolcano plot WT vs KO : FDR < 0.05 & logFC > 1.2\n") + xlab("log2 FC") + ylab("-log10 FDR") +
    scale_x_continuous(limits = c(-xscale,xscale)) + scale_y_continuous(limits = c(0,yscale)) + scale_color_manual(values = color) + 
    geom_hline(yintercept=1.30103, col="black", linetype="dashed") + geom_vline(xintercept=c(-seuil,seuil), col="black", linetype="dashed") +
    theme_classic() + theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), axis.title = element_text(size = rel(1))) 
}
MD_plot      <- function(data, xscale = 6, yscale = 5) {
  data <- data[['tr']]
  if(length(unique(data$table$Expression))==2){color = c("grey","#F8766D")}
  else{color = c("#619CFF","grey","#F8766D")}
  ggplot(data$table) +
    geom_point(aes(x = data$AveLogCPM, y = logFC, colour = Expression)) + 
    #ggrepel::geom_label_repel(data=data %>% filter(abs(data$table$logFC)>1.2), aes(label=data$table$genes, x = data$AveLogCPM, y = data$table$logFC)) +
    ggtitle("\nMean-Difference Plot of Pics Data : \nWT vs KO FDR < 0.05 & logFC > 1.2\n") + xlab("Average log CPM") + ylab("-log10 FC") + 
    scale_x_continuous(limits = c(-1,xscale)) + scale_y_continuous(limits = c(-yscale,yscale)) + scale_color_manual(values = color) + 
    geom_hline(yintercept=c(0.6,-0.6), col="black", linetype="dashed") +
    theme_classic() + theme(plot.title = element_text(size = rel(1.2), hjust = 0.5), axis.title = element_text(size = rel(1))) 
}
Top20        <- function(data,n=20) {
  counts.norm <- edgeR::cpm(data[['y']], log=TRUE, prior.count = 1)
  top_sig <- data[['tr']]$table
  FC <- top_sig %>% rownames_to_column(var = "symbol") %>% dplyr::arrange(FDR)  %>% dplyr::pull(logFC, symbol) %>% head(n=n) #  Order results by padj values
  #FC <- FC[which(names(FC) %in% label)]
  top_sig <- data.frame(counts.norm) %>% rownames_to_column(var = "genes") %>% dplyr::filter(genes %in% names(FC))
  FC <- FC[order(factor(names(FC), levels = top_sig$genes))]
  
  if(!is.null(data[['y']]$genes$symbol)){
    names(FC) <- top_sig$genes <- data[['tr']]$genes[names(FC),"symbol"]
    FC <- FC[!is.na(names(FC))]
    top_sig <- top_sig[!is.na(top_sig$genes),]
  }
  
  top_sig <- top_sig %>% gather(colnames(top_sig)[2:length(colnames(top_sig))], key = "Sample", value = "counts.norm")
  top_sig$Condition <- rep(data[['y']]$samples$Condition, each = length(names(FC)))
  top_sig$FC <- rep(FC,length(unique(top_sig$Sample)))
  
  ggplot(top_sig) +
    geom_point(aes(reorder(x = genes, FC, max), y = counts.norm, color = Condition), cex = 1.5, position=position_jitter(w=0.1,h=0)) +  
    xlab("Genes") + ylab("TMM log Normalized Counts") + ggtitle("Top Significant") + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=15))
}
TopHeatmap   <- function(data, FDR=0.01, FC=0.5, row = T, cluster_cols = T){
  annot = as.data.frame(data[['y']]$samples[,c("Sample", "Condition")])
  signature <- rownames(data[["tr"]]$table)[data[["tr"]]$table$FDR<FDR & abs(data[["tr"]]$table$logFC)>FC]
  data <- data[["y"]]$counts
  tmm <- edgeR::cpm(data, log=TRUE, prior.count = 1)[signature,]
  pheatmap(tmm, show_rownames = row, scale = "row", main = "Heatmap using log(TMM+1) read counts\n", 
           annotation_col = annot, cluster_cols = cluster_cols)
}
PCA_edger    <- function(data, var = "Condition", position = "bottom", labels = "Sample"){
  data <- data[['y']]
  col <- data$samples[[var]]
  #col <- stringr::str_replace(col , as.character(unique(col)[1]), c("#619CFF"))
  #col <- stringr::str_replace(col , as.character(unique(col)[2]), c("#F8766D"))
  
  col <- brewer.pal(8, "Dark2")[1:length(unique(col))]
  
  limma::plotMDS(data, top = 7700, labels = data$samples[[labels]], col=col, pch=20, cex = 0.8) 
  legend(position , col = unique(col), as.character(unique(data$samples[[var]])), pch=1, cex = 0.6)
}
Top_pic      <- function(data, pic){
  result <- data[['result']]$table
  rownames(result) <- paste0(rownames(result),"_",result[,'genes'])
  fc <- result[which(result[,'genes'] %in% pic),c('symbol','logFC')]
  gene <- rownames(fc)
  meta <- data.frame(Sample = rep(data[['y']]$samples[,"Sample"],each=length(gene)),
                     Condition = rep(data[['y']]$samples[,"Condition"],each=length(gene)),
                     FC = rep(fc$logFC,4))
  
  if(length(gene)==1){data <- data.frame(t(data[['y']]$counts[gsub('(.*)_\\w+', '\\1', gene),]))
  }else{               data <- data.frame(data[['y']]$counts[gsub('(.*)_\\w+', '\\1', gene),])}
  
  rownames(data) <- gene
  data <- data %>% rownames_to_column(var = "Pics") %>% gather(colnames(data), key = "Sample", value = "Counts") 
  data$Condition <- meta$Condition
  data$FC <- meta$FC
  
  
  ggplot(data, aes(x = Pics, y = Counts)) + 
    ggtitle(paste(pic, collapse = ", ")) + #geom_line(aes(group=Sample)) + 
    geom_point(data = data, aes(reorder(x = Pics, FC, max),color = Condition), cex = 3) + 
    theme_bw() + theme(plot.title = element_text(size = rel(1.2), hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title = element_text(size = rel(1))) 
}

SC.Enrichissement <- function(singlet, ES, Condition){
  singlet <- AddMetaData(singlet, ES)
  singlet@meta.data$active.idents <- singlet@active.ident
  
  HM <- data.frame(singlet[[colnames(ES)]], Idents(singlet)) ; colnames(HM)[ncol(HM)] <- Condition
  
  HM.Result <- getSignificance(HM, group = Condition, fit = "Wilcoxon") 
  HM.Result$Hallmark <- rownames(HM.Result)
  HM.Result <- HM.Result[order(HM.Result$FDR),]
  HM.Result <- HM.Result[HM.Result$FDR<0.05,]
  HM.Result <- HM.Result %>% format(scientific=T) %>% as_tibble() %>% dplyr::select(-W.statistic) %>% relocate(Hallmark,.before = 1) 
  
  return(HM.Result)
}
SC.Volcano <- function(data, FC, xscale = 1, yscale=150, label=NULL,seuil=log2(FC)){
  if(length(unique(data$Expression))==2){color = c("grey","#F8766D")}else{color = c("#619CFF","grey","#F8766D")}
  filter <- data %>% filter(-log10(data$FDR)>1.2 & abs(data$logFC)>log2(FC))
  
  ggplot(data) +
    geom_point(aes(x = logFC, y = -log10(FDR), colour = Expression)) + 
    ggrepel::geom_label_repel(data=filter, aes(label=gene, x = logFC, y = -log10(FDR))) +
    ggtitle("\nVolcano plot WT vs KO : FDR < 0.05 & logFC > 1.2\n") + xlab("log2 FC") + ylab("-log10 FDR") +
    scale_x_continuous(limits = c(-xscale,xscale)) + scale_y_continuous(limits = c(0,yscale)) + scale_color_manual(values = color) + 
    geom_hline(yintercept=FC, col="black", linetype="dashed") + geom_vline(xintercept=c(-seuil,seuil), col="black", linetype="dashed") +
    theme_classic() + theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), axis.title = element_text(size = rel(1))) 
}
SC.Result <- function(response){
  annotations <- read.csv(paste0(path,"/Rapport-RNAseq/annotation.csv"))
  response <- response %>% rownames_to_column(var="gene") %>% left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% as_tibble(response) 
  response1 <- signif(as.data.frame(response[,c("p_val", "avg_log2FC", "p_val_adj")]),16)
  response <- cbind(response1, response[,c("gene", "pct.1", "pct.2", "description")])
  response <- response[,c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "description")]
  colnames(response) <- c("gene", "p_val", "logFC", "pct.1", "pct.2", "FDR", "description")
  
  response[response$logFC > log2(FC) & response$FDR < 0.05, "Expression"] <- "Up"
  response[response$logFC < -log2(FC) & response$FDR < 0.05, "Expression"] <- "Down"
  response[response$FDR >= 0.05,  "Expression"] <- "Not Sig"
  response[response$logFC > -log2(FC) & response$logFC < log2(FC), "Expression"] <- "Not Sig"
  return(response)
}

DT.table <- function(table){
  DT::datatable(table, class = 'cell-border stripe', rownames = F, 
                extensions = c("Buttons"), 
                selection = "none",
                options = list(dom = "tipfB",
                               select = list(style = 'os', items = 'row'),
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print','colvis')
                )) %>%
    DT::formatSignif(colnames(select_if(table, is.numeric))  , digits=3)
}

# Other -------------------------------------------------------------------
DESeq2 <- function(counts, metadata, FC = 0.6){
  y <- DESeqDataSetFromMatrix(counts, colData = metadata, design = ~ Sample+Condition)
  y <- DESeq(y)#, minReplicatesForReplace=Inf)
  y$samples$Condition <- y$Condition
  
  tr <- result <- list()
  tr$table <- results(y, contrast=c("Condition",unique(as.vector(y$samples$Condition))[2],unique(as.vector(y$samples$Condition))[1]), pAdjustMethod = "BH", alpha=0.05)#, cooksCutoff=FALSE, independentFiltering=FALSE)
  tr$table <- tr$table[complete.cases(tr$table),]
  tr$table <- tr$table[order(tr$table$padj),]
  tr$table$threshold <- tr$table$padj < 0.05
  tr$table <- as.data.frame(tr$table) %>% bind_cols(feature = rownames(tr$table)) %>% mutate(log_padj = - log(.data$padj, base = 10))
  
  colnames(tr$table) <- c("baseMean","logFC", "lfcSE","stat","PValue","FDR","threshold","genes","log_padj")
  tr$table[tr$table$logFC > FC & tr$table$FDR < 0.05, "Expression"] <- "Up"
  tr$table[tr$table$logFC < -FC & tr$table$FDR < 0.05, "Expression"] <- "Down"
  tr$table[tr$table$FDR >= 0.05 , "Expression"] <- "Not Sig"
  tr$table[tr$table$logFC > -FC & tr$table$logFC < FC, "Expression"] <- "Not Sig"
  
  result$table <- tr$table[tr$table$threshold %in% T,]
  
  return(return(list("y"= y,"tr"= tr, "result"=result)))
}
Ontology_Enrichment <- function(condition){
  GO.Up = as.data.frame(topGO(condition,sort="down"))
  GO_custom.Up <- data.frame(
    Term = rep(paste(GO.Up$Ont, ":", GO.Up$Term, "\n", GO.Up$N, "genes - FDR =", GO.Up$FDR.Up),2),
    Regulation = c(rep("up",length(GO.Up$Term)),rep("down",length(GO.Up$Term))),
    Count = c(GO.Up$Up,GO.Up$Down), 
    pValueUp = c(GO.Up$FDR.Up, GO.Up$FDR.Up),
    label_ypos = c(GO.Up$Up-3,GO.Up$Up+GO.Up$Down+2)
  )
  plot1 <- ggplot(data=GO_custom.Up, aes(x=reorder(Term, -pValueUp), y=Count, fill=Regulation), color = Regulation) +
    geom_bar(stat="identity", width = 0.7)+ xlab("GO term") + ylab("Number of DEGs up and down inside the up-regulated term.") + ggtitle("Up-regulated pathways")  +
    geom_text(aes(y=label_ypos, label=Count), color="black", size=6) + scale_fill_brewer(palette="Paired") + theme_minimal() + 
    theme(aspect.ratio = 3/2, axis.text=element_text(size=17), axis.title=element_text(size=15), plot.title = element_text(size=22)) + coord_flip()
  
  
  GO.Down = as.data.frame(topGO(condition,sort="up"))
  GO_custom.Down <- data.frame(
    Term = rep(paste(GO.Down$Ont, ":", GO.Down$Term, "\n", GO.Down$N, "genes - FDR =", GO.Down$FDR.Down),2),
    Regulation = c(rep("up",length(GO.Down$Term)),rep("down",length(GO.Down$Term))),
    Count = c(GO.Down$Up,GO.Down$Down), 
    pValueDown = c(GO.Down$FDR.Down, GO.Down$FDR.Down),
    label_ypos = c(GO.Down$Up-3,GO.Down$Up+GO.Down$Down+2)
  )
  plot2 <- ggplot(data=GO_custom.Down, aes(x=reorder(Term, -pValueDown), y=Count, fill=Regulation), color = Regulation) +
    geom_bar(stat="identity", width = 0.7) + xlab("GO term") + ylab("Number of DEGs up and down inside the down-regulated term.") + ggtitle("Down-regulated pathways") +
    geom_text(aes(y=label_ypos, label=Count),  color="black", size=6) + scale_fill_brewer(palette="Paired") + theme_minimal() + 
    theme(aspect.ratio = 3/2, axis.text=element_text(size=15), axis.title=element_text(size=15), plot.title = element_text(size=22)) + coord_flip()
  
  return(list(plot1, plot2))
}
KEGG_Enrichment <- function(condition){
  GO.Up = as.data.frame(topKEGG(condition,sort="down"))
  GO.Up = GO.Up[GO.Up$FDR.Up<0.6,]
  GO_custom.Up <- data.frame(
    Term = rep(paste(GO.Up$Pathway, "\n", GO.Up$N, "genes - FDR =", GO.Up$FDR.Up),2),
    Regulation = c(rep("up",length(GO.Up$Pathway)),rep("down",length(GO.Up$Pathway))),
    Count = c(GO.Up$Up,GO.Up$Down), 
    pValueUp = c(GO.Up$FDR.Up, GO.Up$FDR.Up),
    label_ypos = c(GO.Up$Up-0.3,GO.Up$Up+GO.Up$Down+0.3)
  )
  plot1 <- ggplot(data=GO_custom.Up, aes(x=reorder(Term, -pValueUp), y=Count, fill=Regulation), color = Regulation) +
    geom_bar(stat="identity", width = 0.7)+ xlab("KEGG term") + ylab("Number of DEGs up and down inside the up-regulated term.") + ggtitle("Up-regulated pathways")  +
    geom_text(aes(y=label_ypos, label=Count), color="black", size=6) + scale_fill_brewer(palette="Paired") + theme_minimal() +
    theme(aspect.ratio = 3/2, axis.text=element_text(size=17), axis.title=element_text(size=15), plot.title = element_text(size=22)) + coord_flip()
  
  
  GO.Down = as.data.frame(topKEGG(condition,sort="up"))
  GO.Down = GO.Down[GO.Down$FDR.Down<0.6,]
  GO_custom.Down <- data.frame(
    Term = rep(paste(GO.Down$Pathway, "\n", GO.Down$N, "genes - FDR =", GO.Down$FDR.Down),2),
    Regulation = c(rep("up",length(GO.Down$Pathway)),rep("down",length(GO.Down$Pathway))),
    Count = c(GO.Down$Up,GO.Down$Down), 
    pValueDown = c(GO.Down$FDR.Down, GO.Down$FDR.Down),
    label_ypos = c(GO.Down$Up-0.3,GO.Down$Up+GO.Down$Down+0.3)
  )
  plot2 <- ggplot(data=GO_custom.Down, aes(x=reorder(Term, -pValueDown), y=Count, fill=Regulation), color = Regulation) +
    geom_bar(stat="identity", width = 0.7) + xlab("KEGG term") + ylab("Number of DEGs up and down inside the down-regulated term.") + ggtitle("Down-regulated pathways") +
    geom_text(aes(y=label_ypos, label=Count),  color="black", size=6) + scale_fill_brewer(palette="Paired") + theme_minimal() +
    theme(aspect.ratio = 3/2, axis.text=element_text(size=15), axis.title=element_text(size=15), plot.title = element_text(size=22)) + coord_flip()
  
  return(list(plot1, plot2))
}

#mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#symbol <- getBM(values = rownames(counts), filters = "hgnc_symbol", mart = mart, attributes = c("entrezgene_id","hgnc_symbol"))

# Enrichment
# go <- goana(tr, species="Mm", geneid = tr$genes$entrezgene_id, FDR = 0.05)
# p.go <- c(go$P.Up, go$P.Down)
# lgo <- length(go$P.Up)
# 
# FDR <- p.adjust(p.go)
# go$FDR.Up <- signif(FDR[1:as.numeric(lgo)],3)
# go$FDR.Down <- signif(FDR[as.numeric(lgo+1) : as.numeric(lgo+lgo)],3)
# 
# ko <- kegga(tr, species="Mm", geneid = tr$genes$entrezgene_id, FDR = 0.05)
# p.ko <- c(ko$P.Up, ko$P.Down)
# lko <- length(ko$P.Up)
# 
# FDR <- p.adjust(p.ko)
# ko$FDR.Up <- signif(FDR[1:as.numeric(lko)],3)
# ko$FDR.Down <- signif(FDR[as.numeric(lko+1) : as.numeric(lko+lko)],3)


# - Nombre de DEGs : `r sum(tr$table$FDR < 0.05)` (FDR > 0.05 & logFC > 1.2)
# - Nombre de gènes up-régulés dans les KO : `r summary(decideTests(tr))[3]`
# - Nombre de gènes down-régulés dans les KO : `r summary(decideTests(tr))[1]`
# - logFC positif = up-régulé dans les KO / down régulé dans les WT
# - logFC négatif = up-régulé dans les WT / down régulé dans les KO

