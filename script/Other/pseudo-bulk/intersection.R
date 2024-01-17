knitr::opts_chunk$set(echo = TRUE)
path <- "/home/boris/Bureau/Flinovo/"
path.meta <- paste0(path,"result/analyse_meta/")
path.patient <- paste0(path,"result/analyse_patient/")
source(file = paste0(path,"fonction/library.R"))
source(file = paste0(path,"fonction/workflow.R"))
source(file = paste0(path,"fonction/aggregate.R"))
work="/home/boris/Documents/lipinskib/boris"
source(paste0(work,"/ATACseq/script/Functions_ATAC_v4.R"))
library(pheatmap)
library(scater)



# RCHOP -------------------------------------------------------------------
condition <- "Post-greffe"
phénotype <- "B" ;
load(file = paste0(path.meta,"All_",condition,"_",phénotype,".RData")) ; singlet <- noig(singlet) ; DefaultAssay(singlet) <- "RNA" ; Idents(singlet) <- "Condition"
sc <- PB(singlet)

sc.EdgeR <- RNAseq(sc[["cluster_counts"]], sc[["cluster_metadata"]], sc[["cluster_metadata"]]$Condition)
sc.EdgeR <- sc.EdgeR[,which(!sc.EdgeR$samples$sample_id == "FL05G0330RCHOP")]
sc.EdgeR <- EdgeR(sc.EdgeR, model.matrix(~ Sample+Condition, data=sc.EdgeR$samples))
signature.EdgeR.PB.RCHOP <- rownames(sc.EdgeR[["tr"]]$table[sc.EdgeR[["tr"]]$table$FDR < 0.05,])

Clustering(sc.EdgeR)
PCA_edger(sc.EdgeR)
result.RCHOP <- signif(as.data.frame(sc.EdgeR[['result']][,c("logFC", "PValue", "FDR")]),3)
DT::datatable(result.RCHOP, class = 'cell-border stripe')
signature.RCHOP <- rownames(result.RCHOP)
Volcano_plot(sc.EdgeR, xscale = 4.5, yscale = 6, FC=0.6, seuil=0.6)
MD_plot(sc.EdgeR, xscale = 15)
Top20(sc.EdgeR,20)
col <- c("FL02G095Excipient","FL08G0293Excipient", "FL05G0330Excipient","FL06G1206Excipient","FL08G0404Excipient","FL08G0431Excipient","FL09C1164Excipient","FL12C1888Excipient", "FL140304Excipient",
         "FL06G1206RCHOP","FL08G0293RCHOP","FL02G095RCHOP","FL08G0404RCHOP","FL08G0431RCHOP","FL09C1164RCHOP","FL12C1888RCHOP","FL140304RCHOP")
TopHeatmap(sc.EdgeR, FDR=0.01, FC=0.5, col = col, cluster_cols = F)
TopHeatmap(sc.EdgeR, FDR=0.1, FC=0.5, col = col, cluster_cols = F)

# norm_counts <- as.data.frame(sc[["cluster_counts"]])
# norm_counts$names <- norm_counts$description <- rownames(sc[["cluster_counts"]])
# norm_counts <- norm_counts[,c("names", "description",col)]
# write.table(norm_counts, file='/home/boris/Bureau/Flinovo/script/Signature/GSEA/RCHOP.gct', quote=FALSE, sep='\t', row.names = F) #1.2  8986    17



# Greffe ------------------------------------------------------------------
condition <- "Pre-greffe-Excipient"
phénotype <- "B" ; load(file = paste0(path.meta,"All_",condition,"_",phénotype,".RData")) ; singlet <- noig(singlet) ; DefaultAssay(singlet) <- "RNA" ; Idents(singlet) <- "Condition"
singlet@meta.data$Condition <- stringr::str_replace(singlet@meta.data$Condition, "Pre-greffe", "Pregreffe")

sc <- PB(singlet)
sc.EdgeR <- RNAseq(sc[["cluster_counts"]], sc[["cluster_metadata"]], sc[["cluster_metadata"]]$Condition)
sc.EdgeR <- EdgeR(sc.EdgeR, model.matrix(~ Sample+Condition, data=sc.EdgeR$samples))
signature.EdgeR.PB.Pregreffe <- rownames(sc.EdgeR[["tr"]]$table[sc.EdgeR[["tr"]]$table$FDR < 0.05,])

Clustering(sc.EdgeR)
PCA_edger(sc.EdgeR)
result.greffe <- signif(as.data.frame(sc.EdgeR[['result']][,c("logFC", "PValue", "FDR")]),3)
DT::datatable(result.greffe, class = 'cell-border stripe')
signature.greffe <- rownames(result.greffe)
Volcano_plot(sc.EdgeR, xscale = 4.5, yscale = 10, FC=1, seuil=0.6)
MD_plot(sc.EdgeR, xscale = 15)
Top20(sc.EdgeR, 30)
TopHeatmap(sc.EdgeR, FDR=0.001, FC=3.5, col = colnames(sc.EdgeR[['y']]))
TopHeatmap(sc.EdgeR, FDR=0.1, FC=0.5, row = F, col = colnames(sc.EdgeR[['y']]))


Seurat::DimPlot(object = singlet, group.by = "Condition", label.size = 5, pt.size = 1,reduction = "pca", label = T) & Seurat::NoLegend() & 
  xlab(label = paste0("PCA 1 : ", round(Seurat::Stdev(singlet[["pca"]])[1],2), " %")) & 
  ylab(label = paste0("PCA 2 : ", round(Seurat::Stdev(singlet[["pca"]])[2],2), " %"))

Seurat::DimPlot(object = singlet, group.by = "Sample", label.size = 5, pt.size = 1,reduction = "pca", label = T) & Seurat::NoLegend() & 
  xlab(label = paste0("PCA 1 : ", round(Seurat::Stdev(singlet[["pca"]])[1],2), " %")) & 
  ylab(label = paste0("PCA 2 : ", round(Seurat::Stdev(singlet[["pca"]])[2],2), " %"))



load("/home/boris/Bureau/Flinovo/script/Signature/intersect_diss.genes_ref")
load("/home/boris/Bureau/Flinovo/script/Signature/intersect_FL18-13_DE_Subtilisin")




# norm_counts <- as.data.frame(sc[["cluster_counts"]])
# norm_counts$names <- norm_counts$description <- rownames(sc[["cluster_counts"]])
# norm_counts <- norm_counts[,c("names", "description",
#                               "FL02G095Excipient","FL08G0293Excipient", "FL05G0330Excipient","FL06G1206Excipient","FL08G0404Excipient","FL08G0431Excipient","FL09C1164Excipient","FL12C1888Excipient", "FL140304Excipient",
#                               "FL06G1206Pregreffe","FL08G0293Pregreffe","FL02G095Pregreffe","FL08G0404Pregreffe","FL08G0431Pregreffe","FL09C1164Pregreffe","FL12C1888Pregreffe","FL140304Pregreffe", "FL05G0330Pregreffe")]
# write.table(norm_counts, file='/home/boris/Bureau/Flinovo/script/Signature/GSEA/Pregreffe.gct', quote=FALSE, sep='\t', row.names = F) #1.2  9469    18


# GSEA --------------------------------------------------------------------
library(fgsea)
library(data.table)
library(ggplot2)

counts <- sc[["cluster_counts"]]

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
symbol <- getBM(values = rownames(counts), filters = "hgnc_symbol", mart = mart, attributes = c("entrezgene_id","hgnc_symbol"))

test <- sc.EdgeR[["tr"]]$table
test <- cbind(test, symbol[, "entrezgene_id"][match(rownames(test), symbol[, "hgnc_symbol"])])
test <- test[!is.na(test$`symbol[, "entrezgene_id"][match(rownames(test), symbol[, "hgnc_symbol"])]`),]
test <- test[!duplicated(test$`symbol[, "entrezgene_id"][match(rownames(test), symbol[, "hgnc_symbol"])]`),]
rownames(test) <- test$`symbol[, "entrezgene_id"][match(rownames(test), symbol[, "hgnc_symbol"])]`
test$`symbol[, "entrezgene_id"][match(rownames(test), symbol[, "hgnc_symbol"])]` <- NULL

ranks <- test$logFC ; names(ranks) <- rownames(test) 

barplot(sort(test$logFC, decreasing = T))
pathways <- fgsea::gmtPathways("/home/boris/Bureau/Flinovo/script/Signature/c6.all.v7.5.1.symbols.gmt") #pathways.hallmark %>% head() %>% lapply(head)
fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize = 500)
head(fgseaRes[order(padj, -abs(NES)), ], n=10)
plotEnrichment(pathways[["PDGF_ERK_DN.V1_DN"]], ranks)

pathways <- fgsea::gmtPathways("/home/boris/Bureau/Flinovo/script/Signature/c7.all.v7.5.1.symbols.gmt") #pathways.hallmark %>% head() %>% lapply(head)
fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize = 500, nPermSimple = 10000)
head(fgseaRes[order(padj, -abs(NES)), ], n=10)
plotEnrichment(pathways[["NAKAYA_PBMC_FLUARIX_FLUVIRIN_AGE_18_50YO_7DY_DN"]], ranks)

topUp <- fgseaRes %>% filter(ES > 0) %>% top_n(10, wt=-padj)
topDown <- fgseaRes %>% filter(ES < 0) %>% top_n(10, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% arrange(-ES)
plotGseaTable(pathways[topPathways$pathway], ranks, fgseaRes, gseaParam = 0.5)

library('enrichplot')
gseaplot2(fgseaRes, geneSetID = 1, title = fgseaRes$Description[1])



GSEA = function(gene_list, GO_file, pval) {
  set.seed(54321)
  library(dplyr)
  library(gage)
  library(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = fgsea::gmtPathways(GO_file)
  
  fgRes <- fgsea::fgsea(pathways = myGO, 
                        stats = gene_list,
                        minSize=15,
                        maxSize=600,
                        nperm=10000) %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval)
  print(dim(fgRes))
  
  ## Filter FGSEA by using gage results. Must be significant and in same direction to keep 
  gaRes = gage::gage(gene_list, gsets=myGO, same.dir=TRUE, set.size =c(15,600))
  
  ups = as.data.frame(gaRes$greater) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
  
  downs = as.data.frame(gaRes$less) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
  
  print(dim(rbind(ups,downs)))
  ## Define up / down pathways which are significant in both tests
  keepups = fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ]
  keepdowns = fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ]
  
  fgRes = fgRes[ !is.na(match(fgRes$pathway, 
                              c( keepups$pathway, keepdowns$pathway))), ] %>% 
    arrange(desc(NES))
  fgRes$pathway = stringr::str_replace(fgRes$pathway, "GO_" , "")
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
  
  
  upcols =  colorRampPalette(colors = c("red4", "red1", "lightpink"))( sum(filtRes$Enrichment == "Up-regulated"))
  downcols =  colorRampPalette(colors = c( "lightblue", "blue1", "blue4"))( sum(filtRes$Enrichment == "Down-regulated"))
  colos = c(upcols, downcols)
  names(colos) = 1:length(colos)
  filtRes$Index = as.factor(1:nrow(filtRes))
  
  g = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_col( aes(fill = Index )) +
    scale_fill_manual(values = colos ) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="GSEA - Biological Processes") + 
    theme_minimal()
  
  
  output = list("Results" = fgRes, "Plot" = g)
  return(output)
}
out <- GSEA(ranks, "/home/boris/Bureau/Flinovo/script/Signature/c7.all.v7.5.1.symbols.gmt", 0.05)

# Single Cell -------------------------------------------------------------
Condition <- "Condition"
item1 <- "RCHOP"
item2 <- "Excipient"
fc_seuil <- 0.25

DE <- noig(singlet) ; Idents(DE) <- Condition ;

#FindConservedMarkers
result <- FindConservedMarkers(DE, ident.1 = item1, ident.2 = item2 , grouping.var = "Sample", logfc.threshold = fc_seuil , assay = "RNA", slot = "data") ; annotations <- read.csv(paste0(path,"document/annotation.csv"))
result <- result %>% rownames_to_column(var="gene") %>% left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% as_tibble(result) %>% dplyr::select("gene", "max_pval", "minimump_p_val", "description")

#FindMarkers
response <- FindMarkers(DE, ident.1 = item1, ident.2 = item2 , logfc.threshold = fc_seuil, assay = "RNA", slot = "data")
response2 <- FindMarkers(DE, ident.1 = item1, ident.2 = item2 , logfc.threshold = fc_seuil, assay = "SCT", slot = "data")

gene <- intersect(rownames(response),result$gene)



# Intersection ------------------------------------------------------------
library(ggVennDiagram)
x <- list(EgdeR.PB = signature.EdgeR.PB, DESeq2.PB = signature.DESeq2.PB, EgdeR.ABV = signature.EdgeR.ABV, DESeq2.ABV = signature.DESeq2.ABV,SC.intersect = rownames(response))
ggVennDiagram(x, category.names = c("EgdeR.PB","DESeq2.PB","EgdeR.ABV","DESeq2.ABV","SC.intersect"), label_alpha = 1, label = "count") +
  scale_fill_gradient(low="white",high = "white") + theme(legend.position = "none") +
  ggtitle("DE list intersection") + theme(plot.title = element_text(hjust = 0.5))


x <- list(EgdeR.PB = signature.EdgeR.PB, DESeq2.PB = signature.DESeq2.PB, EgdeR.ABV = signature.EdgeR.ABV, DESeq2.ABV = signature.DESeq2.ABV)
ggVennDiagram(x, category.names = c("EgdeR.PB","DESeq2.PB","EgdeR.ABV","DESeq2.ABV"), label_alpha = 1, label = "count") +
  scale_fill_gradient(low="white",high = "white") + theme(legend.position = "none") +
  ggtitle("DE list intersection") + theme(plot.title = element_text(hjust = 0.5))
