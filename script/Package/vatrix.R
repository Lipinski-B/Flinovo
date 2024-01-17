# Vartrix -----------------------------------------------------------------
library(Seurat)
library(Matrix)
library(stringr)
library(ggplot2)
library(gplots)

setwd(dir = "/home/boris/Documents/lipinskib/boris/mutect2/")

# Read in the sparse genotype matrix
snv <- readMM("result/Vartrix/Filtered/FL14_variants.txt")
bc <- read.table("data/bam/barcodes.tsv", header = F)
snps <- read.table("result/Vartrix/Unfiltered/SNV.loci.txt", header = F)
dimnames(snv) <- list(snps$V1,bc$V1) ; 
dim(snv) ; snv <- snv[,colSums(snv != 0) != 0] ; dim(snv)
dim(snv) ; snv <- snv[rowSums(snv != 0) != 0,] ; dim(snv)
snv <- as.matrix(snv)

col <- c("No Call","ref/ref","alt/alt","alt/ref")
result <- as.data.frame(setNames(replicate(length(col),numeric(0), simplify = F),col ))

for (i in 1:nrow(snv)) {
  if(length(table(snv[i,]))==4){result[i,] <- table(snv[i,])}
  if(length(table(snv[i,]))==3){result[i,] <- c(table(snv[i,]),0)}  
  if(length(table(snv[i,]))==2){result[i,] <- c(table(snv[i,]),0,0)}
  if(length(table(snv[i,]))==1){result[i,] <- c(table(snv[i,]),0,0,0)}
}

rownames(result) <- rownames(snv)
result <- as.matrix(result) ; rm(snv, col, snps, bc) ; gc(); gc(); gc(); gc(); gc(); gc(); gc();
View(result)
saveRDS(result, file = "/home/boris/Documents/lipinskib/boris/mutect2/result/Vartrix/Filtered/result.Rds")

result <- readRDS(file = "/home/boris/Documents/lipinskib/boris/mutect2/result/Vartrix/Filtered/result.Rds")
target <- 'chr21:8393292'
hotspot <- snv[target,]

load(file = paste0("/home/boris/Bureau/Flinovo/result/analyse_patient/", patient,"/Post-greffe_B_", patient,".RData"))
singlet <- AddMetaData(object = singlet, metadata = hotspot, col.name = "hotspot")
UMAPPlot(object = singlet, pt.size = 1, group.by = "hotspot")
