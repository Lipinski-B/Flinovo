library(SingleCellExperiment)
library(Seurat)
library(SCopeLoomR)
library(SCENIC)

load(file="/home/boris/Documents/analyse/singlet_FL09C1164.RData")
#cellInfo <- data.frame(seuratCluster=Idents(all))
loom <- build_loom("/home/boris/Bureau/scShiny/script/scenic/R/FL09.loom", dgem=all[["RNA"]]@counts, 
                   default.embedding = cbind( all@reductions[["pca"]]@cell.embeddings[,1], all@reductions[["pca"]]@cell.embeddings[,2]) , 
                   default.embedding.name = 'pca')

loom <- add_cell_annotation(loom, all@meta.data)

## -- Running SCENIC ------------------------------------------------------------------------------- ##
## Input
setwd("/home/boris/Bureau/scShiny/script/scenic")

# Expression matrix
#loomPath <- system.file(package="SCENIC", "R/FL09.loom")
#loom <- open_loom(loomPath)
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)
dim(exprMat) # FL09 : gene 17080  / cells : 3308

# Cell info/phenodata
cellInfo$nGene <- colSums(exprMat>0) #head(cellInfo)
cellInfo <- data.frame(cellInfo) #cbind(table(cellInfo$Phénotype))
dir.create("int") ; saveRDS(cellInfo, file="int/cellInfo.Rds")

colVars <- list(CellType=c("Astrocytes"="forestgreen","B-cells"="darkorange","CD4+ T-cells"="magenta4","CD8+ T-cells"="hotpink","DC"="red3","Fibroblasts"="skyblue",
                           "HSC"="darkblue","Macrophages"="grey","Mesangial cells"="green","Monocytes"="blue","Myocytes"="red","NK cells"="orange","Skeletal muscle"="black")) # Color to assign to the variables (same format as for NMF::aheatmap)
  
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$Phénotype)]
saveRDS(colVars, file="/home/boris/Bureau/scShiny/script/scenic/int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))


# Initialize SCENIC settings
data(defaultDbNames)
scenicOptions <- initializeScenic(org="hgnc" , dbDir="/home/boris/Bureau/scShiny/script/scenic/DB", dbs=defaultDbNames[["hgnc"]], datasetTitle="SCENIC FL09", nCores=12) 


# Co-expression network : Gene filter/selection
genesKept <- geneFiltering(exprMat, 
                           scenicOptions=scenicOptions,
                           minCountsPerGene=5*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
interestingGenes <- c("BCL2")
interestingGenes[which(!interestingGenes %in% genesKept)]# any missing?

exprMat_filtered <- exprMat[genesKept, ] 
dim(exprMat_filtered)  # FL09 : gene 6704  / cells : 3308
rm(exprMat)


# Correlation
runCorrelation(exprMat_filtered, scenicOptions)


# GENIE3
exprMat_filtered <- log2(exprMat_filtered+1) 
#runGenie3(exprMat_filtered, scenicOptions)

# GRNBoost (pour les gros datasets)
#exportsForArboreto(exprMat_filtered, scenicOptions, dir = "int") #+ script python


# Importing pySCENIC results : Read table
motifEnrichment <- data.table::fread(file.path("/home/boris/Bureau/scShiny/script/scenic/", "int/1.1_grn_output.tsv"), header = F, sep="\t")
grnres<-read.table('/home/boris/Bureau/scShiny/script/scenic/int/1.1_grn_output.tsv',stringsAsFactors=F,sep="\t")
colnames(grnres) <- c("TF", "Target", "weight")
saveRDS(grnres,"int/1.4_GENIE3_linkList.Rds")

#Build and score the GRN (runSCENIC_…)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) #** Only for toy run!!
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered, skipHeatmap = T )
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


#Exploring/interpreting the results
exprMat_log <- exprMat_filtered # Better if it is logged/normalized
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log) # default t-SNE
savedSelections <- shiny::runApp(aucellApp)

#Projection the AUC and TF expression onto t-SNEs
print(tsneFileName(scenicOptions))
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

# Show TF expression:
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat_log, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("JUN")],], plots="Expression")


# Save AUC as PDF:
Cairo::CairoPDF("output/Step4_BinaryRegulonActivity_tSNE_colByAUC.pdf", width=20, height=15)
par(mfrow=c(4,6))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots="AUC")
dev.off()

# Density plot to detect most likely stable states (higher-density areas in the t-SNE):
library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)

#Show several regulons simultaneously:
#par(bg = "black")
par(mfrow=c(1,2))
regulonNames <- c("RPL41","JUN")
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)

regulonNames <- list(red=c("Sox10", "Sox8"),
                     green=c("Irf1"),
                     blue=c( "Tef"))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="Binary")





#GRN: Regulon targets and motifs
regulons <- loadInt(scenicOptions, "regulons")
regulons[c("FOS", "JUN")]

regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))


regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="FOS" & highConfAnnot==TRUE]
viewMotifs(tableSubset, options=list(pageLength=5)) 


motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="JUN"]
viewMotifs(tableSubset) 



#Regulators for known cell types or clusters
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(cellInfo$Phénotype,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")




# regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "Phénotype"])
rssPlot <- plotRSS(rss)


plotRSS_oneSet(rss, setName = "B-cells")


library(Seurat)
dr_coords <- Embeddings(all, reduction="tsne")

tfs <- c("IRF8_extended","SPIB","SPI1", "ATF5_extended")
par(mfrow=c(2,2))
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, tfs), plots = "AUC")

















## -- Sheet list SCENIC ------------------------------------------------------------------------------- ##
### Load data
loomPath <- system.file(package="SCENIC", "/home/boris/Bureau/scShiny/script/scenic/R/FL09.loom")
loom <- open_loom("/home/boris/Bureau/lab/IgCaller/R/FL09.loom")
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

### Initialize settings
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
for(featherURL in dbFiles) {download.file(featherURL, destfile=paste0("/home/boris/Bureau/scShiny/script/scenic/DB/",basename(featherURL)))} # saved in current dir
  
scenicOptions <- initializeScenic(org="hgnc", dbDir="/home/boris/Bureau/scShiny/script/scenic/DB/", nCores=12)
# scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="/home/boris/Bureau/scShiny/script/scenic/DB/scenicOptions.Rds") 

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

# Optional: Binarize activity
# aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
# savedSelections <- shiny::runApp(aucellApp)
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

# Export:
# saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo")) # Temporary, to add to loom
export2loom(scenicOptions, exprMat)

# To save the current status, or any changes in settings, save the object again:
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### Exploring output 
# Check files in folder 'output'
# Browse the output .loom file @ http://scope.aertslab.org

# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Sox8"]
viewMotifs(tableSubset) 

# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 

# Cell-type specific regulators (RSS): 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)