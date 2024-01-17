library(monocle)
library(ggplot2)
library(cowplot)
library(clues)
library(dplyr)
library(reshape2)
library(tidyr)
library(cellrangerRkit)

# monocle workflow
# 1 - Store Data in a CellDataSet Object
patient <- "hFL_130337"
load(file = paste0("/home/boris/Documents/analyse/singlet_",patient,".RData"))

importCDS <- function(otherCDS, seurat_scale=F, import_all = FALSE){
  if (class(otherCDS)[1] == "Seurat") {
    requireNamespace("Seurat")
    if (!seurat_scale) {
      data <- otherCDS@assays$RNA@counts
    } else {
      data <- otherCDS@assays$RNA@scale.data
    }
    if (class(data) == "data.frame") {
      data <- as(as.matrix(data), "sparseMatrix")
    }
    pd <- tryCatch({
      pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
      pd
    }, error = function(e) {
      pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd <- new("AnnotatedDataFrame", data = pData)
      message("This Seurat object doesn't provide any meta data")
      pd
    })
    if (length(setdiff(colnames(data), rownames(pd))) > 0) {
      data <- data[, rownames(pd)]
    }
    fData <- data.frame(gene_short_name = row.names(data), 
                        row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    #lowerDetectionLimit <- otherCDS@is.expr
    if (all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
      expr <- "negbinomial.size"
    }
    else if (any(data < 0)) {
      expressionFamily <- uninormal()
      expr <- "unimormal"
    }
    else {
      expressionFamily <- tobit()
      expr <- "tobit"
    }
    print(paste0("expressionFamily ",expr))
    # valid_data <- data[, row.names(pd)]
    monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
                                  #lowerDetectionLimit = lowerDetectionLimit,
                                  expressionFamily = expressionFamily)
    if (import_all) {
      if ("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
      }
      else {
        mist_list <- otherCDS
      }
    }
    else {
      mist_list <- list()
    }
    if ("var.genes" %in% slotNames(otherCDS)) {
      var.genes <- setOrderingFilter(monocle_cds, otherCDS@var.genes)
    }
    monocle_cds@auxClusteringData$seurat <- mist_list
  }
  else if (class(otherCDS)[1] == "SCESet") {
    requireNamespace("scater")
    message("Converting the exprs data in log scale back to original scale ...")
    data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset
    fd <- otherCDS@featureData
    pd <- otherCDS@phenoData
    experimentData = otherCDS@experimentData
    if ("is.expr" %in% slotNames(otherCDS)) 
      lowerDetectionLimit <- otherCDS@is.expr
    else lowerDetectionLimit <- 1
    if (all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    }
    else if (any(data < 0)) {
      expressionFamily <- uninormal()
    }
    else {
      expressionFamily <- tobit()
    }
    if (import_all) {
      mist_list <- otherCDS
    }
    else {
      mist_list <- list()
    }
    monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
                                  lowerDetectionLimit = lowerDetectionLimit, expressionFamily = expressionFamily)
    monocle_cds@auxOrderingData$scran <- mist_list
  }
  else {
    stop("the object type you want to export to is not supported yet")
  }
  return(monocle_cds)
} # Importer object seurat (sous forme de matrix) dans monocle: CellDataSet (CDS)
data <- importCDS(singlet, import_all = TRUE)
data <- estimateSizeFactors(data)         # Estimate size factors 
data <- estimateDispersions(data)         # and dispersions 


## -- Filtering low-quality cells -- ##
#valid_cells <- row.names(subset(pData(data), Cells.in.Well == 1, Control == FALSE & Clump == FALSE & Debris == FALSE & Mapped.Fragments > 1000000))
#data <- data[,valid_cells]
pData(data)$Total_mRNAs <- Matrix::colSums(exprs(data))
data <- data[,pData(data)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(data)$Total_mRNAs)) + 2*sd(log10(pData(data)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(data)$Total_mRNAs)) - 2*sd(log10(pData(data)$Total_mRNAs)))
qplot(Total_mRNAs, data = pData(data),  geom = "density") + geom_vline(xintercept = lower_bound) + geom_vline(xintercept = upper_bound)

data <- data[,pData(data)$Total_mRNAs > lower_bound & pData(data)$Total_mRNAs < upper_bound]
data <- detectGenes(data, min_expr = 0.1)


# 2 - Classify cells with known marker genes ; ## -- Classifying and Counting Cells -- ##
#by type
#marqueur des B
#CD19_id <- row.names(subset(fData(data), gene_short_name == "CD19"))
#MS4A1_id <- row.names(subset(fData(data), gene_short_name == "MS4A1"))
#CD79A_id <- row.names(subset(fData(data), gene_short_name == "CD79A"))
#CD79B_id <- row.names(subset(fData(data), gene_short_name == "CD79B"))

#cth <- newCellTypeHierarchy()
#cth <- addCellType(cth, "T cells", classify_func = function(x){ x[MS4A1_id,] < 1 & x[CD19_id,] < 1 & x[CD79A_id,] < 1 & x[CD79B_id,] < 1})
#cth <- addCellType(cth, "B cells", classify_func = function(x){ x[MS4A1_id,] > 1 & x[CD19_id,] > 1})# & x[CD79A_id,] > 1 & x[CD79B_id,] > 1})
#data <- classifyCells(data, cth, 0.1)
#table(pData(data)$CellType)

#pie <- ggplot(pData(data), aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
#pie + coord_polar(theta = "y") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())


# 3 - Cluster your cells
#without marker gene
disp_table <- dispersionTable(data)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
data <- setOrderingFilter(data, unsup_clustering_genes$gene_id)
plot_ordering_genes(data)

# data@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(data, return_all = F) # norm_method='log'

data <- reduceDimension(data, max_components = 2, num_dim = 6, reduction_method = 'tSNE', verbose = T)
data <- clusterCells(data, num_clusters = 2)
plot_cell_clusters(data, 1, 2, color = "CellType", markers = c("CD19", "MS4A1"))

data <- reduceDimension(data, max_components = 2, num_dim = 2, reduction_method = 'tSNE', residualModelFormulaStr = "~num_genes_expressed", verbose = T)
data <- clusterCells(data, num_clusters = 2)
plot_cell_clusters(data, 1, 2, color = "CellType")

data <- clusterCells(data, num_clusters = 2)
plot_cell_clusters(data, 1, 2, color = "Cluster") + facet_wrap(~CellType)


#using marker genes
expressed_genes <- row.names(subset(fData(data), num_cells_expressed >= 10))
marker_diff <- markerDiffTable(data[expressed_genes,], cth, residualModelFormulaStr = "~CellType + num_genes_expressed", cores = 1)

candidate_clustering_genes <- row.names(subset(marker_diff, qval < 0.01))
marker_spec <- calculateMarkerSpecificity(data[candidate_clustering_genes,], cth)
head(selectTopMarkers(marker_spec, 3))

semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
data <- setOrderingFilter(data, semisup_clustering_genes)
plot_ordering_genes(data)

plot_pc_variance_explained(data, return_all = F)

data <- reduceDimension(data, max_components = 2, num_dim = 3, norm_method = 'log', reduction_method = 'tSNE', residualModelFormulaStr = "~CellType + num_genes_expressed", verbose = T)
data <- clusterCells(data, num_clusters = 2)
plot_cell_clusters(data, 1, 2, color = "CellType")


#Imputing cell type
data <- clusterCells(data, num_clusters = 2, frequency_thresh = 0.1, cell_type_hierarchy = cth)
plot_cell_clusters(data, 1, 2, color = "CellType", markers = c("MS4A1", "CXCR4", "CD83"))

pie <- ggplot(pData(data), aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())






# 4 - Order cells in pseudotime along a trajectory ; ## -- Constructing Single Cell Trajectories -- ##
# The ordering workflow : a - choosing genes that define progress
head(fData(data))
head(pData(data))

expressed_genes <- row.names(subset(fData(data), num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(data[expressed_genes,], fullModelFormulaStr = "~HTO_maxID") #find all genes that are differentially expressed in response to the switch from the model
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))

data <- setOrderingFilter(data, ordering_genes) # methods to select genes that require no knowledge of the design of the experiment at all
plot_ordering_genes(data)

# b - reducing the dimensionality of the data 
data <- reduceDimension(data, max_components = 2, method = 'DDRTree')

# c - ordering the cells in pseudotime 
data <- orderCells(data)
plot_cell_trajectory(data, color_by = "HTO_maxID")
plot_cell_trajectory(data, color_by = "SingleR.calls")
plot_cell_trajectory(data, color_by = "State") + facet_wrap(~State, nrow = 1)
plot_cell_trajectory(data, color_by = "Pseudotime")
plot_cell_trajectory(data, color_by = "Cluster")




#jitter plot to pick figure out which state corresponds to rapid proliferation
blast_genes <- row.names(subset(fData(data), gene_short_name %in% c("CD19", "MS4A1")))
plot_genes_jitter(data[blast_genes,], grouping = "State", min_expr = 0.1)


#confirm that the ordering is correct 
data_expressed_genes <-  row.names(subset(fData(data), num_cells_expressed >= 10))
data_filtered <- data[data_expressed_genes,]
my_genes <- row.names(subset(fData(data_filtered), gene_short_name %in% c("CD3E","MS4A1", "CD19")))
cds_subset <- data_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "SingleR.calls")


# 5 - Perform differential expression analysis
# 6 - Analyse of branch : The branches occur because cells execute alternative gene expression programs. 
plot_cell_trajectory(data, color_by = "HTO_maxID")
BEAM_res <- BEAM(data, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]


#visualize modules of genes that have similar lineage-dependent expression patterns. 
plot_genes_branched_heatmap(lung[row.names(subset(BEAM_res, qval < 1e-4)),], branch_point = 1, num_clusters = 4, cores = 1, use_gene_short_name = T, show_rownames = T)



lung_genes <- row.names(subset(fData(lung), gene_short_name %in% c("CD3E","MS4A1", "CD19")))
plot_genes_branched_pseudotime(lung[lung_genes,], branch_point = 1, color_by = "HTO_maxID", ncol = 1)

















###########################


### A partir de là c'est plus la partie du blog où il clusterise lui en fontion de cluster
### Je me suis pas vraiment penchée dessus, parce que je sais pas trop si c'est pertinent
############################## Constructing Single Cell Trajectories ##############################
### For the trajectory analysis in this post, I will use only a subset of all cells (cluster 1) and genes that are expressed in at least 10 cells.
expressed_genes <- row.names(subset(fData(data), num_cells_expressed >= 10))

# my_colour is TRUE if a cell belongs to either cluster 3, 4, 7 or 10
data_subset <- data[expressed_genes, pData(data)$my_colour]
data_subset


# keeping only genes expressed in greater than 5% of cells
data_subset <- detectGenes(data_subset, min_expr = 0.1)
fData(data_subset)$use_for_ordering <- fData(data_subset)$num_cells_expressed > 0.05 * ncol(data_subset)

# how many genes are used?
table(fData(data_subset)$use_for_ordering)

plot_pc_variance_explained(data_subset, return_all = FALSE)

data_subset <- reduceDimension(data_subset, max_components = 2,norm_method = 'log',num_dim = 10,reduction_method = 'tSNE',verbose = TRUE)
data_subset <- clusterCells(data_subset, verbose = FALSE)
plot_rho_delta(data_subset, rho_threshold = 2, delta_threshold = 10)

data_subset <- clusterCells(data_subset,rho_threshold = 2,delta_threshold = 10,skip_rho_sigma = T,verbose = FALSE)
table(pData(data_subset)$Cluster)
plot_cell_clusters(data_subset)

clustering_DEG_genes <- differentialGeneTest(data_subset,fullModelFormulaStr = '~Cluster',cores = 8)
dim(clustering_DEG_genes)

clustering_DEG_genes %>% arrange(qval) %>% head()

my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
data_subset <- setOrderingFilter(data_subset, ordering_genes = my_ordering_genes)
data_subset <- reduceDimension(data_subset, method = 'DDRTree')

# the warnings were for use of deprecated code
data_subset <- orderCells(data_subset)

plot_cell_trajectory(data_subset, color_by = "Cluster")

# pseudotime is now a column in the phenotypic data as well as the cell state
head(pData(data_subset))

my_pseudotime_de <- differentialGeneTest(data_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 8)
my_pseudotime_de %>% arrange(qval) %>% head()

# save the top 6 genes
my_pseudotime_de %>% arrange(qval) %>% head(6) %>% select(gene_short_name) -> my_pseudotime_gene
my_pseudotime_gene <- my_pseudotime_gene$gene_short_name

plot_genes_in_pseudotime(data_subset[my_pseudotime_gene,])

# cluster the top 50 genes that vary as a function of pseudotime
my_pseudotime_de %>% arrange(qval) %>% head(50) %>% select(status) -> gene_to_cluster
gene_to_cluster <- gene_to_cluster$status

my_pseudotime_cluster <- plot_pseudotime_heatmap(data_subset[gene_to_cluster,], num_clusters = 3, cores = 8, show_rownames = TRUE, return_heatmap = TRUE)