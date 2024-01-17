source(file = "/home/boris/Bureau/Flinovo/fonction/subset.R")
source(file = "/home/boris/Bureau/Flinovo/library.R")
# Load package
library(ParetoTI)
library(ggplot2)
library(SingleCellExperiment)

# If package does not load because "py_pcha is not found" make sure you do not have
# a python environment already loaded in this R session (e.g. restart R and try loading again).
# Install python dependencies (like py_pcha) into python conda environment, and (optionally) install *extra_packages*.
ParetoTI::install_py_pcha(method = "conda", extra_packages = c("tensorflow", "tensorflow-probability","pandas", "keras", "h5py","geosketch", "pydot", "scikit-learn==0.20","umap-learn"))
reticulate::py_discover_config("py_pcha") # Finally, check that py_pcha library is successfully installed and discoverable
reticulate::use_condaenv("reticulate_PCHA", conda = "auto", required = TRUE) # To make sure R uses the correct conda enviroment you can run this when you start R: set TRUE to force R to use reticulate_PCHA

##################################################################################################################################
siege <- c( "FL09C1164", "FL02G095",  "FL05G0330",  "FL12C1888","FL08G0293")#, "FL140304", "FL08G0431")
patient <- "FL06G1206"
parti <- function(singlet, condition){  
  sce <- as.SingleCellExperiment(singlet, assay = 'RNA')
  
  ##################################################################################################################################
  ## -- 1. Load data from GEO and filter as described in the paper, normalise and PCs for finding polytopes -- ##
  # Annotation MT genes
  go_annot = map_go_annot(taxonomy_id = 9606, keys = rownames(sce), columns = c("GOALL"), keytype = "SYMBOL", ontology_type = c("CC")) # look at nuclear-encoded MT genes (find those genes using GO annotations)
  mt_gene = unique(go_annot$annot_dt[GOALL == "GO:0005739", SYMBOL])
  sce$all_mito_genes = colSums(counts(sce[mt_gene, ])) / colSums(counts(sce)) #qplot(sce$perc.mito, sce$all_mito_genes, geom = "bin2d")
  
  # Filtering
  sce = sce[,colSums(counts(sce)) > 1000 & colSums(counts(sce)) < 30000] # remove cells with more less than 1000 or more than 30000 UMI
  sce = sce[rowMeans(counts(sce) > 0) > 0.05,]                           # remove genes with too many zeros (> 95% cells)
  sce = sce[,colMeans(counts(sce) > 0) > 0.15]                           # remove cells with too many zeros (> 85%)
  
  # Normalise gene expression by cell sum factors and log-transform
  sce = scran::computeSumFactors(sce)
  sce = scater::logNormCounts(sce)
  sce = scater::runPCA(sce, ncomponents = 7, exprs_values = "logcounts") # Find principal components
  #scater::plotReducedDim(sce, ncomponents = 3, dimred = "PCA")           # Plot PCA colored by batch
  PCs4arch = t(reducedDim(sce, "PCA"))                                   # extract PCs (centered at 0 with runPCA())
  
  gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; 
  ##################################################################################################################################
  ## -- 2. Examine the polytope with best k & look at known markers of subpopulations -- ##
  # Find archetypes
  arc_ks = k_fit_pch(PCs4arch, ks = 2:8, check_installed = T,bootstrap = T, bootstrap_N = 200, maxiter = 100, bootstrap_type = "s", seed = 2543, 
                     volume_ratio = "none", delta=0, conv_crit = 1e-04, order_type = "align", sample_prop = 0.75) # set to "none" if too slow
  #plot_arc_var(arc_ks, type = "varexpl", point_size = 2, line_size = 1.5) + theme_bw()     # Show variance explained by a polytope with each k (cumulative)
  #plot_arc_var(arc_ks, type = "res_varexpl", point_size = 2, line_size = 1.5) + theme_bw() # Show variance explained by k-vertex model on top of k-1 model (each k separately)
  #plot_arc_var(arc_ks, type = "total_var", point_size = 2, line_size = 1.5) + theme_bw() + ylab("Mean variance in position of vertices") # Show variance in position of vertices obtained using bootstraping - use this to find largest k that has low variance
  
  # fit a polytope with bootstraping of cells to see stability of positions
  if(length(arc_ks$summary$varexpl[which(arc_ks$summary$varexpl < 0.8)]) > 1){k = length(arc_ks$summary$varexpl[which(arc_ks$summary$varexpl < 0.8)])} 
  else{k = 2}
  
  arc = fit_pch_bootstrap(PCs4arch, n = 200, sample_prop = 0.75, seed = 235, noc = k, delta = 0, conv_crit = 1e-04, type = "s")
  #plotly::layout(plot_arc(arc_data = arc, data = PCs4arch, which_dimensions = 1:3, line_size = 3.5, text_size = 40, data_size = 3) , title = "PCA") # You can also check which cells have high entropy of logistic regression predictions when classifying all cells in a tissue into cell types. These could have been misclassified by the method and wrongly assigned to sce, or these could be doublets.
  
  # Find archetypes on all data (allows using archetype weights to describe cells)
  arc_1 = fit_pch(PCs4arch, volume_ratio = "none", maxiter = 500, noc = k, delta = 0, conv_crit = 1e-04)
  singlet <- subset(singlet, cell = colnames(PCs4arch))
  
  gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; 
  ##################################################################################################################################
  ## -- 3.a Find genes and gene sets enriched near vertices -- ## CC
  # Map GO annotations and measure activities
  activ = measure_activity(as.matrix(sce@assays@data$logcounts), # row names are assumed to be gene identifiers c("BP", "MF", "CC")
                           which = "CC", return_as_matrix = F, taxonomy_id = 9606, keytype = "SYMBOL", lower = 20, upper = 1000, activity_method ="AUCell",
                           aucell_options = list(aucMaxRank = nrow(sce) * 0.1, binary = F, nCores = 6, plotStats = FALSE))

  # Merge distances, gene expression and gene set activity into one matrix
  data_attr = merge_arch_dist(arc_data = arc_1, data = PCs4arch, feature_data = as.matrix(logcounts(sce)), colData = activ, dist_metric = c("euclidean", "arch_weights")[1], colData_id = "cells", rank = F) 
  
  # Use Wilcox test to find genes maximally expressed in 10% closest to each vertex
  enriched_genes = find_decreasing_wilcox(data_attr$data, data_attr$arc_col,features = data_attr$features_col,bin_prop = 0.1, method = "BioQC")
  enriched_sets = find_decreasing_wilcox(data_attr$data, data_attr$arc_col,features = data_attr$colData_col,bin_prop = 0.1, method = "BioQC")
  
  # Take a look at top genes and functions for each archetype
  labs <- get_top_decreasing(summary_genes = enriched_genes, summary_sets = enriched_sets,cutoff_genes = 0.01, cutoff_sets = 0.05, cutoff_metric = "wilcoxon_p_val", p.adjust.method = "fdr",order_by = "mean_diff", order_decreasing = T, min_max_diff_cutoff_g = 0.4, min_max_diff_cutoff_f = 0.03)
  activ = cbind(activ, singlet@meta.data)
  save(enriched_genes, enriched_sets, labs, arc, PCs4arch, activ, file=paste0("/home/boris/Bureau/scShiny/datasets/ParTI/",condition,"_",patient,"_parti_CC.RData"))

  gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; 
  ##################################################################################################################################
  ## -- 3.b Find genes and gene sets enriched near vertices -- ## MF
  activ = measure_activity(as.matrix(sce@assays@data$logcounts), which = "MF", return_as_matrix = F, taxonomy_id = 9606, keytype = "SYMBOL", lower = 20, upper = 1000, activity_method ="AUCell",
                           aucell_options = list(aucMaxRank = nrow(sce) * 0.1, binary = F, nCores = 6, plotStats = FALSE))
  data_attr = merge_arch_dist(arc_data = arc_1, data = PCs4arch, feature_data = as.matrix(logcounts(sce)), colData = activ, dist_metric = c("euclidean", "arch_weights")[1], colData_id = "cells", rank = F) 
  enriched_genes = find_decreasing_wilcox(data_attr$data, data_attr$arc_col,features = data_attr$features_col,bin_prop = 0.1, method = "BioQC")
  enriched_sets = find_decreasing_wilcox(data_attr$data, data_attr$arc_col,features = data_attr$colData_col,bin_prop = 0.1, method = "BioQC")
  labs <- get_top_decreasing(summary_genes = enriched_genes, summary_sets = enriched_sets,cutoff_genes = 0.01, cutoff_sets = 0.05, cutoff_metric = "wilcoxon_p_val", p.adjust.method = "fdr",order_by = "mean_diff", order_decreasing = T, min_max_diff_cutoff_g = 0.4, min_max_diff_cutoff_f = 0.03)
  activ = cbind(activ, singlet@meta.data)
  save(enriched_genes, enriched_sets, labs, arc, PCs4arch, activ, file=paste0("/home/boris/Bureau/scShiny/datasets/ParTI/",condition,"_",patient,"_parti_MF.RData"))
  gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; 
  ##################################################################################################################################
  ## -- 3.c Find genes and gene sets enriched near vertices -- ## BP
  activ = measure_activity(as.matrix(sce@assays@data$logcounts),which = "BP", return_as_matrix = F, taxonomy_id = 9606, keytype = "SYMBOL", lower = 20, upper = 1000, activity_method ="AUCell",
                           aucell_options = list(aucMaxRank = nrow(sce) * 0.1, binary = F, nCores = 6, plotStats = FALSE))
  data_attr = merge_arch_dist(arc_data = arc_1, data = PCs4arch, feature_data = as.matrix(logcounts(sce)), colData = activ, dist_metric = c("euclidean", "arch_weights")[1], colData_id = "cells", rank = F) 
  enriched_genes = find_decreasing_wilcox(data_attr$data, data_attr$arc_col,features = data_attr$features_col,bin_prop = 0.1, method = "BioQC")
  enriched_sets = find_decreasing_wilcox(data_attr$data, data_attr$arc_col,features = data_attr$colData_col,bin_prop = 0.1, method = "BioQC")
  labs <- get_top_decreasing(summary_genes = enriched_genes, summary_sets = enriched_sets,cutoff_genes = 0.01, cutoff_sets = 0.05, cutoff_metric = "wilcoxon_p_val", p.adjust.method = "fdr",order_by = "mean_diff", order_decreasing = T, min_max_diff_cutoff_g = 0.4, min_max_diff_cutoff_f = 0.03)
  activ = cbind(activ, singlet@meta.data)
  save(enriched_genes, enriched_sets, labs, arc, PCs4arch, activ, file=paste0("/home/boris/Bureau/scShiny/datasets/ParTI/",condition,"_",patient,"_parti_BP.RData"))
  gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; gc() ; 
}
for(patient in siege){
  load(file = paste0("/home/boris/Bureau/Flinovo/result/analyse_patient/",patient,"/All_All_",patient,".RData")) 
  
  singlet@assays$RNA@data <- singlet@assays$RNA@data[rownames(singlet@assays$RNA@data)[which(!str_detect(rownames(singlet@assays$RNA@data), "IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]"))],]
  singlet@assays$RNA@scale.data <- singlet@assays$RNA@scale.data[rownames(singlet@assays$RNA@scale.data)[which(!str_detect(rownames(singlet@assays$RNA@scale.data), "IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]"))],]
  singlet@assays$RNA@counts <- singlet@assays$RNA@counts[rownames(singlet@assays$RNA@counts)[which(!str_detect(rownames(singlet@assays$RNA@counts), "IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]"))],]
  singlet <- seurat_subset(singlet, "Phénotype", "B cells")
  
  singlet0 <- seurat_subset(singlet, "Condition", "Pré-greffe")
  parti(singlet0, "Pré-greffe")
  rm(singlet0); gc(); gc(); gc(); gc()

  singlet1 <- seurat_subset(singlet, "Condition", "Excipient")
  parti(singlet1, "Excipient")
  rm(singlet1); gc(); gc(); gc(); gc()

  singlet2 <- seurat_subset(singlet, "Condition", "RCHOP")
  parti(singlet2, "RCHOP")
  rm(singlet2); gc(); gc(); gc(); gc()

  singlet3 <- seurat_subset(singlet, "Condition", c("Excipient","RCHOP"))
  parti(singlet3, "Post-greffe")
  rm(singlet3); gc(); gc(); gc(); gc()
  
  singlet4 <- seurat_subset(singlet, "Condition", c("Excipient","Pré-greffe"))
  parti(singlet4, "Pré-greffe-Excipient")
  rm(singlet4); gc(); gc(); gc(); gc()
  
}


##################################################################################################################################
final<- list()

for(condition in c("RCHOP","Excipient")){final[[condition]] <- list()
  for(database in c("BP","CC","MF")){gene <- list(); fonction <- list()
    for(patient in c("FL09C1164", "FL02G095",  "FL05G0330",  "FL08G0293","FL12C1888", "FL140304", "FL08G0431")){
      load(file=paste0("/home/boris/Bureau/scShiny/datasets/ParTI/",condition,"_",patient,"_parti_",database,".RData"))
      
      info <- unique(stringr::str_split(labs$lab$arch_lab,"\n\n"))
      info <- lapply(info,"[", -1)
      info <- lapply(1:length(info), function(x) info[[x]][info[[x]] != ""])
      
      a <- stringr::str_split(lapply(info,"[", 1),"\n")
      a <- lapply(1:length(a), function(x) a[[x]][a[[x]] != ""])
      gene[[patient]] <- a[!is.na(a)]
      
      b <- stringr::str_split(lapply(info,"[", 2),"\n")
      b <- lapply(1:length(b), function(x) b[[x]][b[[x]] != ""])
      fonction[[patient]] <- b[!is.na(b)]
      
      data <- data.frame(names = stringr::str_sub(names(unlist(fonction)),1,-2), fonction = unlist(fonction))
      data <- as.data.frame(table(data))
      data<-data[which(data$Freq != 0),][,-c(3)]
      result <- list()
      for(i in 1:nrow(data)){result[[as.character(data$fonction[i])]] <- c(result[[as.character(data$fonction[i])]],as.character(data$names[i]))}
    }
    final[[condition]][[database]] <- result
  }
}

##################################################################################################################################
p_pca = plot_arc(arc_data = arc, data = PCs4arch, which_dimensions = 1:3, data_lab = activ[[labs$enriched_sets$y_name[1]]], line_size = 3.5, text_size = 40, data_size = 3)
plotly::layout(p_pca, title = paste(labs$enriched_sets$y_name[1], "activity"))

# Project to UMAP coordinates (3D -> 2D)
arc_umap = arch_to_umap(arc_data = arc, data = PCs4arch, which_dimensions = 1:2, method = c("naive", "umap-learn")) # implemented in R and slow, requires python module
plot_arc(arc_data = arc_umap$arc_data, data = arc_umap$data, which_dimensions = 1:2) + theme_bw()

arc_tsne = arch_to_tsne(arc_data = arc, data = PCs4arch, which_dimensions = 1:2) # Project to tSNE coordinates (3D -> 2D, requires Rtsne package)
plot_arc(arc_data = arc_tsne$arc_data, data = arc_tsne$data, which_dimensions = 1:2) + theme_bw()
