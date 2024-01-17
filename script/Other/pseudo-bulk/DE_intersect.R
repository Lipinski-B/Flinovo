source(file = "/home/boris/Bureau/scShiny/script/work_tools.R")
patient <- siege[6]
load(file = paste0("/home/boris/Documents/analyse/Patient/", patient,"/singlet_", patient,".RData"))

################################################################################################################################################################################################################################################################################################################################################
## -- DE intersect -- ## liste = read.table("/home/boris/Bureau/dissoc.txt", header = T) ; liste <- liste$x
DE <- list()
for (patient in siege) {
  load(file = paste0("/home/boris/Bureau/scShiny/datasets/Patient/", patient,".RData"))
  
  DE[[patient]]$DE_RE <- singlet@tools$DE_RE[which(singlet@tools$DE_RE$p_val_adj < 0.05),]
}
partiel  <- Reduce(intersect, list(rownames(DE[["FL140304"]]$DE_RE),rownames(DE[["FL12C1888"]]$DE_RE))) ##Partielle  FL08G0293
complete <- Reduce(intersect, list(rownames(DE[["FL02G095"]]$DE_RE),rownames(DE[["FL08G0431"]]$DE_RE),rownames(DE[["FL09C1164"]]$DE_RE))) #rownames(DE[["FL05G0330"]]$DE_RE),
Reduce(intersect, list(partiel,complete))


partiel  <- Reduce(intersect, list(
  rownames(DE[["FL05G0330"]]$DE_RE),
  rownames(DE[["FL08G0293"]]$DE_RE),
  rownames(DE[["FL140304"]]$DE_RE),
  rownames(DE[["FL02G095"]]$DE_RE),
  rownames(DE[["FL12C1888"]]$DE_RE),
  rownames(DE[["FL08G0431"]]$DE_RE),
  rownames(DE[["FL09C1164"]]$DE_RE)))


library("ggVennDiagram")
partiel <- list(
  FL140304 = rownames(DE[["FL140304"]]$DE_RE),
  FL12C1888 = rownames(DE[["FL12C1888"]]$DE_RE),
  FL08G0293 = rownames(DE[["FL08G0293"]]$DE_RE))
ggVennDiagram(partiel, label_alpha = 1, label = "count") + scale_fill_gradient(low="white",high = "white") + theme(legend.position = "none") + ggtitle("DE RCHOP vs Excipient : Réponse partielle") + theme(plot.title = element_text(hjust = 0.5))


complete <- list(
  FL02G095 = rownames(DE[["FL02G095"]]$DE_RE),
  FL08G0431 = rownames(DE[["FL08G0431"]]$DE_RE),
  FL09C1164 = rownames(DE[["FL09C1164"]]$DE_RE))
ggVennDiagram(complete, label_alpha = 1, label = "count") + scale_fill_gradient(low="white",high = "white")+ theme(legend.position = "none") + ggtitle("DE RCHOP vs Excipient : Réponse complète") + theme(plot.title = element_text(hjust = 0.5))



partiel <- list(
  FL140304 = rownames(DE[["FL140304"]]$DE_RE),
  FL12C1888 = rownames(DE[["FL12C1888"]]$DE_RE),
  #FL08G0293 = rownames(DE[["FL08G0293"]]$DE_RE),
  #FL02G095 = rownames(DE[["FL02G095"]]$DE_RE),
  FL08G0431 = rownames(DE[["FL08G0431"]]$DE_RE),
  FL09C1164 = rownames(DE[["FL09C1164"]]$DE_RE)#,
  #FL05G0330 = rownames(DE[["FL05G0330"]]$DE_RE)
)
ggVennDiagram(partiel, label_alpha = 1, label = "count") + scale_fill_gradient(low="white",high = "white") + theme(legend.position = "none") + ggtitle("DE RCHOP vs Excipient : Réponse partielle") + theme(plot.title = element_text(hjust = 0.5))


vector <- list()
options("digits"=3)
for (patient in siege) {
  load(file = paste0("/home/boris/Bureau/scShiny/datasets/", patient,".RData"))
  DE[[patient]]$DE_RE <- singlet@tools$DE_RE[which(singlet@tools$DE_RE$p_val_adj < 0.05),]
  DE[[patient]]$DE_RE <- DE[[patient]]$DE_RE[,-c(1,3,4)]
  
  if(patient %in% c("FL140304","FL12C1888","FL08G0293")){DE[[patient]]$DE_RE$Groupe <- "Partielle" } else { DE[[patient]]$DE_RE$Groupe <- "Complete"}
  
  
  for (gene in partiel) {
    finale <- as.matrix(DE[[patient]]$DE_RE[gene,])
    if (!is.na(finale)) {
      rownames(finale) <- patient
      vector[[gene]] <- rbind(vector[[gene]],finale)
    }
  }
  for (gene in complete) {
    finale <- as.matrix(DE[[patient]]$DE_RE[gene,])
    if (!is.na(finale)) {
      rownames(finale) <- patient
      vector[[gene]] <- rbind(vector[[gene]],finale)
    }
  }
  
  
}






df <- data.frame(p_val_adj = vector[[1]][,1], patient = rep('Partielle',5), gene = rep('GAPDH',5))
ggplot(df, aes(x=gene,y=p_val_adj)) + geom_dotplot(binaxis='y', stackdir='center')

vector[[2]] %>% DT::datatable(rownames = T)

DE_PE_candidat <- setdiff(DE_PE,liste)

DE_PE_FC <- data.frame(matrix(ncol = 40, nrow = 6))
rownames(DE_PE_FC) <- siege ; colnames(DE_PE_FC) <- DE_PE_candidat
for (patient in siege) {DE_PE_FC[patient,] <- DE[[patient]]$DE_PE[rownames(DE[[patient]]$DE_PE) %in% DE_PE_candidat,]$avg_log2FC}

write.table(DE_PE_FC, "/home/boris/Bureau/DE_PE_FC.txt")
write.table(DE_PE_candidat_FC, "/home/boris/Bureau/DE_PE_FC.txt")


DE_PE <- Reduce(intersect, list(rownames(DE[["FL140304"]]$DE_PE),rownames(DE[["FL12C1888"]]$DE_PE)[1:X],rownames(DE[["FL09C1164"]]$DE_PE)[1:X],rownames(DE[["FL08G0293"]]$DE_PE)[1:X],rownames(DE[["FL02G095"]]$DE_PE)[1:X],rownames(DE[["FL05G0330"]]$DE_PE)[1:X]))

