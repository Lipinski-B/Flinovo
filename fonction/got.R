## -- GoT -- ##
got <- function(hotspot, colname, patient){
  BC <- read.table(paste0("/home/boris/Documents/lipinskib/flinovo/data/",patient,"/HTO/barcodes.txt"))
  result <- cbind(hotspot[,17],hotspot[,18],hotspot[,19], hotspot[,19])
  rownames(result) <- unique(unlist(str_split(hotspot[,1],";")))
  colnames(result) <- c("WT","MUT", colname,"SUM")
  
  for (i in 1:length(rownames(result))){
    result[i,4] <- sum(as.integer(result[i,1]),as.integer(result[i,2]),as.integer(result[i,3]))
    if(result[i,1]=="0" & result[i,2]!="0"){result[i,3]  <- colnames(result)[2]}
    else if(result[i,1]!="0" & result[i,2]=="0"){ result[i,3]  <- colnames(result)[1]} 
    else if(result[i,1]=="0" & result[i,2]=="0"){ result[i,3]  <-  "AMB"} 
    else {
      WT <- as.double(result[i,1]) / (as.double(result[i,1]) +  as.double(result[i,2]))
      MUT <- as.double(result[i,2]) / (as.double(result[i,1]) +  as.double(result[i,2]))
      genotype <- max(WT,MUT)
      
      if(genotype  < 0.8){result[i,3]  <- "AMB"} 
      else {
        if(genotype == WT){result[i,3]  <- colnames(result)[1]
        } else {result[i,3]  <- colnames(result)[2]}
      }
    }
  }
  final <- list()
  final[["tomerge"]] <- data.frame(tomerge=result[,3])
  final[["sum"]] <- data.frame(sum=result[,4], hotspot=rep(colname, length(result[,4])))
  rownames(final[["sum"]])<-NULL
  a <- data.frame(tomerge=setdiff(BC$V1, rownames(final[["tomerge"]])))
  rownames(a) <- a$tomerge
  a[,1] <- NA
  
  final[["tomerge"]] <- rbind(final[["tomerge"]],a)
  
  return(final)
}
vlnp <- function(hotspot,mutation, patient){
  gt <- got(read.table(paste0("/home/boris/Documents/lipinskib/flinovo/result/", patient, "/GOT/result/",mutation,".out/",mutation,".summTable.txt"), header = T), mutation, patient)
  
  figure1 <- read.table(paste0("/home/boris/Documents/lipinskib/flinovo/result/", patient,"/GOT/figure_GoT/", hotspot,"/", hotspot,"_nb_umi_figure1.txt"), header = T)
  figure2 <- read.table(paste0("/home/boris/Documents/lipinskib/flinovo/result/", patient,"/GOT/figure_GoT/", hotspot,"_hotspot/", hotspot,"_hotspot_nb_umi_figure1.txt"), header = T)
  figure3 <- data.frame(name=rep(hotspot,length(gt[["sum"]][,1])), sum=as.numeric(gt[["sum"]][,1]))
  
  m <- data.frame(
    name=c(rep("10X",length(figure1[,1])),rep("10X Targeted",length(figure2[,1])),rep("GoT",length(figure3[,1]))), 
    sum=c(as.numeric(figure1[,2]),as.numeric(figure2[,2]),as.numeric(figure3$sum)),
    patien=rep(paste0(patient,": ",round((length(figure3$name)/length(gt$tomerge[,1]))*100,2)," % genotyping"),length(c(as.numeric(figure1[,2]),as.numeric(figure2[,2]),as.numeric(figure3$sum)))),
    hotspot=rep(paste0(mutation,": ",round((length(figure3$name)/length(gt$tomerge[,1]))*100,2)," % genotyping"), length(c(as.numeric(figure1[,2]),as.numeric(figure2[,2]),as.numeric(figure3$sum))))
  )
  return(m)
}