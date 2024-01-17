library(edgeR)
library(ggplot2)

annotation <- read.csv("/media/boris/Elements/RNAseq_Chloé_Patients SMZL_32_mRNAs/Analyses_Sylvain/AHCYL2_QC.tar/AHCYL2_QC/featureCounts/annotation.csv")
counts <- readRDS("/media/boris/Elements/RNAseq_Chloé_Patients SMZL_32_mRNAs/Analyses_Sylvain/AHCYL2_QC.tar/AHCYL2_QC/featureCounts/all_counts.rds")
dge <- DGEList(counts=counts, genes=annotation)
dge <- calcNormFactors(dge, method="TMM")
cpm <- cpm(dge)


#ENSG00000141510.17 = p53 
#ENSG00000124762.13 = CDKN1A = p21
#ENSG00000158467.16 = AHCYL2


result <- t(rbind(cpm["ENSG00000141510.17",], cpm["ENSG00000124762.13",], cpm["ENSG00000158467.16",]))
colnames(result) <- c('p53','p21','AHCYL2')
View(result)

cor.test(result[,'p21'], result[,'AHCYL2'], method="pearson")
cor.test(result[,'p53'], result[,'AHCYL2'],method="pearson")

write.csv(result, file = "/home/boris/Bureau/correlation.csv")

result <- as.data.frame(result)



ggplot(result, aes(x=AHCYL2, y=p53)) + 
  geom_point()+
  geom_smooth(method=lm, color="black") + theme_classic() + geom_text(x = 28, y = 100, label =  paste("y =", round(m$coefficients[2],3), "x +", round(m$coefficients[1],3), "; r² =", r2)  , parse = F)

lm_eqn <- function(result){
  x <- result$AHCYL2
  y <- result$p53
  m <- lm(y ~ x, result);
  
  r2 = format(summary(m)$r.squared, digits = 3)

  eq <- paste("y =", round(m$coefficients[2],3), "+", round(m$coefficients[1],3), "; r² =", r2)  
}



