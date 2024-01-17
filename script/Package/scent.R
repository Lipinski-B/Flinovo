source(file = "/home/boris/Bureau/Flinovo/script/work_tools.R")
library(SCENT);
data(net13Jun12);
print(dim(net13Jun12.m));
View(head(net13Jun12.m))

data(dataChu)
print(dim(scChu.m));
View(head(scChu.m))

print(summary(factor(phenoChu.v)));

lscChu.m <- log2(scChu.m+1.1);
range(lscChu.m);

lscChu0.m <- log2(scChu.m+1);
range(lscChu0.m);


integ.l <- DoIntegPPI(exp.m = lscChu.m, ppiA.m = net13Jun12.m)
str(integ.l)

sr.o <- CompSRana(integ.l, local = FALSE, mc.cores = 12)

boxplot(srChu.v ~ phenoChu.v, main = "SR potency estimates", xlab = "Cell Type", ylab = "SR")