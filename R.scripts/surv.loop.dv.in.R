library(survival)
library(reshape2)
library(plyr)

setwd("/mnt/scratch/DMP/MYEGRP/johnson/meth.os/Jan17.os")
data.os <- read.table("/mnt/scratch/DMP/MYEGRP/johnson/meth.os/Jan17.os/data.os.jan.in.csv",sep=",", header = TRUE, row.names = 1)
osdemo <- read.table("/mnt/scratch/DMP/MYEGRP/johnson/meth.os/Jan17.os/Meth.array.89_PFS_OS_161201.txt", sep="\t", header = TRUE, row.names = 1)
kmdata <- read.table("/mnt/scratch/DMP/MYEGRP/johnson/meth.os/Jan17.os/kmdata.jan.csv",sep=",", header = TRUE, row.names = 1)
annt3 <- read.table ("/mnt/scratch/DMP/MYEGRP/johnson/meth.os/Jan17.os/hm450.annot.GRCh38.chromHMM.txt",sep="\t", header = TRUE, row.names = 1)
colnames(data.os) = sub("X","",colnames(data.os))
methsamps <- as.character(colnames(data.os))
osdata <- osdemo[rownames(osdemo) %in% methsamps,]
tdata <- as.data.frame(t(data.os))
rnames <- row.names(data.os)
tdata <- tdata[, !sapply(tdata, function(x) { sd(x) == 0} )]

perl.data.in <- merge(osdata, tdata, by=0, all=TRUE)
perl.data.in <- data.frame(perl.data.in[,-1], row.names=perl.data.in[,1])

df.m <- melt(perl.data.in, id.vars=names(perl.data.in)[1:11],variable.name="gene")

GeneSet_chisq <- ddply(df.m, .(gene), function(x){
  a <- survdiff(Surv(PFS.months, PFS.status=='1') ~ value,x)
  a$chisq})

GeneSet_n <- ddply(df.m, .(gene), function(x){
  a <- survdiff(Surv(PFS.months, PFS.status=='1') ~ value,x)
  a$n})

GeneSet_merge <-merge(GeneSet_chisq,GeneSet_n, by="gene")
GeneSet_merge$p <- (1 - pchisq(GeneSet_merge$V1, 1))

rownames(GeneSet_merge) <- GeneSet_merge$gene

result1 <- merge(GeneSet_merge, kmdata, by=0, all=TRUE)
result1 <- data.frame(result1[,-1], row.names=result1[,1])

result2 <- merge(result1, annt3, by=0, all=TRUE)
result2 <- data.frame(result2[,-1], row.names=result2[,1])
write.table(result2, file="all.pfs.nocovars.feb.in.GRCh38.csv", sep=",")