library(affy)
library(simpleaffy)
library(AffyRNADegradation)
library(affyPLM)
library(gcrma)
library(RColorBrewer)
library(gdata)
library(panp)
library(made4)
library(hgu133plus2cdf)
library(sva)

setwd("/GEP") 

##read in Data

exp <- read.csv2("/exp.txt", sep="\t", as.is=T)
View(exp)
cel1 <- exp$Cel.name
data1 <- ReadAffy(filenames=cel1)

##RawQC
myPLM1 <- fitPLM(data1)
deg1 <- AffyRNAdeg(data1)
QC1 <- qc(data1) ##Uses probes used in mas5 norm

pdf("exp.Jul.pdf", width=40, height=25)
par(mar=c(20,4,4,2))
boxplot(myPLM1, main="NUSE - normalisesd unscaled standard errors", outline=F, col="lightblue", las=3, whisklty=0, staplelty=0, ylim=c(0.8,1.2))
Mbox(myPLM1, main="RLE - Relative Log Expression", outline=F, col="mistyrose", las=3, whisklty=0, staplelty=0, ylim=c(-1,1))
dev.off()

##Normalisation QC
x.rma <- call.exprs(data1,"rma")
x.mas <- call.exprs(data1,"mas5")

pdf("overview_mas5.pdf", width=40, height=25)
par(mar=c(20,4,4,2))
overview(exprs(x.mas5)) # overview of the data
dev.off()

pdf("overview_mas5.pdf", width=40, height=25)
par(mar=c(20,4,4,2))
overview(exprs(x.mas)) # overview of the data
dev.off()


##PCA plots 
pc = prcomp(t(exprs(x.mas)))
pc1 = prcomp(t(exprs(x.rma)))

family <- as.factor(exp[,3])
plot( pc$x[ , 1:2 ], pch=19, col=family, main="PCA plot-mas5")
legend("topright", col=1:3, legend=c("KMS11", "L363", "RPMI"), pch = 20, bty='n')

plot( pc1$x[ , 1:2 ], pch=19, col=family, main="PCA plot-rma")
legend("topright", col=1:3, legend=c("KMS11", "L363", "RPMI"), pch = 20, bty='n')


