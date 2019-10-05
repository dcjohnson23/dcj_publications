library(minfi)
library(sva)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(FlowSorted.Blood.450k)
library(RColorBrewer)
library(limma)
library(Gviz)

setwd("/mnt/scratch/DMP/MYEGRP/johnson/Meth.450/500_Exomes_Project/Meth.450.CL/")

## -------read in epic files - place sheet in data.dir-----------------------------------------------------------------##
data.dir <- "/mnt/scratch/DMP/MYEGRP/johnson/Meth.450/500_Exomes_Project/epic.IDATs"
epic.CL <- read.metharray.sheet(data.dir, pattern="MM_Cell_lines.epic.csv")
rownames(epic.CL) <- epic.CL$ID
RGset.epic.CL <- read.metharray.exp(targets=epic.CL)
RGset.epic.CL@annotation <- c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b2.hg19")

## --------read in 450K files----------------------------------------------------------------##
data.dir2 <- "/mnt/scratch/DMP/MYEGRP/johnson/Meth.450/500_Exomes_Project/450k.IDATs"
m450k.CL <- read.metharray.sheet(data.dir2, pattern="MM_Cell_lines.450.csv")
rownames(m450k.CL) <- m450k.CL$ID
RGset.m450k.CL <- read.metharray.exp(targets=m450k.CL)
RGset.m450k.CL@annotation <- c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19")

## ------------------QC ------------------------------------------------------##
qcReport(RGset.epic.CL, sampNames=epic.CL$Sample_Name, sampGroups=epic.CL$Cell, pdf="qcReport.CL.epic.cell.pdf")
qcReport(RGset.epic.CL, sampNames=epic.CL$Sample_Name, sampGroups=epic.CL$Tnx, pdf="qcReport.CL.epic.tnx.pdf")
qcReport(RGset.m450k.CL, sampNames=m450k.CL$Sample_Name, sampGroups=m450k.CL$Cell, pdf="qcReport.CL.450.cell.pdf")
qcReport(RGset.m450k.CL, sampNames=m450k.CL$Sample_Name, sampGroups=m450k.CL$Origin, pdf="qcReport.CL.450.origin.pdf")
qcReport(RGset.m450k.CL, sampNames=m450k.CL$Sample_Name, sampGroups=m450k.CL$Tnx, pdf="qcReport.CL.450.tnx.pdf")

###############################################################################################
## ------------------preprocessRaw.epic------------------------------------------------------##
###############################################################################################

MSetRaw.epic <- preprocessRaw(RGset.epic.CL) 

## ------------------Remove poor EPIC QC arrays------------------------------------------------------##
epic.qc <- getQC(MSetRaw.epic)
epic.detP <- detectionP (RGset.epic.CL)
epic.qc$Mult <- epic.qc$mMed*epic.qc$uMed
epic.qc$detP <- colMeans(epic.detP)
write.table(epic.qc, file="epic.qc.all.csv", sep=",")
goodQC.epic <- rownames(epic.qc[ which(epic.qc$Mult > 109.5),])
RGset.epic.CL<- RGset.epic.CL[,goodQC.epic]
epic.CL <- epic.CL[goodQC.epic,]

## ------------------Remove poor EPIC background arrays------------------------------------------------------##
keep.epic <- colMeans(epic.detP) < 0.05
RGset.epic.CL <- RGset.epic.CL[,keep.epic]
epic.CL <- epic.CL[keep.epic,]
sampleNames(RGset.epic.CL) <- epic.CL$ID
MSetRaw.epic <- preprocessRaw(RGset.epic.CL) 
epic.qc.post <- getQC(MSetRaw.epic)
write.table(epic.qc.post, file="epic.qc.post.qc.csv", sep=",")
RSetRaw.epic <- ratioConvert(MSetRaw.epic, what = "both", keepCN = TRUE)
GRsetRaw.epic <- mapToGenome(RSetRaw.epic)

## ------------------sex and snps epic------------------------------------------------------##
pdf("plotsex.epic.CL.pdf")
plotSex(getSex(GRsetRaw.epic, cutoff = -2))
dev.off()
predSex <- getSex(GRsetRaw.epic, cutoff = -2)$predictedSex
write.table(predSex, file="sex.epic.CL.csv", sep=",")
snps1 <- getSnpInfo(GRsetRaw.epic)
write.table(snps1, file="snps.epic.CL.csv", sep=",")
snpsBeta1 <- getSnpBeta(RGset.epic.CL)
write.table(snpsBeta1, file="snpsBeta.epic.CL.csv", sep=",")

###############################################################################################
## ------------------preprocessRaw.450k------------------------------------------------------##
###############################################################################################

MSetRaw.450k <- preprocessRaw(RGset.m450k.CL) 

## ------------------Remove poor QC arrays------------------------------------------------------##
m450k.qc <- getQC(MSetRaw.450k)
m450k.detP <- detectionP(RGset.m450k.CL)
m450k.qc$Mult <- m450k.qc$mMed*m450k.qc$uMed
m450k.qc$detP <- colMeans(m450k.detP)
write.table(m450k.qc, file="m450k.qc.all.csv", sep=",")
goodQC.450k <- rownames(m450k.qc[ which(m450k.qc$Mult > 109.5),])
RGset.m450k.CL <- RGset.m450k.CL[,goodQC.450k]
m450k.CL <- m450k.CL[goodQC.450k,]

## ------------------Remove poor QC arrays------------------------------------------------------##
keep.m450k <- colMeans(m450k.detP) < 0.05
RGset.m450k.CL <- RGset.m450k.CL[,keep.m450k]
m450k.CL <- m450k.CL[keep.m450k,]
MSetRaw.450k <- preprocessRaw(RGset.m450k.CL) 
m450k.qc.post <- getQC(MSetRaw.450k)
write.table(m450k.qc.post, file="m450k.qc.post.qc.csv", sep=",")
RSetRaw.450k <- ratioConvert(MSetRaw.450k, what = "both", keepCN = TRUE)
GRSetRaw.450k <- mapToGenome(RSetRaw.450k)

## ------------------sex and snps 450k------------------------------------------------------##
pdf("plotsex.450.CL.pdf")
plotSex(getSex(GRSetRaw.450k, cutoff = -2))
dev.off()
predSex <- getSex(GRSetRaw.450k, cutoff = -2)$predictedSex
write.table(predSex, file="sex.m450.CL.csv", sep=",")
snps2 <- getSnpInfo(GRSetRaw.450k)
write.table(snps2, file="snps.m450.CL.csv", sep=",")
snpsBeta2 <- getSnpBeta(RGset.m450k.CL)
write.table(snpsBeta2, file="snpsBeta.m450.CL.csv", sep=",")

## ---------------------------------------QC ------------------------------------------------------##
pdf("good.epic.qc.CL.pdf")
plotQC(epic.qc)
dev.off()

pdf("good.m450k.CL.pdf")
plotQC(m450k.qc)
dev.off()

pdf("detP.qc.epic.CL.pdf")
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(epic.detP), col=pal[factor(epic.CL$Origin)], las=2,
cex.names=0.2, ylim = c(0,0.06),ylab="Mean detection p-values")
abline(h=0.01,col="red")
legend("topleft", legend=levels(factor(epic.CL$Origin)), fill=pal,
bg="white")
barplot(colMeans(epic.detP), col=pal[factor(epic.CL$Origin)], las=2,
cex.names=0.2, ylim = c(0,0.004), ylab="Mean detection p-values")
legend("topleft", legend=levels(factor(epic.CL$Origin)), fill=pal,
bg="white")
dev.off()

pdf("detP.qc.m450k.CL.pdf")
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(m450k.detP), col=pal[factor(m450k.CL$Origin)], las=2,
cex.names=0.2, ylim = c(0,0.06),ylab="Mean detection p-values")
abline(h=0.01,col="red")
legend("topleft", legend=levels(factor(m450k.CL$Origin)), fill=pal,
bg="white")
barplot(colMeans(m450k.detP), col=pal[factor(m450k.CL$Origin)], las=2,
cex.names=0.2, ylim = c(0,0.004), ylab="Mean detection p-values")
legend("topleft", legend=levels(factor(m450k.CL$Origin)), fill=pal,
bg="white")
dev.off()

###############################################################################
####-----------------------------Combine-----------------------################
###############################################################################

## ------------------combine RGsets------------------------------------------------------##
RGset.epic.m450k.CL <- combineArrays(RGset.m450k.CL, RGset.epic.CL, outType = "IlluminaHumanMethylation450k")
epic.m450k.CL <- pData(RGset.epic.m450k.CL)
MSetRaw.epic.450k <- preprocessRaw(RGset.epic.m450k.CL) 
RSetRaw.epic.450k <- ratioConvert(MSetRaw.epic.450k, what = "both", keepCN = TRUE)
GRSetRawepic450k <- mapToGenome(RSetRaw.epic.450k)



############################################################################################################
#####--------------------------------------Normalise---------------------------------------#################
############################################################################################################

##------------------------normalisation-----------------------------------------------------------------##
MSet.epic.450k.ssnoob <- preprocessNoob(RGset.epic.m450k.CL)
RSetNorm.epic.450k <- ratioConvert(MSet.epic.450k.ssnoob, what = "both", keepCN = TRUE)
GRSetNorm.epic.450k <- mapToGenome(RSetNorm.epic.450k)
annotation.norm <- getAnnotation(GRSetNorm.epic.450k)
##------------------------visualise what the data looks like before and after normalisation---------------##
pdf("RawvsNorm.epic.m450k.pdf")
par(mfrow=c(1,2))
densityPlot(getBeta(MSetRaw.epic.450k), sampGroups=epic.m450k.CL$Origin, main="Raw", legend=FALSE)
legend("top", legend = levels(factor(epic.m450k.CL$Origin)),
 text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(MSet.epic.450k.ssnoob), sampGroups=epic.m450k.CL$Origin,
 main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(epic.m450k.CL$Origin)),
 text.col=brewer.pal(8,"Dark2"))
dev.off()

pdf("RawvsNorm.epic.m450k.Tnx.pdf")
par(mfrow=c(1,2))
densityPlot(getBeta(MSetRaw.epic.450k), sampGroups=epic.m450k.CL$Tnx, main="Raw", legend=FALSE)
legend("top", legend = levels(factor(epic.m450k.CL$Tnx)),
 text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(MSet.epic.450k.ssnoob), sampGroups=epic.m450k.CL$Tnx,
 main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(epic.m450k.CL$Tnx)),
 text.col=brewer.pal(8,"Dark2"))
dev.off()

pdf(file="Norm.combined.bean.pdf", width=8, height=60)
par(mfrow=c(1,1), oma=c(2,10,1,1))
densityBeanPlot(MSet.epic.450k.ssnoob, sampGroups = epic.m450k.CL$Origin, sampNames=epic.m450k.CL$Sample_Name)
legend("topright", legend = levels(factor(epic.m450k.CL$Origin)),
text.col=brewer.pal(8,"Dark2"))
dev.off()

####################################
####----Filter-----#################
####################################

detP <- detectionP(RGset.epic.m450k.CL)
detP <- detP[match(featureNames(GRSetNorm.epic.450k),rownames(detP)),]

##--------filter poor probes - probes that have failed in one or more samples-------------------------##
keep.good <- rowSums(detP < 0.01) == ncol(GRSetNorm.epic.450k)
table(keep.good)
tab.good <- table(keep.good)
write.table(tab.good, file="keep.good.combined.csv", sep=",")
GRSetNorm.epic.450k <- GRSetNorm.epic.450k[keep.good,]

##-------------------filter probes on sex chromosomes-----------------------------------------##
keep.sex <- !(featureNames(GRSetNorm.epic.450k) %in% annotation.norm$Name[annotation.norm$chr %in% c("chrX","chrY")])
tab.sex <- table(keep.sex)
write.table(tab.sex, file="keep.sex.combined.csv", sep=",")
GRSetNorm.epic.450k <- GRSetNorm.epic.450k[keep.sex,]

##-------------------filter probes with SNPs---------------------------------------------##
GRSetNorm.epic.450k <- dropLociWithSnps(GRSetNorm.epic.450k, snps=c("SBE","CpG"), maf=0)

##--------filter cross reactive probes-------------------------##
xReactiveProbes <- read.csv("/mnt/scratch/DMP/MYEGRP/johnson/Meth.450/annotation/48639-non-specific-probes-Illumina450k.csv", sep=",", stringsAsFactors=FALSE)
keep.cross <- !(featureNames(GRSetNorm.epic.450k) %in% xReactiveProbes$TargetID)
tab.cross <- table(keep.cross)
write.table(tab.cross, file="keep.cross.combined.csv", sep=",")
GRSetNorm.epic.450k <- GRSetNorm.epic.450k[keep.cross,]

#########################################
##-------MDS plots---------------########
#########################################

pdf("MDSNorm.epic.m450.Origin.pdf")
par(mfrow=c(1,2))
plotMDS(getM(MSet.epic.450k.ssnoob), top=1000, labels=epic.m450k.CL$Sample_Name, cex = 0.6, gene.selection="common",
col=pal[factor(epic.m450k.CL$Origin)])
legend("top", legend=levels(factor(epic.m450k.CL$Origin)), text.col=pal,
bg="white", cex=0.2)
plotMDS(getM(MSet.epic.450k.ssnoob), top=1000, labels=epic.m450k.CL$Sample_Name, cex = 0.6, gene.selection="common",
col=pal[factor(epic.m450k.CL$Origin)])
legend("top", legend=levels(factor(epic.m450k.CL$Origin)), text.col=pal,
bg="white", cex=0.2)
dev.off()

pdf("MDSNorm.higher.epic.m450.Origin.pdf")
par(mfrow=c(1,3))
plotMDS(getM(MSet.epic.450k.ssnoob), top=1000, labels=epic.m450k.CL$Sample_Name, cex = 0.6, gene.selection="common",
col=pal[factor(epic.m450k.CL$Origin)], dim=c(1,3))
legend("top", legend=levels(factor(epic.m450k.CL$Origin)), text.col=pal,
cex=0.2, bg="white")
plotMDS(getM(MSet.epic.450k.ssnoob), top=1000, labels=epic.m450k.CL$Sample_Name, cex = 0.6, gene.selection="common",
col=pal[factor(epic.m450k.CL$Origin)], dim=c(2,3))
legend("topleft", legend=levels(factor(epic.m450k.CL$Origin)), text.col=pal,
cex=0.2, bg="white")
plotMDS(getM(MSet.epic.450k.ssnoob), top=1000, labels=epic.m450k.CL$Sample_Name, cex = 0.6, gene.selection="common",
col=pal[factor(epic.m450k.CL$Origin)], dim=c(3,4))
legend("topright", legend=levels(factor(epic.m450k.CL$Origin)), text.col=pal,
cex=0.2, bg="white")
dev.off()

pdf("MDS.Norm.epic.m450.Tnx.pdf")
par(mfrow=c(1,2))
plotMDS(getM(MSet.epic.450k.ssnoob), top=1000, labels=epic.m450k.CL$Sample_Name, cex = 0.6, gene.selection="common",
col=pal[factor(epic.m450k.CL$Tnx)])
legend("top", legend=levels(factor(epic.m450k.CL$Tnx)), text.col=pal,
bg="white", cex=0.2)
plotMDS(getM(MSet.epic.450k.ssnoob), top=1000, labels=epic.m450k.CL$Sample_Name, cex = 0.6, gene.selection="common",
col=pal[factor(epic.m450k.CL$Tnx)])
legend("top", legend=levels(factor(epic.m450k.CL$Tnx)), text.col=pal,
bg="white", cex=0.2)
dev.off()

pdf("MDS.higher.Norm.epic.m450.Tnx.pdf")
par(mfrow=c(1,3))
plotMDS(getM(MSet.epic.450k.ssnoob), top=1000, labels=epic.m450k.CL$Sample_Name, cex = 0.6, gene.selection="common",
col=pal[factor(epic.m450k.CL$Tnx)], dim=c(1,3))
legend("top", legend=levels(factor(epic.m450k.CL$Tnx)), text.col=pal,
cex=0.2, bg="white")
plotMDS(getM(MSet.epic.450k.ssnoob), top=1000, labels=epic.m450k.CL$Sample_Name, cex = 0.6, gene.selection="common",
col=pal[factor(epic.m450k.CL$Tnx)], dim=c(2,3))
legend("topleft", legend=levels(factor(epic.m450k.CL$Tnx)), text.col=pal,
cex=0.2, bg="white")
plotMDS(getM(MSet.epic.450k.ssnoob), top=1000, labels=epic.m450k.CL$Sample_Name, cex = 0.6, gene.selection="common",
col=pal[factor(epic.m450k.CL$Tnx)], dim=c(3,4))
legend("topright", legend=levels(factor(epic.m450k.CL$Tnx)), text.col=pal,
cex=0.2, bg="white")
dev.off()

#########################################
##----Outputs------------------##########
#########################################
cellCounts <- estimateCellCounts(RGset.epic.m450k.CL)
rownames(cellCounts) <- epic.m450k.CL$Sample_Name
write.table(cellCounts, file="Pts.CL.cellCounts.csv", sep=",")

##-----write Betas------------------------------------------------------------------------------------##
beta.norm <- getBeta(GRSetNorm.epic.450k)
M.norm <- getM(GRSetNorm.epic.450k)
CN.norm <- getCN(GRSetNorm.epic.450k)
sampleNames.norm <- sampleNames(GRSetNorm.epic.450k)
probeNames.norm <- featureNames(GRSetNorm.epic.450k)
pheno.norm <- pData(GRSetNorm.epic.450k)
gr.norm <- granges(GRSetNorm.epic.450k)
annotation.norm <- getAnnotation(GRSetNorm.epic.450k)
write.table(beta.norm, file="betaNorm.csv", sep=",")
write.table(M.norm, file="MValNorm.csv", sep=",")
## --------------Save Image----------------------------------------------------------##
save.image(file="/mnt/scratch/DMP/MYEGRP/johnson/Meth.450/500_Exomes_Project/Meth.450.CL/Cell.line.pipeline.may.2017.Rdata")