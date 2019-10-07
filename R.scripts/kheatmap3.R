library(ComplexHeatmap)
library(circlize)
library(hgu133plus2hsentrezgprobe)
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2probe)
library(hgu133plus2cdf)
library(hgu133plus2.db)
library(limma)
library(sva)
library(annotate)

#experiment labels
KMexp <- read.csv2("/exp.txt", sep="\t", as.is=T)
setwd("/celfiles")

##read in Data
data1 <- ReadAffy(filenames=KMexp$Cel.name)

#Normalise probe level
#expr.mas <- mas5(data1)
expr.gcrma <- gcrma(data1)
expr.rma <- rma(data1)
#expr.mas <- exprs(x.mas)
expr.gcrma <- exprs(expr.gcrma)
colnames(expr.gcrma) <- as.character(KMexp$ID)

##ComBAT##
setwd("/GEP")
variance = apply(expr.gcrma,MARGIN=1,var)
summary(variance<=0)
edata = as.matrix(expr.gcrma[rownames(expr.gcrma[variance>0,]),])
samplenames <- data.frame(colnames(edata))
modcombat = model.matrix(~1,data=samplenames)
batch = as.numeric(KMexp$Batch)
combat_normexprs = ComBat(dat=edata, mod=modcombat, batch=batch, par.prior=FALSE, prior.plots=FALSE)

#Limma test for differences probe 
SibShip <- factor(KMexp$Expt)
Treatment <- factor(KMexp$Treat, levels=c("Control","Treat"))
design <- model.matrix(~0+SibShip+Treatment, data=KMexp)
colnames(design) <- gsub("SibShip","",colnames(design))
colnames(design) <- gsub("Treatment","",colnames(design))
samplenames1 <- as.character(samplenames$colnames.edata.)
row.names(design) <- samplenames1
combat_normexprs <- as.matrix(combat_normexprs)

fit <- lmFit(combat_normexprs, design)
fit1 <- eBayes(fit)
fit1$genes$Symbol <- getSYMBOL(rownames(fit1), "hgu133plus2")
write.table(topTable(fit1,coef="Treat",adjust="fdr",number=nrow(combat_normexprs)),file="fit1table_allcells_combat.txt",sep="\t",quote=FALSE)

CellLine <- factor(KMexp$CellLine)
design2 <- model.matrix(~0+CellLine+Treatment, data=KMexp)
colnames(design2) <- gsub("CellLine","",colnames(design2))
colnames(design2) <- gsub("Treatment","",colnames(design2))
row.names(design2) <- samplenames1

fit0 <- lmFit(combat_normexprs, design2)
fit2 <- eBayes(fit0)
fit2$genes$Symbol <- getSYMBOL(rownames(fit2), "hgu133plus2")
write.table(topTable(fit2,coef="Treat",adjust="fdr",number=nrow(combat_normexprs)),file="fit2table_allcells_combat.txt",sep="\t",quote=FALSE)

fit3 <- lmFit(combat_normexprsRMA, design2)
fit4 <- eBayes(fit3)
fit4$genes$Symbol <- getSYMBOL(rownames(fit2), "hgu133plus2")
write.table(topTable(fit4,coef="Treat",adjust="fdr",number=nrow(combat_normexprsRMA)),file="fit4table_allcells_combat.txt",sep="\t",quote=FALSE)

############
##VolcPlot##
############

library(dplyr)
library(ggplot2)
library(ggrepel)

Resultsfit1mut = mutate(Resultsfit1, sig=ifelse(Resultsfit1$adj.P.Val<0.00005, "FDR<0.00005", "Not Sig"))
p = ggplot(Resultsfit1mut, aes(logFC, -log10(P.Value))) +
  geom_point(aes(col=sig)) +
  scale_color_manual(values=c("red", "black"))
p
p+theme_bw()+geom_text_repel(data=filter(Resultsfit1mut, adj.P.Val<0.00005), aes(label=Symbol))


Resultsfit2mut = mutate(Resultsfit2, sig=ifelse(Resultsfit2$adj.P.Val<0.0000005, "FDR<0.0000005", "Not Sig"))
p = ggplot(Resultsfit2mut, aes(logFC, -log10(P.Value))) +
  geom_point(aes(col=sig)) +
  scale_color_manual(values=c("red", "black"))
p
p+theme_bw()+geom_text_repel(data=filter(Resultsfit2mut, adj.P.Val<0.0000005), aes(label=Symbol))

###########
##HeatMap##
###########
UPDOWNK <- read.table ("G:/K GEP/UPDOWNfinal_fit2table_Ks_combat.txt",sep="\t", header = TRUE, row.names = 1)
UPDOWNKProbes <- rownames(UPDOWNK)
silmKMupdown <- expr.gcrma[row.names(expr.gcrma) %in% UPDOWNKProbes,]
UPDOWNKO <- UPDOWNK[order(rownames(UPDOWNK)),]
row.names(silmKMupdown) <- as.character(UPDOWNKO$Symbol)
annK <- KMexp[c(2,3,4,5,6)]
annK <- data.frame(annK[,-1], row.names=annK[,1])
ka1 = HeatmapAnnotation(df=annK, col = list(CellLine = c("KMS11" = "ORANGE","RPMI" = "MAGENTA", "L363" ="BLUE"), Rep = c("rep1" = "BLACK", "rep2"= "DARKGREY", rep3="gray90"), Treat =c("Treat"="RED", "Control"= "GREEN")))

Heatmap(silmKMupdown, name ="log2Exp", top_annotation = ka1, row_names_gp = gpar(fontsize = 2), column_names_gp = gpar(fontsize = 4), clustering_distance_columns = "pearson", top_annotation_height = unit(1.4, "cm"))

###########################
####GENE LEVEL ANALYSIS####
###########################
data2 <- data1
#change probes designations to entrez
data2@cdfName <- "hgu133plus2hsentrezg"

#Normalise gene level
#gene.mas <- mas5(data1)
gene.gcrma <- gcrma(data2)
#exprG.mas <- exprs(gene.mas)
exprG.gcrma <- exprs(gene.gcrma)
colnames(exprG.gcrma) <- as.character(KMexp$ID)

variance1 = apply(exprG.gcrma,MARGIN=1,var)
summary(variance1<=0)
edata1 = as.matrix(exprG.gcrma[rownames(exprG.gcrma[variance1>0,]),])
combat_normexprsG = ComBat(dat=edata1, mod=modcombat, batch=batch, par.prior=FALSE, prior.plots=FALSE)





#Limma test for differences GeneLeveL 
SibShip <- factor(KMexp$Expt)
Treatment <- factor(KMexp$Treat, levels=c("Control","Treat"))
design <- model.matrix(~SibShip+Treatment)
fit <- lmFit(exprG.maslog2, design)
fit <- eBayes(fit)
topTable(fit, coef="TreatT")
topTable(fit, coef="TreatmentTreat")

#basic 
resutlsTable <- topTable(fit, coef="TreatmentTreat", number=20414)

#adjust 
resutlsTable2 <- topTable(fit, coef="TreatmentTreat", adjust="BH", number=20414)
TopresultsTable2 <- subset (resutlsTable2, resutlsTable2$adj.P.Val >=0.05)
TopresultsTable3 <- subset (resutlsTable2, resutlsTable2$adj.P.Val <=0.005)


#top 5000 variable genes 
variance = apply(exprG.mas,MARGIN=1,var)
hist(variance,breaks=10)$breaks[2]
summary(variance<=0.05)
slimG=exprG.maslog2[rev(order(variance))[1:500],]

#Or
#mads=apply(exprG.maslog2,MARGIN=1,mad) #median absolute deviation
#hist(slimG,breaks=10)$breaks[2]
#slimG=exprG.maslog2[rev(order(mads))[1:500],]

#labels for heatmap 
names(exp)
exp1 <- exp[c(3,4,5,6)]
colnames(exp1) <- c("CellLine","Expt","Rep","Treatment")
ha1 = HeatmapAnnotation(df=exp1)

#heatmap of expts with no contrasts
Heatmap(slim, name ="Mas5.gExp", top_annotation = ha1, row_names_gp = gpar(fontsize = 1), column_names_gp = gpar(fontsize = 1), clustering_distance_columns = "pearson")





#List of differing probes
diffprobes <- as.character(rownames(TopresultsTable3))

#Subset expression with differing probes 
exprG.maslog2Diff <- as.matrix(exprG.maslog2[rownames(exprG.maslog2) %in% diffprobes,])







#heatmap of expts with no contrasts
Heatmap(exprG.maslog2Diff, name ="Mas5.gExp", top_annotation = ka1, row_names_gp = gpar(fontsize =2), column_names_gp = gpar(fontsize = 6), clustering_distance_columns = "pearson")

