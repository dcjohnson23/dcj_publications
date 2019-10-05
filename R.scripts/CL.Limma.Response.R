library(limma)

##Read,subset and order samples##
m450k.CL <- read.table("C:/Users/johnson/Desktop/Meth.analysis/Results - Cell Lines/MM_Cell_lines.M450.csv", sep=",", header=T, na.strings="NA", row.names = 1)
epic.CL <- read.table("C:/Users/johnson/Desktop/Meth.analysis/Results - Cell Lines/MM_Cell_lines.EPIC.csv", sep=",", header=T, na.strings="NA", row.names = 1)
total.CL <- rbind(epic.CL, m450k.CL)
total.CL <- total.CL[order(rownames(total.CL)),]

totalMM.CL <- subset(total.CL, total.CL$Expt=="standard")
totalMM.CL <- totalMM.CL[order(rownames(totalMM.CL)),]
totalMM.CL$IMID <- factor(totalMM.CL$IMID)
totalMM.CL$Lab_ID <- factor(totalMM.CL$Lab_ID)
totalMM.CL$Replicate <- factor(totalMM.CL$Replicate)
MM.CL <- rownames(totalMM.CL)

totalMM.CL.Min <- subset(totalMM.CL, totalMM.CL$Replicate==1)
totalMM.CL.Min <- totalMM.CL.Min[order(rownames(totalMM.CL.Min)),]
totalMM.CL.Min$IMID <- factor(totalMM.CL.Min$IMID)
totalMM.CL.Min$Lab_ID <- factor(totalMM.CL.Min$Lab_ID)
MM.CL.Min <- rownames(totalMM.CL.Min)

MVal.cL.ncorr <- read.table("MValNorm.CL.csv", sep=",", header = TRUE, row.names = 1)
colnames(MVal.cL.ncorr) = sub("X","",colnames(MVal.cL.ncorr))
MVal.cL.ncorrMM <- MVal.cL.ncorr[,colnames(MVal.cL.ncorr) %in% MM.CL]
MVal.cL.ncorrMM <- MVal.cL.ncorrMM[,order(colnames(MVal.cL.ncorrMM))]

MVal.cL.ncorrMM.Min <- MVal.cL.ncorrMM[,colnames(MVal.cL.ncorrMM) %in% MM.CL.Min]
MVal.cL.ncorrMM.Min <- MVal.cL.ncorrMM.Min[,order(colnames(MVal.cL.ncorrMM.Min))]

all(diff(match((colnames(MVal.cL.ncorrMM.Min)),(row.names(totalMM.CL.Min)))) > 0)
#TRUE

TREAT <- factor(totalMM.CL.Min$IMID, levels=c("NONRES","SEN"))
design <- model.matrix(~0+TREAT, data=totalMM.CL.Min)
colnames(design) <- gsub("TREAT","",colnames(design))
fit <- lmFit(MVal.cL.ncorrMM.Min, design)
#summary(decideTests(fit))
#   NONRES    SEN
#-1 190480 173398
#0  116771 117950
#1   89793 105696
contMatrix <- makeContrasts(NONRES-SEN,levels=design)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
MMCL_LenNONRES_SEN <- topTable(fit2, adjust="BH", confint=TRUE, number=nrow(MVal.cL.ncorrMM.Min))
MMCL_LenNONRES_SEN <- merge (MMCL_LenNONRES_SEN, ChromHMM, by="row.names", all.x=TRUE, sort=FALSE)
MMCL_LenNONRES_SEN <- data.frame(MMCL_LenNONRES_SEN[,-1], row.names=MMCL_LenNONRES_SEN[,1])
write.table(MMCL_LenNONRES_SEN,file="C:/Users/johnson/Desktop/Meth.analysis/Results - Cell Lines/MMCL_LenNONRES_SEN.Nov17.txt",sep="\t")

