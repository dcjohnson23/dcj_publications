setwd("C:/")

library(made4)
library(affyPLM)
library(MatrixEQTL)
library(SNPRelate)
library(gdata)
library(affy)
library(hgu133plus2hsentrezg.db)


load('gcrma_mas5_panp_gwas_validataion_UK_022013.Rdata')

data.file <- read.csv2("UK_GEP_David_186_co.txt", as.is=T, sep="\t") 
data2.file <- data.file[-grep("MM094.cel", data.file$mm_id), ] 
data2.file <- data2.file[-grep("MM003.cel", data2.file$mm_id), ]
data2.file <- data2.file[-grep("MM224.cel", data2.file$mm_id), ]


expr <- exprs(gcrma.gwas.uk)


colnames(expr)<- data.file$Cel.name

expr <- expr[,match(data2.file$Cel.name, colnames(expr))]


zero.var <-  which(apply(expr,1,var)==0)
expr.red <- expr[-match(names(zero.var), rownames(expr)), ]

cut.var <-  which(apply(expr.red,1,var)<0.1)
cut.3.5 <- which((apply(expr.red,1,function(x) length(which(x>=3.5)))/(dim(expr)[2])*100) <= 10)

cut.affx <- grep("AFFX", rownames(expr.red))

cut <- unique(c(cut.var, cut.3.5, cut.affx))

expr.cut <- expr.red[-cut,]

annot <- read.csv("HGU133Plus2_Hs_ENTREZG_desc.annot", header=F, as.is=T, sep="\t")
desc <- read.csv("HGU133Plus2_Hs_ENTREZG_desc.txt", header=T, as.is=T, sep="\t")
mapping <- read.csv("HGU133Plus2_Hs_ENTREZG_mapping.txt", header=T, as.is=T, sep="\t")

gstart <- mget(annot$V1, env=hgu133plus2hsentrezgCHRLOC)
gstart <- unlist(sapply(gstart, "[[",1))

gend <- mget(annot$V1, env=hgu133plus2hsentrezgCHRLOCEND)
gend <- unlist(sapply(gend, "[[",1))

gchr <- mapping$Chr[match(annot$V1, mapping$Probe.Set.Name)]
geneloc <- data.frame(geneid = annot$V1,
                      chr = gchr,
                      s1 = abs(gstart),
                      s2 = abs(gend)
)
write.table(geneloc, "geneloc.txt", quote=F, row.names=FALSE)

yy <- geneloc$geneid[which(geneloc$chr=="X" | geneloc$chr=="Y" )]

expr.fin <- expr.cut[-na.omit(match(yy, rownames(expr.cut))),]

x <- names(which(is.na(gstart)))

expr.fin <- expr.fin[-na.omit(match(x,rownames(expr.fin))),]

gep = expr.fin

bed.fn <- "my9_illumina_5_13_b37_183.bed"
bim.fn <- "my9_illumina_5_13_b37_183.bim"
fam.fn <- "my9_illumina_5_13_b37_183.fam"

snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "UK.gds")

snpgdsSummary("UK.gds")

genofile <- openfn.gds("UK.gds")

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
snp.id <- read.gdsn(index.gdsn(genofile, "snp.id"))
snp.allele <- read.gdsn(index.gdsn(genofile, "snp.allele"))
snp.sex <- read.gdsn(index.gdsn(genofile, "sample.annot"))$sex
snp.position <- read.gdsn(index.gdsn(genofile, "snp.position"))
snp.chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
snp.allele.a1 <- unlist(strsplit(snp.allele, split="/"))[seq(1,length(snp.allele)*2,by=2)]
snp.allele.a2 <- unlist(strsplit(snp.allele, split="/"))[seq(2,length(snp.allele)*2,by=2)]
g <- read.gdsn(index.gdsn(genofile, "genotype"))
g <- t(g)
colnames(g) <- sample.id
rownames(g) <- snp.id

g <- g[,match(data2.file$IID, sample.id)]
snp.sex <- gsub("F", 1, snp.sex)
snp.sex <- gsub("M", 2, snp.sex)

gep <- gep[,match(data2.file$Cel.name, colnames(gep))]

colnames(gep) <- colnames(g)

snploc <- data.frame(SNP = snp.id,
                     chr = snp.chr,
                     pos = snp.position
)

write.table(snploc, "snpsloc.txt", quote=F, row.names=FALSE)

# remove snps & genes on x,y chr
xx <- snploc$SNP[which(snploc$chr==23 | snploc$chr==24 | snploc$chr==25)]


g <- g[-na.omit(match(xx, rownames(g))),]
snp.allele.a1 <- snp.allele.a1[-na.omit(match(xx, rownames(g)))]
snp.allele.a2 <- snp.allele.a2[-na.omit(match(xx, rownames(g)))]


write.table(rownames(gep), "probesets_UK.txt", quote=F, col.names=NA)

write.table(g, "SNP.txt", quote=F, col.names=NA)
write.table(gep, "GE.txt", quote=F, col.names=NA)

co <- data.frame(id = data2.file$IID[match(colnames(gep), data2.file$IID)],
                 sex = snp.sex)

write.table(t(co), "covariates_UK.txt", quote=F, col.names=FALSE, sep="\t")

useModel = modelLINEAR # modelANOVA or modelLINEAR or modelLINEAR_CROSS 
SNP_file_name = "SNP.txt"
expression_file_name = "GE.txt"

# in case of no covariates set the variable covariates_file_name to character()
covariates_file_name = "UKpeercot2a.txt" # character()
output_file_name = "eQTL_UK.txt"
output_file_name_cis = "eQTL_UK_cis.txt"
output_file_name_tra = "eQTL_UK_tra.txt"

# Setting the threshold to a high value for a large dataset may cause excessively large output (many gigabytes)
pvOutputThreshold = 5e-4
pvOutputThreshold_cis = 1e-3
pvOutputThreshold_tra = 5e-7

# Finally, define the covariance matrix for the error term. This a seldom used parameter. 
# If the matrix is a multiple of identity, set the parameter to numeric()
errorCovariance = numeric()

snpspos = read.table("snpsloc.txt", header = TRUE, stringsAsFactors = FALSE)
genepos = read.table("geneloc3.txt", header = TRUE, stringsAsFactors = FALSE)

# cis dist = 1Mb
cisDist = 1e4

snps = SlicedData$new()
snps$fileDelimiter = " " # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values;
snps$fileSkipRows = 1 # one row of column labels
snps$fileSkipColumns = 1 # one column of row labels
snps$fileSliceSize = 2000 # read file in pieces of 2,000 rows
snps$LoadFile(SNP_file_name)

# load gep data
gene = SlicedData$new()
gene$fileDelimiter = " " # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values;
gene$fileSkipRows = 1 # one row of column labels
gene$fileSkipColumns = 1 # one column of row labels
gene$fileSliceSize = 2000 # read file in pieces of 2,000 rows
gene$LoadFile(expression_file_name);

# Load covariates
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t" # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values;
cvrt$fileSkipRows = 1 # one row of column labels
cvrt$fileSkipColumns = 1 # one column of row labels
cvrt$fileSliceSize = 2000 # read file in one piece
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name)
}

training.co = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = 10,
  min.pv.by.genesnp = TRUE
)  

training <- training.co
# ----------------------------------------------------------------------
# CIS/TRANS Tables
# ----------------------------------------------------------------------

# compute test statistiks for sig. cis/trans effects
teststat <- function(snp,gene) {
  f <- lm(as.numeric(gep[gene, ]) ~ as.numeric(g[snp,]))
  fbeta <- summary(f)$coefficients[2,1] # BETA
  fse <- summary(f)$coefficients[2,2] # SE
  fr2 <- summary(f)$r.squared # R2
  nmiss <- length(which(g[snp,]==3))
  return(list(fbeta,fse,fr2,nmiss))
}

# cis table		
cis <- training$cis$eqtls
cis$chr <- genepos$chr[match(cis$gene, genepos$geneid)]
cis$chr_snp <- snpspos$snp.chr[match(cis$snps, snpspos$SNP)]
cis$bp <- snpspos$snp.postion[match(cis$snps, snpspos$SNP)]
cis$desc <- desc$Description[match(cis$gene, desc$Probe.Set.Name)]
cis$symbol <- as.vector(unlist(mget(as.vector(cis$gene), env=hgu133plus2hsentrezgSYMBOL)))
cis$TSS <- genepos$s1[match(cis$gene, genepos$geneid)]
cis$SNPpos <- snpspos$pos[match(cis$snps, snpspos$SNP)]
# cis$SNPlocation <- gwas.anno$location[match(cis$snps, gwas.anno$rs_id)]

res <- numeric(length=4*dim(cis)[1])
resseq <- seq(1,4*dim(cis)[1],by=4)
for(i in 1:dim(cis)[1]) {
  x <- unlist(teststat(as.vector(cis$snps[i]), as.vector(cis$gene[i])))
  res[resseq[i]] <- x[1]
  res[resseq[i]+1] <- x[2]
  res[resseq[i]+2] <- x[3]
  res[resseq[i]+3] <- x[4]
}

cis$BETA <- res[seq(1,4*dim(cis)[1],by=4)]
cis$SE <- res[seq(2,4*dim(cis)[1],by=4)]
cis$R2 <- res[seq(3,4*dim(cis)[1],by=4)]
cis$NMISS <- dim(gep)[2]-res[seq(4,4*dim(cis)[1],by=4)]

cis.uq <- cis[match(unique(cis$gene), cis$gene),]

cis.meta <- data.frame(SNP = paste(cis$snps,cis$gene,sep="."),
                       OR = cis$BETA,
                       SE = cis$SE,
                       P = cis$pvalue,
                       CHR = cis$chr,
                       BP = cis$SNPpos,
                       A1 = snp.allele.a1[match(cis$snps, rownames(g))],
                       A2 = snp.allele.a2[match(cis$snps, rownames(g))],
                       N = cis$NMISS
)

write.table(cis.meta, "2014.02.13.cis.meta.UK_test.assoc", col.names=TRUE, sep="\t", quote=F, row.names=FALSE)


# trans table		
trans <- training$trans$eqtls
trans$chr <- genepos$chr[match(trans$gene, genepos$geneid)]
trans$chr_snp <- snpspos$snp.chr[match(trans$snps, snpspos$SNP)]
trans$bp <- snpspos$snp.position[match(trans$snps, snpspos$SNP)]
trans$desc <- desc$Description[match(trans$gene, desc$Probe.Set.Name)]
trans$symbol <- as.vector(mget(as.vector(trans$gene), env=hgu133plus2hsentrezgSYMBOL))
trans$TSS <- genepos$s1[match(trans$gene, genepos$geneid)]
trans$SNPpos <- snpspos$pos[match(trans$snps, snpspos$SNP)]
# trans$SNPlocation <- gwas.anno$location[match(trans$snps, gwas.anno$rs_id)]

res <- numeric(length=4*dim(trans)[1])
resseq <- seq(1,4*dim(trans)[1],by=4)
for(i in 1:dim(trans)[1]) {
  x <- unlist(teststat(as.vector(trans$snps[i]), as.vector(trans$gene[i])))
  res[resseq[i]] <- x[1]
  res[resseq[i]+1] <- x[2]
  res[resseq[i]+2] <- x[3]
  res[resseq[i]+3] <- x[4]
}

trans$BETA <- res[seq(1,4*dim(trans)[1],by=4)]
trans$SE <- res[seq(2,4*dim(trans)[1],by=4)]
trans$R2 <- res[seq(3,4*dim(trans)[1],by=4)]
trans$NMISS <- dim(gep)[2]-res[seq(4,4*dim(trans)[1],by=4)]

trans.uq <- trans[match(unique(trans$gene), trans$gene),]

trans.meta <- data.frame(SNP = paste(trans$snps,trans$gene,sep="."),
                         OR = trans$BETA,
                         SE = trans$SE,
                         P = trans$pvalue,
                         CHR = trans$chr,
                         BP = trans$SNPpos,
                         A1 = snp.allele.a1[match(trans$snps, rownames(g))],
                         A2 = snp.allele.a2[match(trans$snps, rownames(g))],
                         N = trans$NMISS
)

write.table(trans.meta, "2014.02.13.trans.meta.UK_test.assoc", col.names=TRUE, sep="\t", quote=F, row.names=FALSE)



# tables
write.table(as.matrix(cis), "2014.02.13.matrix_eQTL_UK_cis.txt", quote=F, col.names=NA, sep="\t")
write.table(as.matrix(cis.uq), "2014.02.13.matrix_eQTL_UK_cis_uq.txt", quote=F, col.names=NA, sep="\t")

write.table(as.matrix(trans), "2014.02.13.matrix_eQTL_UK_trans.txt", quote=F, col.names=NA, sep="\t")
write.table(as.matrix(trans.uq), "2014.02.13.matrix_eQTL_UK_trans_uq.txt", quote=F, col.names=NA, sep="\t")





save.image("matrix_eqtl_UK_customCDF_12_06_2013.Rdata")

