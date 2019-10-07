library(minfi)
library(sva)
library(conumee)
library(CopyNumber450kData)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(FlowSorted.Blood.450k)
library(RColorBrewer)
library(limma)
library(Gviz)
setwd("/500_Exomes_Project/conumee.epic/")

## -------read in epic files - place sheet in data.dir-----------------------------------------------------------------##
data.dir <- "/500_Exomes_Project/epic.IDATs"
epic.CL <- read.metharray.sheet(data.dir, pattern="MM_pts.epic.csv")
RGsetTest <- read.metharray.exp(targets=epic.CL)
RGsetTest@annotation <- c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b2.hg19")

data.dir2 <- "/500_Exomes_Project/epic.IDATs"
epic.Con <- read.metharray.sheet(data.dir2, pattern="controls.epic.csv")
RGsetCon <- read.metharray.exp(targets=epic.Con)
RGsetCon@annotation <- c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b2.hg19")

MsetTest <- preprocessIllumina(RGsetTest)
data(RGcontrolSetEx)
MsetControls <- preprocessIllumina(RGcontrolSetEx)

MsetControls2 <- preprocessIllumina(RGsetCon)
data(exclude_regions)

MLPA <- read.table("/500_Exomes_Project/conumee.450.nov/MLPA5.100000.Igs.csv", sep=",",header = TRUE)
MLPA.gr <- as.data.frame.matrix(MLPA)
seqnames <- MLPA.gr$seqnames
start <- MLPA.gr$start
end <- MLPA.gr$end
names <- MLPA.gr$names
score <- MLPA.gr$score
strand <- MLPA.gr$strand
gr <- GRanges(seqnames=seqnames, ranges=IRanges(start=start,end=end), names=names, score=score, strand=strand)
anno <- CNV.create_anno(exclude_regions = exclude_regions, array_type = "overlap", detail_regions = gr)
minfi.data <- CNV.load(MsetTest)
controls.data <- CNV.load(MsetControls2)
minfi.controls <- pData(MsetTest)$status == "normal"

MMg <- read.table("/500_Exomes_Project/conumee.450.nov/MM.genes.01.csv", sep=",",header = TRUE)
MM.gr <- as.data.frame.matrix(MMg)
seqnames <- MM.gr$seqnames
start <- MM.gr$start
end <- MM.gr$end
names <- MM.gr$names
score <- MM.gr$score
strand <- MM.gr$strand
gr2 <- GRanges(seqnames=seqnames, ranges=IRanges(start=start,end=end), names=names, score=score, strand=strand)
anno2 <- CNV.create_anno(exclude_regions = exclude_regions, array_type = "EPIC", detail_regions = gr2)

gr.con <- GRanges(seqnames=c("chr2","chr2","chr4","chr4","chr10","chr10","chr12","chr12","chr14","chr16", "chr18","chr18","chr20","chr21"),
                  ranges=IRanges(start=c(61244812,189839099,42410392,110661848,11962021,111767711,3186521,117645947,74946643,87684443,3451591,21111463,35520227,40547372),end=c(61279125,189877472,42659122,110723335,12084840,111895323,3395730,117742857,74960084,88878432,3458406,21166581,35580246,40555440)),
                  names=c("PEX13","COL3A1","ATP8A1","CFI","UPF2","ADD3","TSPAN9","NOS1","NPC2","ABAT","TGIF1","NPC1","SAMHD1","PSMG1"), 
                  scores=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),
                  strand=c("+","+","-","-","-","+","+","-","-","+","+","-","-","-")
)
anno3 <- CNV.create_anno(exclude_regions = exclude_regions, array_type = "EPIC", detail_regions = gr.con)

clist <- epic.CL$Sample_Name

lapply(clist, function(z) {
x <- CNV.fit(minfi.data[z], controls.data, anno)
x <- CNV.bin(x)
x <- CNV.detail(x)
x <- CNV.segment(x)
y <- CNV.fit(minfi.data[z], controls.data, anno2)
y <- CNV.bin(y)
y <- CNV.detail(y)
y <- CNV.segment(y)
q <- CNV.fit(minfi.data[z], controls.data, anno3)
q <- CNV.bin(q)
q <- CNV.detail(q)
q <- CNV.segment(q)

pdf(file=(paste(z,".pdf")), height = 9, width = 18)
CNV.genomeplot(y)
CNV.genomeplot(y, chr = "chr1")
CNV.genomeplot(y, chr = "chr2")
CNV.genomeplot(y, chr = "chr3")
CNV.genomeplot(y, chr = "chr4")
CNV.genomeplot(y, chr = "chr5")
CNV.genomeplot(y, chr = "chr6")
CNV.genomeplot(y, chr = "chr7")
CNV.genomeplot(y, chr = "chr8")
CNV.genomeplot(y, chr = "chr9")
CNV.genomeplot(y, chr = "chr10")
CNV.genomeplot(y, chr = "chr11")
CNV.genomeplot(y, chr = "chr12")
CNV.genomeplot(y, chr = "chr13")
CNV.genomeplot(y, chr = "chr14")
CNV.genomeplot(y, chr = "chr15")
CNV.genomeplot(y, chr = "chr16")
CNV.genomeplot(y, chr = "chr17")
CNV.genomeplot(y, chr = "chr18")
CNV.genomeplot(y, chr = "chr19")
CNV.genomeplot(y, chr = "chr20")
CNV.genomeplot(y, chr = "chr21")
CNV.genomeplot(y, chr = "chr22")
CNV.detailplot(x, name = "IGL-IMGT-HG19")
CNV.detailplot(x, name = "IGK-IMGT-HG19")
CNV.detailplot(x, name = "IGH-IMGT-HG19")
CNV.detailplot_wrap(y)
dev.off()
CNV.write(y, what = "detail",file=(paste(z, "CNVdetail.450.MMgenes.adj.txt")))
CNV.write(x, what = "detail",file=(paste(z, "CNVdetail.450.mlpa.adj.txt")))})



