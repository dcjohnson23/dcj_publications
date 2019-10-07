library(snpStats)
mydata <- read.table ("/Clin.demo/UK.bone.Nov.14/Meta.all/bone.meta.all.log.4.tab.head.txt", header = T, sep="")
png("Bonedis.meta.90.2.png")
qq.chisq(-2 * log(mydata$P), df = 2, pvals = TRUE, overdisp = TRUE, slope.one=TRUE)
dev.off()

library(qqman)
png("Bonedis.meta.mann.plot.90.2.png")
manhattan(mydata80, main = "Manhattan Plot", cex = 0.5, cex.axis = 0.8)
dev.off()