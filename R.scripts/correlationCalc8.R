data1 <- try(read.table("VAF.txt"))
data3 <- subset(data1, data1$V1 <= 0.6)
data <- subset(data1, data1$V1 <= 0.30)
data <- subset(data, data$V1 >= 0.12)
data2 <- try(read.table("NAME.txt"))
data2 <- unique(data2)
if(inherits(data, "try-error")) {
   write("NA", file = "correlations8.txt", append = TRUE)
} else {

data$V2<- 1/data$V1
data <- data[order(data$V2),] 
data$V3 <- 1
data$Mf <- cumsum(data$V3)
r2 <- cor(data$V2,data$Mf)
noc <- length(data$V1)
be <- noc/5.3
ae <- noc-(be*8.33)
bob <- paste(data2$V1)
xAxis2 <- c("0.33","0.25", "0.2", "0.17", "0.14", "0.125")
pdf(file=(paste(data2$V1,".pdf")), height = 9, width = 18)
par(mfrow=c(2,1))
hist(data3$V1, breaks=seq(0,0.6,0.01), main=data2$V1, xlim=c(0, 0.6), xlab="Allelic frequency", ylab="Number of mutations")
plot(data$V2, data$Mf, xlab="Allelic frequency (1/f)", ylab="cumulative subclonal mutations", xaxt="n", ylim=c(0,noc), xlim=c(3, 8.5) )
axis(2)
axis(1, at=3:8, labels = xAxis2)
abline(a=ae, b=noc/5.3, col = "red")
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
text(x = 8.3, y = 1, labels = mylabel)
dev.off()
write(noc, file = "mutationsNumbers8.txt", append = TRUE)
write(bob, file = "names8.txt", append = TRUE)
write(r2, file = "correlations8.txt", append = TRUE)
}
