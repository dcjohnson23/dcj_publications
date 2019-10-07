hh #Barplot of R2 values
r2 <- table(perl.data2$R.2)
barplot(r2, main="R2 values", xlab="Number of cases")

# Grouped Bar Plot
counts2 <- table(perl.data$Mutations, perl.data$R.2.98bin)
barplot(counts2, main="MM pts counts by Neutral and Mutations",
ylab="Number of cases", col=c("darkblue","red"),
legend = rownames(counts2), cex.names=0.65, ylim=c(0,100), args.legend = list(x = "topright", cex=.5))

port2 <- prop.table(counts2,2)
barplot(port2, main="MM pts proportion by Neutral and Chr17",
        ylab="Percentage of cases", col=c("darkblue","red"),
        legend = rownames(port2),cex.names=0.65,las=2, args.legend = list(x = "topright", cex=.5))

counts1 <- table(perl.data.436$R.2.98, perl.data.436$age.group)
port1 <- prop.table(counts1,2)
barplot(counts1, main="Number of MM pts by Neutral and age group",
        ylab="Number of cases", col=c("black","red"))

barplot(port1, main="MM pts proportion by Neutral and age group",
        ylab="Percentage of cases", col=c("black","red"))


counts4 <- table(perl.data$R.2.98, perl.data$chr.abs.bin2)
barplot(counts3, main="MM pts proportion by Neutral and Response",
        ylab="Number of cases", col=c("ivory4","ivory3"),
        legend = c("Neutral","Non-neutral"), ylim=c(0,380),cex.names=0.65,las=2,args.legend = list(x = "topright", cex=.5))
port4 <- prop.table(counts4,2)

counts3 <- table(perl.data$R2Call, perl.data$Response.induction)
port3 <- prop.table(counts3,2)
barplot(port3, main="MM pts proportion by Neutral and Response - non intensive",
        ylab="Percentage of cases", col=c("black","red"),
        legend = c("Neutral","Non-neutral"),cex.names=0.65,las=2,args.legend = list(x = "topright")

barplot(port4, main="MM pts distribution by R.2 versus any chr abs (not 1q21 gain",
        ylab="Percentage of cases", col=c("darkblue","red"),
        legend = rownames(counts4),cex.names=0.65,las=2,args.legend = list(x = "topright", cex=.5))

response.all=factor(perl.data$Response.induction, levels=levels(perl.data$Response.induction)[c(1,6,2,5,3,4)])
response.in=factor(perl.data.int$Response.induction.1, levels=levels(perl.data.int$Response.induction.1)[c(1,5,4,2,3)])
response.no.int=factor(perl.data.no.int$Response.induction.1, levels=levels(perl.data.no.int$Response.induction.1)[c(1,5,4,2,3)])
counts3 <- table(response.no.int, perl.data.no.int$R298)
port3 <- prop.table(counts3,2)
barplot(counts4, main="MM pts counts by Neutral and Response non-intensive",
        ylab="Counts of cases", col=c("darkblue","red","green","yellow","grey","orange"),
        legend = rownames(counts4),ylim=c(0,200), args.legend = list(x = "topright", cex=.5))

barplot(port4, main="MM pts proportion by Neutral and Response non-intensive",
        ylab="Percentage of cases", col=c("darkblue","red","green","yellow","grey","orange"),
        legend = rownames(port4), args.legend = list(x = "topright", cex=.5))

par(xpd=TRUE)
barplot(sig.table, col=c("darkblue","grey","mistyrose2","orange","blue","lightblue","red","green","lightpink1","cyan","gray71"))
par(mar=c(10.2, 8.2, 8.2, 16.2), xpd=TRUE)
legend("topright", inset=c(-0.4,0), legend=row.names(sig.table), pch=c(1,3), title="Gene signatures", par(xpd=TRUE))


bb <- barplot(tt, col=c("grey60","grey80"))
text(bb,tt[1,]-4,labels=tt[1,],cex=.8)
text(bb,colSums(tt)-4,labels=tt[2,],cex=0.8)


#Multivariant analysis - Quick
#myvars <- c("R.2","R.2.group", "chr17", "myc.tnx", "Intensive", "age", "Male", "os.months", "os.status", "pfs.months", "pfs.status","tnx", "iss")
#perl.data.mv <- perl.data[myvars]
#perl.data.noNAs <- na.omit(perl.data.mv)
#my11.os.R2.fit <- coxph(Surv(os.months,os.status)~R.2.group+Intensive+age+Male+iss,data=perl.data.noNAs, method="exact")
#stepAIC(my11.os.R2.fit, direction=c("backward"))

# 2x2 Factorial MANOVA with dependent variables
#Y <- cbind(y1,y2,y3)
#fit <- manova(Y~A*B)
#summary(fit, test="Pillai")


#cor <- cor(perl.data.noNAs, use="complete.obs", method="kendall")