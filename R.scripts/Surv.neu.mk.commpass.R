# Load survival library and set temporary directory
library(survival)
library(KMsurv)
library(coin)

# Create custom survival function.  This thing was originally written by 
# "Xiaochun Li" at the Dana Farber in Harvard - thanks, Google!
smed <- function(x) { 
  ox <- capture.output(print(x)) 
  n <- length(ox) 
  tmp <- t(sapply(ox[4:n],           
                  function(l) strsplit(l, split=' +')[[1]])) 
  nres <- strsplit(ox[3],split=' +')[[1]][2:6] 
  res <- matrix(as.numeric(tmp[,2:6]), ncol=5, 
                dimnames=list(tmp[,1], nres)) 
  res 
}

# Read in the correctly formatted data file (created from the database and written to a file named
# using the current date and time in the temporary directory.  I renamed it to make it clearer).
# Numbers are crunched and the pdf output file is then instantiated.  To work without writing to a file, 
# you could just comment out the pdf line.

perl.data.new <- read.table("C:/Users/johnson/Desktop/CoMMpass/R2.coMMpass.2.csv", sep=",",header = TRUE, na.strings="NA")
perl.data <- subset(perl.data.new, !is.na(R2Call))
perl.data.no.int <- subset(perl.data, pathway=="NonIntensive")
perl.data.int <- subset(perl.data, pathway=="Intensive")
perl.data.no.int <- subset(perl.data.no.int, Intensity=="standard")

perl.surv<- Surv(os.months, os.status) ~ R2Call
perl.diff <- survdiff(perl.surv, data=perl.data.no.int, na.action=na.omit)
perl.fit <- survfit(perl.surv, data=perl.data.no.int, na.action=na.omit)
this.sf1 <- smed(perl.fit)

pdf("Neu.vs.No_Neu.R2.98.no.int.os.pdf" , paper='a4', width=6.4 , height=6.4 , pointsize=10)
# The plot is produced and logrank p calculated, then pasted on to the plot, along with a title and
# other data, like how many samples there are in each group.  Also, a legend is produced.
# Finally, "dev.off()" disconnects the output device (the pdf file).

plot(perl.fit,col=c("red", "black", "dark blue"," purple", "orange", "yellow"), lwd=c("3","3"),
     pch= c("0","1"),mark="|",
     lty =c(1,1), xmax=100, xlab='Survival (months)',ylab='Proportion of Patients')
mydf<-length(perl.diff$obs)-1
this_x<-1-pchisq(perl.diff$chisq,mydf)
opvals = 0
if(this_x<0.001 && opvals!=1){mysub<-paste("Logrank P","<0.001",sep=" ")
}else{
mysub<-paste("Logrank P= ",sprintf("%5.4f",1-pchisq(perl.diff$chisq,mydf)),sep=" ")}
title(main=" Overall Survival: R2>0.98, no-int", sub=mysub)
this.medians<-paste("Medians: ","Neutral = ",this.sf1[2,3]," ",": Non-Neutral = ",this.sf1[1,3]," "," months",sep="")
mtext(paste("No of Patients: Neutral = ",perl.diff$n[2]," ",": Non-Neutral = ",perl.diff$n[1]," ",sep=""),side=4)
mtext(this.medians, line=1, side=4)
legend(0.2, 0.15, border = "white",legend=c("Neutral","Non-Neutral"), bty="n", col=c( "black", "red", "dark blue"," purple", "orange", "yellow"), lwd=c("3","3"), lty=c(1,1), cex=1.2)
dev.off()