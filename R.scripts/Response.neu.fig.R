perl.data.new <- read.table("C:/Users/johnson/Documents/Neutral cancers/Jan.analysis/my11_neu_subtype.plus4.jan.csv", sep=",",header = TRUE, na.strings="NA")
perl.data <- subset(perl.data.new, !is.na(R2Call))
perl.data.no.int <- subset(perl.data, Intensive==0)
perl.data.int <- subset(perl.data, Intensive==1)


res.perl.data.int <- subset(perl.data.int, !is.na(res))

response.int <- factor(res.perl.data.int$res, levels=levels(res.perl.data.int$res)[c(1,5,4,2,3)])
counts5 <- table(res.perl.data.int$R2Call, response.int)
port5 <- prop.table(counts5,2)
barplot(counts5, main="Number of MM pts by Neutral and Response - int",
        ylab="Number of cases", col=c("black","red"))
res.int <- barplot(port5, main="MM pts proportion by Neutral and Response - int",
        ylab="Percentage of cases", col=c("black","red"))


res.perl.data.no.int <- subset(perl.data.no.int, !is.na(res))

response.no.int <- factor(res.perl.data.no.int$res, levels=levels(res.perl.data.no.int$res)[c(1,5,4,2,3)])
counts6 <- table(res.perl.data.no.int$R2Call, response.no.int)
port6 <- prop.table(counts6,2)
barplot(counts6, main="Number of MM pts by Neutral and Response - non int",
        ylab="Number of cases", col=c("black","red"))
barplot(port6, main="MM pts proportion by Neutral and Response - non int",
        ylab="Percentage of cases", col=c("black","red"))