library(metafor)

# Load input.tyx
IMIDFP <- read.table("FP.input.IMID.txt",as.is=T, header=T, sep="\t")
IMIDnoFP <- read.table("FP.input.nonIMID.OS.txt",as.is=T, header=T, sep="\t")

# # If OR's and 95% CI's, run these commented out lines
##data3$beta <- log(data3$OR)
##data3$se   <- (log(data3$U95)-log(data3$L95))/(2*1.96)
IMIDFP$beta <- log(IMIDFP$OR)
IMIDnoFP$beta <- log(IMIDnoFP$OR)


# Assign values for plotting
labs <- IMIDFP$Study
labs2 <- IMIDnoFP$Study

# Combine data into summary estimate
res  <- rma(yi=IMIDFP$beta, sei=IMIDFP$SE, method="FE")
res2 <- rma(yi=IMIDnoFP$beta, sei=IMIDnoFP$SE, method="FE")

summary(res)
summary(res2)

##Combined plot
par(mfrow=c(2,1))
forest(res, transf=exp, refline=2.00, xlab="Odds Ratio (95%CI)",alim=c(0,4), xlim=c(-3,6),  slab=labs, mlab="Meta")
segments(1, -1.8, 1, res$k+1)
mtext(expression(plain("IMID-based trials")),adj=0)

text(-2, -1.5, pos=4, cex=0.75, bquote(paste("RE Model for Subgroup (Q = ",.(formatC(res$QE, digits=2, format="f")), ", df = ", .(res$k - res$p), ", p = ", .(formatC(res$QEp, digits=2, format="f")), "; ", I^2, " = ",.(formatC(res$I2, digits=1, format="f")), "%)")))

text(4, -1.5, pos=4, cex=0.75, bquote(paste("Meta FE (p = ", .(formatC(res$pval, digits=6, format="f")), ")")))

forest(res2, transf=exp, refline=1.04, xlab="Odds Ratio (95%CI)",alim=c(0,4), xlim=c(-3,6),  slab=labs2, mlab="Meta")
segments(1, -1.8, 1, res2$k+1)
mtext(expression(plain("Combination or non IMID-based trials")),adj=0)

text(-2, -1.5, pos=4, cex=0.75, bquote(paste("RE Model for Subgroup (Q = ",.(formatC(res2$QE, digits=2, format="f")), ", df = ", .(res2$k - res2$p), ", p = ", .(formatC(res2$QEp, digits=2, format="f")), "; ", I^2, " = ",.(formatC(res2$I2, digits=1, format="f")), "%)")))

text(4, -1.5, pos=4, cex=0.75, bquote(paste("Meta FE (p = ", .(formatC(res2$pval, digits=3, format="f")),")")))


# Load input.Examples
mtext(expression(paste("Cu"^"2+","at EC50",sep="")))
mtext(expression(paste( plain("Cu") ^ plain("2+"), plain(" at EC50") )), side=2, line = 4, padj=1, at=30, cex=1.2)
mtext(expression(italic("Cu"^"2+","at EC50",sep="")))
mtext(expression( italic("h") ^ italic("2")), side=2, line = 2)

