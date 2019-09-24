library(survivalROC)


#Parameters
#Stime Event time or censoring time for subjects
#status Indicator of status, 1 if death or event, 0 otherwise
#marker Predictor or marker value
#entry Entry time for the subjects
#predict.time Time point of the ROC curve
#cut.values marker values to use as a cut-off for calculation of sensitivity and specificity
#method Method for fitting joint distribution of (marker,t), either of KM or NNE, the
#default method is NNE
#lambda smoothing parameter for NNE
#span Span for the NNE, need either lambda or span for NNE
#window window for NNE, either of symmetric or asymmetric

#create ROC calculation
SurvROC= survivalROC(My11.demo.int.cpg$OS.months, status=My11.demo.int.cpg$OS.status, 
            marker=My11.demo.int.cpg$cg00216138, 
            predict.time=24.0,
            method = "KM")

SurvROC= survivalROC(My11.demo.int.cpg$OS.months, status=My11.demo.int.cpg$OS.status, 
                     marker=My11.demo.int.cpg$cg00216138, 
                     predict.time=24.0,
                     method = "NNE")




str(SurvROC)

cutoff<-SurvROC$cut.values[SurvROC$TP-SurvROC$FP==min((SurvROC$TP-SurvROC$FP))]

str(cutoff)

SurvROC.1 = survivalROC(Stime=P425_01_pfs_os_intensive$os_months, status=P425_01_pfs_os_intensive$os_status, 
                     marker=P425_01_pfs_os_intensive$TP53.probe.08304.L21074_01, 
                     predict.time=24.0,
                     method = "KM")
str(SurvROC.1)

#plot
plot(SurvROC$FP, SurvROC$TP, type = "l", xlim=c(0,1), ylim=c(0,1),
     xlab=paste ("FP", "\n", "AUC =",round(SurvROC$AUC,3)),
     ylab="TP", main="CoSeg, Method = Kaplan-Meier\n 18 Months")
abline(0,1)


#plot 2 lines
plot(SurvROC$FP, SurvROC$TP, col="red", type = "l", xlim=c(0,1), ylim=c(0,1),
     xlab=paste ("FP", "\n", "AUC =",round(SurvROC$AUC,3)),
     ylab="TP", main="CoSeg, Method = Kaplan-Meier\n 18 Months")
abline(0,1)
lines (SurvROC.1$FP, SurvROC.1$TP, col ="blue",lty =2)
legend(0.45, 0.4, legend = c("CDKN2C deletion; AUC =", round(SurvROC$AUC,3), 
                            "CDKN2C & COLL11A1; AUC =", round(SurvROC.1$AUC,3)), col = c("red","blue"),
       lty = c(1,2), bty ="n")


#ROC over time #OS
AUC_CoSeg = NULL; AUC_CoSegInt = NULL
for(t in 1:80) {
  AUC_CoSeg = c(AUC_CoSeg, survivalROC(Stime=P425_01_pfs_os_intensive$os_months, 
                                       status=P425_01_pfs_os_intensive$os_status, 
                                       marker=P425_01_pfs_os_intensive$X..CDKN2C.probe.14652.L16304..01.051.212271_01,  
                                       predict.time=t, method = "KM")$AUC )
  AUC_CoSegInt = c(AUC_CoSegInt, survivalROC(Stime=P425_01_pfs_os_intensive$os_months, 
                                       status=P425_01_pfs_os_intensive$os_status, 
                                       marker=P425_01_pfs_os_intensive$FAM46C.probe.18949.L24912..01.117.966975_01,  
                                       predict.time=t, method = "KM")$AUC )
}

#plot
plot(1:80, AUC_CoSeg, type= "l", xlab = "Time", ylab ="AUC",
     ylim = c(0.50,0.60), col = "red", main="AUC (OS) over time \n Method: KM")
lines (1:80, AUC_CoSegInt, col ="blue",lty =2)
legend(35, 0.60, legend = c("CDKN2C", "FAM46C"), col = c("red","blue"),
       lty = c(1,2), bty ="n")

#plot
plot(min(poisson$OS_months):36, AUC_CoSeg, type= "l", xlab = "Time", ylab ="AUC",
     ylim = c(0.55,0.75), col = "red")
lines (min(poisson$OS_months):36, AUC_CoSegInt, col ="blue",lty =2)
legend(28, 0.75, legend = c("CoSeg", "CoSegInt"), col = c("red","blue"),
       lty = c(1,2), bty ="n")

