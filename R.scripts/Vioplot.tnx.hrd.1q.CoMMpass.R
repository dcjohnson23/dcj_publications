library(vioplot)

IA9_fig_vio <- IA9_long2_tnx[c(1,6,7,36,37,38,39,40,41,42)]
IA9_fig_vio <- subset(IA9_fig_vio, !is.na(Any.Tnx))
IA9_fig_vio2 <- subset(IA9_fig_vio, !is.na(Any.Tnx.FNAH)) 

x1 <- IA9_fig_vio$New.R2[IA9_fig_vio$tnx.1q==11]
x2 <- IA9_fig_vio$New.R2[IA9_fig_vio$tnx.1q=="11.1qgain"]
x3 <- IA9_fig_vio$New.R2[IA9_fig_vio$tnx.1q==4]
x4 <- IA9_fig_vio$New.R2[IA9_fig_vio$tnx.1q=="4.1qgain"]
x5 <- IA9_fig_vio$New.R2[IA9_fig_vio$tnx.1q==20]
x6 <- IA9_fig_vio$New.R2[IA9_fig_vio$tnx.1q==16]
x7 <- IA9_fig_vio$New.R2[IA9_fig_vio$tnx.1q==6]
x8 <- IA9_fig_vio$New.R2[IA9_fig_vio$tnx.1q=="HRD"]
x9 <- IA9_fig_vio$New.R2[IA9_fig_vio$tnx.1q=="HRD.11q25.gain"]
x10 <- IA9_fig_vio$New.R2[IA9_fig_vio$tnx.1q=="HRD.1q.gain"]
x11 <- IA9_fig_vio$New.R2[IA9_fig_vio$tnx.1q=="NONE"]
x12 <- IA9_fig_vio$New.R2[IA9_fig_vio$tnx.1q=="NONE.1qgain"]

vioplot(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12, names=c("t11;14", "t11;14+1q","t4;14","t4;14+1q","t14;20", "t14;16","t6;14","Hrd","Hrd+11q","Hrd+1q","None","None+1q"), col="gold",  horizontal=FALSE)
abline(h=0.98,col=1,lty=2)
abline(h=0.97,col=2,lty=2)
title(expression("Violin plot of R^2 versus IgH translocation"), ylab="R^2", xlab="Translocation and copy number status")

####figure 1 part B ####
x10 <- IA9_CNV_BL$R2_amp[IA9_fig_vio$any.tnx=="t11_14" | IA9_fig_vio$any.tnx=="t12;14" | IA9_fig_vio$any.tnx=="t14_16" | IA9_fig_vio$any.tnx=="t14_20" | IA9_fig_vio$any.tnx=="t4_14" | IA9_fig_vio$any.tnx=="t6;14" | IA9_fig_vio$any.tnx=="t8_14"]
x11 <- IA9_CNV_BL$R2_amp[IA9_fig_vio$any.tnx=="none" & IA9_CNV_BL$Hyperdiploid_Call==1]
vioplot(x10,x11, names=c("IgH translocations", "No translocation with Hyperdipoidy"),ylim=c(0.60,1), col="gold",  horizontal=FALSE)
abline(h=0.98,col=1,lty=2)
title(expression("Violin plot of R^2 versus IgH translocation - CoMMpass IA9"), ylab="R^2", xlab="Translocation and copy number status")

x10 <- IA9_fig_vio$New.R2[IA9_fig_vio$any.tnx=="t11_14" | IA9_fig_vio$any.tnx=="t12;14" | IA9_fig_vio$any.tnx=="t14_16" | IA9_fig_vio$any.tnx=="t14_20" | IA9_fig_vio$any.tnx=="t4_14" | IA9_fig_vio$any.tnx=="t6;14" | IA9_fig_vio$any.tnx=="t8_14"]
x11 <- IA9_fig_vio$New.R2[IA9_fig_vio$any.tnx=="none" & IA9_CNV_BL$Hyperdiploid_Call==1]
vioplot(x10,x11, names=c("IgH translocations", "No translocation with Hyperdipoidy"),ylim=c(0.60,1), col="darkorange",  horizontal=FALSE)
abline(h=0.98,col=1,lty=2)
title(expression("Violin plot of R^2 versus IgH translocation - CoMMpass IA9"), ylab="R^2", xlab="Translocation and copy number status")


IA9_fig_vio$any.Tnx[IA9_fig_vio$Tnx=="t11_14" | IA9_fig_vio$Tnx=="t12;14" | IA9_fig_vio$Tnx=="t14_16" | IA9_fig_vio$Tnx=="t14_20" | IA9_fig_vio$Tnx=="t4_14" | IA9_fig_vio$Tnx=="t6;14" | IA9_fig_vio$Tnx=="t8_14"] <- 1
IA9_fig_vio$any.Tnx[IA9_fig_vio$Tnx=="none"] <- 0

IA9_fig_vioTnx <- subset(IA9_fig_vio, !is.na(New.R2))
IA9_fig_vioTnx <- subset(IA9_fig_vioTnx, !is.na(any.tnx))
prim <- factor(IA9_fig_vioTnx$any.tnx)
counts <- table(IA9_fig_vioTnx$R2Call, prim)

x15 <- IA9_fig_vio$New.R2[IA9_fig_vio$any.tnx==1]
x16 <- IA9_fig_vio$New.R2[IA9_fig_vio$any.tnx==0 & IA9_fig_vio$HRD=="HRD"]
vioplot(x15,x16, names=c("IgH translocations (24/109)", "No translocation with Hyperdiploidy (39/221)"), ylim=c(0.70,1), col="darkorange1",  horizontal=FALSE)
abline(h=0.98,col=1,lty=2)
title(expression("Violin plot of R^2 versus IgH translocation - MyXI", ylab="R^2"))

####figure 1 part B ####
perl.data.fig <- subset(perl.data, tnx.hrd==11 | tnx.hrd==4 | tnx.hrd==16 | tnx.hrd=="HRD")
x1 <- perl.data.fig$New.R2[perl.data.fig$tnx.1q.fig==11]
x2 <- perl.data.fig$New.R2[perl.data.fig$tnx.1q.fig=="11.1qgain"]
x3 <- perl.data.fig$New.R2[perl.data.fig$tnx.1q.fig==4]
x4 <- perl.data.fig$New.R2[perl.data.fig$tnx.1q.fig=="4.1qgain"]
x5 <- perl.data.fig$New.R2[perl.data.fig$tnx.1q.fig=="HRD"]
x6 <- perl.data.fig$New.R2[perl.data.fig$tnx.1q.fig=="HRD.1qgain"]
x7 <- perl.data.fig$New.R2[perl.data.fig$tnx.1q.fig==16]

vioplot(x3,x4,x7,x1,x2,x5,x6, names=c("t4;14 (3/21)", "t4;14+1q (6/22)", "t14;16 (4/14)","t11;14 (10/37)","t11;14+1q (2/13)","HRD (20/121)", "HRD+1q (13/62)"), col="gold",  horizontal=FALSE)
abline(h=0.98,col=1,lty=2)
title(expression("Violin plot of R^2 versus IgH translocation"), ylab="R^2")

####All tnxs ####
perl.data.fig <- subset(perl.data, tnx.hrd==11 | tnx.hrd==4 | tnx.hrd==16 | tnx.hrd=="HRD")x1 <- perl.data.fig$New.R2[perl.data.fig$tnx.1q.fig==11]
x1 <- IA9_fig_vio2$R2_amp[IA9_fig_vio2$Tnx.x=="t11_14"]
x2 <- IA9_fig_vio2$R2_amp[IA9_fig_vio2$Tnx.x=="t12_14"]
x3 <- IA9_fig_vio2$R2_amp[IA9_fig_vio2$Tnx.x=="none" & IA9_fig_vio2$HRD==1]
x5 <- IA9_fig_vio2$R2_amp[IA9_fig_vio2$Tnx.x=="t14_16"]
x6 <- IA9_fig_vio2$R2_amp[IA9_fig_vio2$Tnx.x=="t14_20"]
x7 <- IA9_fig_vio2$R2_amp[IA9_fig_vio2$Tnx.x=="t4_14"]
x8 <- IA9_fig_vio2$R2_amp[IA9_fig_vio2$Tnx.x=="t6_14"]
x9 <- IA9_fig_vio2$R2_amp[IA9_fig_vio2$Tnx.x=="t8_14"]
x10 <- IA9_fig_vio2$R2_amp[IA9_fig_vio2$Tnx.x=="t4_14_t14_16"]

vioplot(x8,x2,x7,x5,x6,x9,x1,x3, names=c("t6;14 (2/9)","t12;14 (3/8)","t4;14 (13/68)","t14;16 (7/21)","t14;20 (2/9)","t8;14 (0/6)","t11;14 (16/83)", "HRD (37/260)"), col="darkorange1",  horizontal=FALSE)
abline(h=0.98,col=1,lty=3)
title(expression("Violin plot of R^2 versus IgH translocation - CoMMpass IA9 - Amp - FNA"), ylab="R^2")

vioplot(x8,x2,x7,x5,x6,x9,x10,x1,x3, names=c("t6;14","t12;14","t4;14","t4;14+t14;16","t14;16","t14;20","t8;14","t11;14","HRD"), col="darkorange1",  horizontal=FALSE)
abline(h=0.98,col=1,lty=2)
title(expression("Violin plot of R^2 versus IgH translocation - CoMMpass IA9 - Amp - FNA"), ylab="R^2")

####All tnxs and 1q ####
x1 <- IA9_CNV_BL$R2_long[IA9_fig_vio$any.tnx=="t11_14" & IA9_CNV_BL$X1q21_20percent==1]
x2 <- IA9_CNV_BL$R2_long[IA9_fig_vio$any.tnx=="t11_14" & IA9_CNV_BL$X1q21_20percent==0]
x3 <- IA9_CNV_BL$R2_long[IA9_fig_vio$any.tnx=="none" & IA9_CNV_BL$Hyperdiploid_Call==1 & IA9_CNV_BL$X1q21_20percent==1]
x4 <- IA9_CNV_BL$R2_long[IA9_fig_vio$any.tnx=="none" & IA9_CNV_BL$Hyperdiploid_Call==1 & IA9_CNV_BL$X1q21_20percent==0]
x5 <- IA9_CNV_BL$R2_long[IA9_fig_vio$any.tnx=="t14_16" & IA9_CNV_BL$X1q21_20percent==1]
x6 <- IA9_CNV_BL$R2_long[IA9_fig_vio$any.tnx=="t14_16" & IA9_CNV_BL$X1q21_20percent==1]
x7 <- IA9_CNV_BL$R2_long[IA9_fig_vio$any.tnx=="t4_14" & IA9_CNV_BL$X1q21_20percent==1]
x8 <- IA9_CNV_BL$R2_long[IA9_fig_vio$any.tnx=="t4_14" & IA9_CNV_BL$X1q21_20percent==1]

vioplot(x7,x8,x5,x6,x1,x2,x3,x4, names=c("t6;14 (3/8)","t12;14 (2/7)","t4;14 (8/76)","t14;16 (4/23)","t14;20 (2/7)","t8;14 (0/7)","t11;14 (18/104)", "none (12/51)","HRD (31/258)"), col="gold",  horizontal=FALSE)
abline(h=0.98,col=1,lty=3)
title(expression("Violin plot of R^2 versus IgH translocation- CoMMpass IA9"), ylab="R^2")

vioplot(x7,x8,x5,x6,x1,x2,x3,x4, names=c("t4;14+1q","t4;14-1q","t14;16+1q","t14;16-1q","t11;14+1q","t11;14-1q","HRD+1q","HRD-1q"), col="gold",  horizontal=FALSE)
abline(h=0.98,col=1,lty=2)
title(expression("Violin plot of R^2 versus IgH translocation +/-1q21 - CoMMpass IA9"), ylab="R^2")

####figure 1 part A ####
IA9_CNV_R2 <- subset(IA9_CNV_R2, !is.na(R2_long))
IA9_CNV_BL <- subset(IA9_CNV_R2, rep.y==1)
IA9_CNV_BL <- subset(IA9_CNV_R2, !is.na(Tnx))


x1 <- IA9_CNV_BL$R2_long[IA9_fig_vio$any.tnx=="t11_14"]
x2 <- IA9_CNV_BL$R2_long[IA9_fig_vio$any.tnx=="t4_14"]
x3 <- IA9_CNV_BL$R2_long[IA9_fig_vio$any.tnx=="none" & IA9_CNV_BL$Hyperdiploid_Call==1]
x4 <- IA9_CNV_BL$R2_long[IA9_fig_vio$any.tnx=="t14_16"]

vioplot(x2,x4,x1,x3, names=c("t4;14 (8/72)", "t14;16 (4/23)","t11;14 (18/104)","HRD (31/258)"), col="gold", horizontal=FALSE)
abline(h=0.98,col=1,lty=2)
title(expression("Violin plot of R^2 versus IgH translocation - CoMMpass IA9"), ylab="R^2", xlab="Translocation and HRD")

prim <- factor(IA9_fig_vio$any.tnx, levels=levels(IA9_fig_vio$any.tnx)[c(3,1,2,4,5,6,7,8)])
counts <- table(IA9_CNV_BL$R2Ca_long, prim)
counts
prim
##t12;14 none t11_14 t14_16 t14_20 t4_14 t6;14 t8_14
##0      5  266     86     19      5    64     5     7
##1      2   43     18      4      2     8     3     0


IA9_CNV_BL$any.Tnx[IA9_fig_vio$any.tnx=="t11_14" | IA9_fig_vio$any.tnx=="t12;14" | IA9_fig_vio$any.tnx=="t14_16" | IA9_fig_vio$any.tnx=="t14_20" | IA9_fig_vio$any.tnx=="t4_14" | IA9_fig_vio$any.tnx=="t6;14" | IA9_fig_vio$any.tnx=="t8_14"] <- 1
IA9_CNV_BL$any.Tnx[IA9_fig_vio$any.tnx=="none"] <- 0

#test.wil.any.tnx <- wilcox.test(R2_amp ~ any.Tnx, data=IA9_CNV_BL)



