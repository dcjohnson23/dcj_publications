library(survival)
library(survminer)

#Intensive OS
perl.surv<- Surv(OS.months, OS.status) ~ cg00216138
Sperl.diff <- survdiff(perl.surv, data=My11.demo.annt.int.cpg)
Sperl.fit <- surv_fit(perl.surv, data=My11.demo.annt.int.cpg)
perl.surv<- Surv(OS.months, OS.status) ~ cg00216138
perl.diff <- survdiff(perl.surv, data=My11.demo.annt.int.cpg, na.action=na.omit)
perl.fit <- survfit(perl.surv, data=My11.demo.annt.int.cpg, na.action=na.omit)
sthis.sf1 <- smed(perl.fit)
cox.modelOSint<-coxph(Surv(OS.months,OS.status)~cg00216138+rand+Sex+age.x+k10status,data=My11.demo.annt.int.cpg)

tiff(file = "Surminer_methOS.int.nov19.tif",width = 11.7, height=8.3,unit="in", res=600)
ggsurv <- ggsurvplot(Sperl.fit, data = My11.demo.annt.int.cpg, conf.int = FALSE, pval = FALSE, fun = "pct", xlim = c(0,70),  surv.median.line = "hv", risk.table = TRUE, size = 1, linetype = "strata", palette = c("red","black"),legend = "bottom", legend.title = "DNA methylation - cg00216138", legend.labs = c("High","Low"), break.time.by = 3)
ggsurv$plot <- ggsurv$plot+ ggplot2::annotate("text", 
                                              x = 12, y = 25, # x and y coordinates of the text
                                              label = "HR= 0.28 (0.18-0.43) \n Logrank P < 4.0E-13 \n Coxph P < 1.73E-08 \n No of Patients: High = 66, Low = 299 \n Medians: High = 36, Low = 64 months", size = 6)
ggsurv
dev.off()

#Intensive PFS
perl.surv<- Surv(PFS.months, PFS.status) ~ cg00216138
Sperl.diff <- survdiff(perl.surv, data=My11.demo.annt.int.cpg)
Sperl.fit <- surv_fit(perl.surv, data=My11.demo.annt.int.cpg)
perl.surv<- Surv(PFS.months, PFS.status) ~ cg00216138
perl.diff <- survdiff(perl.surv, data=My11.demo.annt.int.cpg, na.action=na.omit)
perl.fit <- survfit(perl.surv, data=My11.demo.annt.int.cpg, na.action=na.omit)
sthis.sf1 <- smed(perl.fit)
cox.modelPFSint<-coxph(Surv(PFS.months,PFS.status)~cg00216138+rand+Sex+age.x+k10status,data=My11.demo.annt.int.cpg)

tiff(file = "Surminer_methPFS.int.nov19.tif",width = 11.7, height=8.3,unit="in", res=600)
ggsurv <- ggsurvplot(Sperl.fit, data = My11.demo.annt.int.cpg, conf.int = FALSE, pval = FALSE, fun = "pct", xlim = c(0,70),  surv.median.line = "hv", risk.table = TRUE, size = 1, linetype = "strata", palette = c("red","black"),legend = "bottom", legend.title = "DNA methylation - cg00216138", legend.labs = c("High","Low"), break.time.by = 3)
ggsurv$plot <- ggsurv$plot+ ggplot2::annotate("text", 
                                              x = 12, y = 25, # x and y coordinates of the text
                                              label = "HR= 0.38 (0.27-0.53) \n Logrank P < 8.46E-10 \n Coxph P < 3.6E-08 \n No of Patients: High = 66, Low = 299 \n Medians: High = 16.7, Low = 34.3 months", size = 6)
ggsurv
dev.off()

#Non-Intensive OS
perl.surv<- Surv(OS.months, OS.status) ~ cg00216138
Sperl.diff <- survdiff(perl.surv, data=My11.demo.annt.noint.cpg)
Sperl.fit <- surv_fit(perl.surv, data=My11.demo.annt.noint.cpg)
perl.surv<- Surv(OS.months, OS.status) ~ cg00216138
perl.diff <- survdiff(perl.surv, data=My11.demo.annt.noint.cpg, na.action=na.omit)
perl.fit <- survfit(perl.surv, data=My11.demo.annt.noint.cpg, na.action=na.omit)
sthis.sf1 <- smed(perl.fit)
cox.modelOSnoint<-coxph(Surv(OS.months,OS.status)~cg00216138+rand+Sex+age.x+k10status,data=My11.demo.annt.noint.cpg)

tiff(file = "Surminer_methOS.noint.nov19.tif",width = 11.7, height=8.3,unit="in", res=600)
ggsurv <- ggsurvplot(Sperl.fit, data = My11.demo.annt.noint.cpg, conf.int = FALSE, pval = FALSE, fun = "pct", xlim = c(0,64),  surv.median.line = "hv", risk.table = TRUE, size = 1, linetype = "strata", palette = c("red","black"),legend = "bottom", legend.title = "DNA methylation - cg00216138", legend.labs = c("High","Low"), break.time.by = 3)
ggsurv$plot <- ggsurv$plot+ ggplot2::annotate("text", 
                                              x = 12, y = 25, # x and y coordinates of the text
                                              label = "HR= 0.48 (0.25-0.92) \n Logrank P < 0.11 \n Coxph P < 0.03 \n No of Patients: High = 23, Low = 192 \n Medians: High = 33.1, Low = 43.3 months", size = 6)
ggsurv
dev.off()

#Non-Intensive PFS
perl.surv<- Surv(PFS.months, PFS.status) ~ cg00216138
Sperl.diff <- survdiff(perl.surv, data=My11.demo.annt.noint.cpg)
Sperl.fit <- surv_fit(perl.surv, data=My11.demo.annt.noint.cpg)
perl.surv<- Surv(PFS.months, PFS.status) ~ cg00216138
perl.diff <- survdiff(perl.surv, data=My11.demo.annt.noint.cpg, na.action=na.omit)
perl.fit <- survfit(perl.surv, data=My11.demo.annt.noint.cpg, na.action=na.omit)
sthis.sf1 <- smed(perl.fit)
cox.modelPFSnoint<-coxph(Surv(PFS.months,PFS.status)~cg00216138+rand+Sex+age.x+k10status,data=My11.demo.annt.noint.cpg)

tiff(file = "Surminer_methPFS.noint.nov19.tif",width = 11.7, height=8.3,unit="in", res=600)
ggsurv <- ggsurvplot(Sperl.fit, data = My11.demo.annt.noint.cpg, conf.int = FALSE, pval = FALSE, fun = "pct", xlim = c(0,60),  surv.median.line = "hv", risk.table = TRUE, size = 1, linetype = "strata", palette = c("red","black"),legend = "bottom", legend.title = "DNA methylation - cg00216138", legend.labs = c("High","Low"), break.time.by = 3)
ggsurv$plot <- ggsurv$plot+ ggplot2::annotate("text", 
                                              x = 46, y = 75, # x and y coordinates of the text
                                              label = "HR= 0.38 (0.22-0.66) \n Logrank P < 0.006 \n Coxph P < 6E-04 \n No of Patients: High = 23 : Low = 192 \n Medians: High = 10.6 : Low = 16.2 months", size = 6)
ggsurv
dev.off()






