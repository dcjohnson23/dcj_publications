##install.packages(c("survival","reshape2","plyr","survminer","rms","glmnet","caTools"))

library(survival)
library(reshape2)
library(plyr)

data.my9 <- read.table("/ph3_UK_all_GWAS_univariant_11_my9.csv",sep=",", header = TRUE, row.names = 1)
data.my11 <- read.table("/ph3_UK_all_GWAS_my11.csv",sep=",", header = TRUE, row.names = 1)
data.hdb <- read.table("/ph3_final.hdb.2014.Jan14.03.os.fish.csv",sep=",", header = TRUE, row.names = 1)
data.us <- read.table("/ph3.US_GWASData_demo.iss.PCA.csv",sep=",", header = TRUE, row.names = 1)

head(data.my9)
head(data.my11)
head(data.hdb)
head(data.us)

colnames(data.my9)
colnames(data.my11)
colnames(data.hdb)
colnames(data.us)

my9.vars <- c("chemothe", "age", "sex", "iss_com", "os_months", "os_status", "pfs_months", "pfs_status")
data.my9 <- data.my9[my9.vars]
my11.vars <- c("ind.treat", "Age", "sex", "iss_com", "os_months", "os_status", "pfs_months", "pfs_status")
data.my11 <- data.my11[my11.vars]
hdb.vars <- c("treatment", "age", "sex", "iss_com", "os_months", "os_status", "pfs_months", "pfs_status")
data.hdb <- data.hdb[hdb.vars]
us.vars <- c("Study", "Age", "Sex", "iss_com", "os_months", "os_status", "pfs_months", "pfs_status")
data.us <- data.us[us.vars]

data.my9$sex <-factor(data.my9$sex, levels = c("0", "1"), labels= c("f", "m"))
data.us$Sex <- factor(data.us$Sex, levels = c("0","1"), labels= c("f", "m"))

summary(data.my9)
summary(data.my11)
summary(data.hdb)
summary(data.us)

library("survminer")
#my9 curve#
cxmod_t9 <- coxph(Surv(os_months, os_status) ~ chemothe, data = data.my9)
cxmod_t9d <- survfit(Surv(os_months, os_status) ~ chemothe, data = data.my9)
newdat1 <- expand.grid(chemothe = levels(data.my9$chemothe))
rownames(newdat1) <- letters[1:7]
cxsf <- survfit(cxmod_t9, data = data.my9, newdata = newdat1, conf.type = "none")
surv_cxmod0 <- surv_summary(cxsf)
surv_cxmod <- cbind(surv_cxmod0, newdat1[as.character(surv_cxmod0$strata), ])
colnames(surv_cxmod)[10] <- "chemothe"
ggsurvplot_df(surv_cxmod, palette = "chemothe", legend.title = NULL, legend.labs = newdat1$chemothe,censor = FALSE, title = "My9 induction treatment")

#my11 curve#
cxmod_t11 <- coxph(Surv(os_months, os_status) ~ ind.treat, data = data.my11)
cxmod_t11d <- survfit(Surv(os_months, os_status) ~ ind.treat, data = data.my11)
newdat2 <- expand.grid(ind.treat = levels(data.my11$ind.treat))
rownames(newdat2) <- letters[1:4]
cxsf1 <- survfit(cxmod_t11, data = data.my11, newdata = newdat2, conf.type = "none")
surv_cxmod10 <- surv_summary(cxsf1)
surv_cxmod1 <- cbind(surv_cxmod10, newdat2[as.character(surv_cxmod10$strata), ])
colnames(surv_cxmod1)[10] <- "ind.treat"
ggsurvplot_df(surv_cxmod1, palette = "ind.treat", legend.title = NULL, legend.labs = newdat2$ind.treat ,censor = FALSE, title = "My11 induction treatment")

#hdb curve#
cxmod_hdb <- coxph(Surv(os_months, os_status) ~ treatment, data = data.hdb)
cxmod_hdbd <- survfit(Surv(os_months, os_status) ~ treatment, data = data.hdb)
newdat3 <- expand.grid(treatment = levels(data.hdb$treatment))
rownames(newdat2) <- letters[1:4]
cxsf2 <- survfit(cxmod_hdb, data = data.hdb, newdata = newdat3, conf.type = "none")
surv_cxmod20 <- surv_summary(cxsf2)
surv_cxmod2 <- cbind(surv_cxmod20, newdat3[as.character(surv_cxmod20$strata), ])
colnames(surv_cxmod2)[10] <- "treatment"
ggsurvplot_df(surv_cxmod2, palette = "treatment", legend.title = NULL, legend.labs = newdat3$treatment, censor = FALSE, title = "Hdb induction treatment")

#US curve#
cxmod_us <- coxph(Surv(os_months, os_status) ~ Study, data = data.us)
cxmod_usd <- survfit(Surv(os_months, os_status) ~ Study, data = data.us)

newdat4 <- expand.grid(Study = levels(data.us$Study))
rownames(newdat4) <- letters[1:5]
cxsf3 <- survfit(cxmod_us, data = data.us, newdata = newdat4, conf.type = "none")
surv_cxmod30 <- surv_summary(cxsf3)
surv_cxmod3 <- cbind(surv_cxmod30, newdat4[as.character(surv_cxmod30$strata), ])
colnames(surv_cxmod3)[10] <- "Study"
ggsurvplot_df(surv_cxmod3, palette = "Study", legend.title = NULL, legend.labs = newdat4$Study, censor = FALSE, title =  "US induction treatment")

##########################
###Medians################
##########################
med9d <- surv_median(cxmod_t9d)
med11d <- surv_median(cxmod_t11d)
medhdb <- surv_median(cxmod_hdbd)
medus <- surv_median(cxmod_usd)

trialmedians <- rbind.data.frame(med9d, med11d, medhdb, medus)

##########################
### Estimated Medians#####
##########################
cxmod_t9r <- survreg(Surv(os_months, os_status) ~ chemothe, data = data.my9)
cxmod_t11r <- survreg(Surv(os_months, os_status) ~ ind.treat, data = data.my11)
cxmod_hdbr <- survreg(Surv(os_months, os_status) ~ treatment, data = data.hdb)
cxmod_usr <- survreg(Surv(os_months, os_status) ~ Study, data = data.us)
summary(cxmod_t9r)
summary(cxmod_t11r)
summary(cxmod_hdbr)
summary(cxmod_usr)
#my9#
Est.med.my9.CTD <-predict(cxmod_t9r, newdata=data.frame(chemothe="CTD"), type="quantile", p=0.5)
Est.med.my9.CTDa <-predict(cxmod_t9r, newdata=data.frame(chemothe="CTDa"), type="quantile", p=0.5)
Est.med.my9.CVAD <-predict(cxmod_t9r, newdata=data.frame(chemothe="CVAD"), type="quantile", p=0.5)
Est.med.my9.MP <-predict(cxmod_t9r, newdata=data.frame(chemothe="MP"), type="quantile", p=0.5)
Est.med.my9.Noin <-predict(cxmod_t9r, newdata=data.frame(chemothe="No_induction"), type="quantile", p=0.5)
#my11#
Est.med.my11.CTD <-predict(cxmod_t11r, newdata=data.frame(ind.treat="CTD"), type="quantile", p=0.5)
Est.med.my11.CTDa <-predict(cxmod_t11r, newdata=data.frame(ind.treat="CTDa"), type="quantile", p=0.5)
Est.med.my11.RCD <-predict(cxmod_t11r, newdata=data.frame(ind.treat="RCD"), type="quantile", p=0.5)
Est.med.my11.RCDa <-predict(cxmod_t11r, newdata=data.frame(ind.treat="RCDa"), type="quantile", p=0.5)
#hdb#
Est.med.hdb.Other <-predict(cxmod_hdbr, newdata=data.frame(treatment="OTHER"), type="quantile", p=0.5)
Est.med.hdb.VAD <-predict(cxmod_hdbr, newdata=data.frame(treatment="VAD"), type="quantile", p=0.5)
Est.med.hdb.PAD <-predict(cxmod_hdbr, newdata=data.frame(treatment="PAD"), type="quantile", p=0.5)
Est.med.hdb.TAD <-predict(cxmod_hdbr, newdata=data.frame(treatment="TAD"), type="quantile", p=0.5)
#us#
Est.med.us.TT2 <-predict(cxmod_usr, newdata=data.frame(Study="TT2"), type="quantile", p=0.5)
Est.med.us.TT3 <-predict(cxmod_usr, newdata=data.frame(Study="TT3"), type="quantile", p=0.5)
Est.med.us.TT3B <-predict(cxmod_usr, newdata=data.frame(Study="TT3B"), type="quantile", p=0.5)
Est.med.us.TT4 <-predict(cxmod_usr, newdata=data.frame(Study="TT4"), type="quantile", p=0.5)
Est.med.us.TTnt <-predict(cxmod_usr, newdata=data.frame(Study="TTnt"), type="quantile", p=0.5)
trialmedians$est.medians <- c(71.1,35.1,72.9,31.1,12,"NA","NA", 81.9,38.1,102,45.8,93.1,68.4,280,86.2,110,166,143,136,69)

###############################################################
###Covariant selection - backward elimination on covariates####
###Fast Backward Variable Selection###########################
##############################################################

library(rms)
##my9##
my.survival.object=Surv(data.my9$os_months, data.my9$os_status)
fit=cph(formula=my.survival.object~chemothe + age + sex + iss_com
, data=data.my9, x=TRUE, y=TRUE)
validate(fit, method="boot", B=400, bw=T, rule="aic",
         type="residual", sls=.1, aics=0, force=NULL, estimates=FALSE,
         pr=FALSE, dxy=TRUE, u, tol=1e-9)
##my11##
my.survival.object11=Surv(data.my11$os_months, data.my11$os_status)
fit11=cph(formula=my.survival.object11~ind.treat+ Age + sex + iss_com
, data=data.my11, x=TRUE, y=TRUE)
validate(fit11, method="boot", B=400, bw=T, rule="aic",
         type="residual", sls=.1, aics=0, force=NULL, estimates=FALSE,
         pr=FALSE, dxy=TRUE, u, tol=1e-9)
##hdb##
my.survival.objectH=Surv(data.hdb$os_months, data.hdb$os_status)
fitH=cph(formula=my.survival.objectH~treatment + age + sex + iss_com
, data=data.hdb, x=TRUE, y=TRUE)
validate(fitH, method="boot", B=400, bw=T, rule="aic",
         type="residual", sls=.1, aics=0, force=NULL, estimates=FALSE,
         pr=FALSE, dxy=TRUE, u, tol=1e-9)
##us##
my.survival.objectU=Surv(data.us$os_months, data.us$os_status)
fitU=cph(formula=my.survival.objectU~Study + Age + Sex + iss_com
, data=data.us, x=TRUE, y=TRUE)
validate(fitU, method="boot", B=400, bw=T, rule="aic",
         type="residual", sls=.1, aics=0, force=NULL, estimates=FALSE,
         pr=FALSE, dxy=TRUE, u, tol=1e-9)

coef(fit)
coef(fit11)
coef(fitH)
coef(fitU)


############
##Lasso#####
############



#create complete dataset#
library(glmnet)
#library(XGBoost)
library(dygraphs)
library(coxnet)
data2.my11 <- read.table("C:/Users/Utilisateur/Desktop/CV and references/Informa/ph3_UK_all_GWAS_my11.csv",sep=",", header = TRUE, row.names = 1)
data2.my11.X <- data2.my11[c(-8,-10,-11,-19,-24,-28,-29,-32)]
data2.my11.X <- data2.my11.X[complete.cases(data2.my11.X), ]
data2.my11.Xc <- data2.my11.X[c(-9,-10)]


#spilt off survival columns from potential covariant columns#
condemosCol <-
survCols <- 
dim(demosCol)

#boxplot for each attribute 
par(mfrow=c(1, 8))
	for(i in 1:4) {
	boxplot(x[,i], main=names(iris)[i])
	}

#scatterplot matrix
featurePlot(x=x, y=y, plot=Ellei


data2.my11.Xfactors <- model.matrix(~ sex + iss_com + WHO + pathway + ind.treat + BD.y1_n0 + ptype + lightc + Tx + HRD + CDKN2C_homo_del + CCND1_gain_bothprobes_morthan_11q25 + MYC_gainoramp_eitherprobe + Gain11q25 + Del1p32_3 + Gain_or_Amp_1q_CKS1B_binary + Del13q + Del17p_P425, data2.my11.Xc)
data2.my11.Xcc <- as.matrix(data.frame(data2.my11.Xc$Age, data2.my11.Xc$b2mloc, data2.my11.Xc$screat,data2.my11.Xfactors))

cv.fit = cv.glmnet(data2.my11.Xcc, Surv(data2.my11.X$os_months, data2.my11.X$os_status),
                   family="cox", alpha=1, nlambda= 100)
plot(cv.fit)
plot(cv.fit$glmnet.fit, xvar="lambda", label = TRUE)
coef(cv.fit, s="lambda.min")
cv.fit$lambda.min
gl.fit = glmnet(data2.my11.Xcc, Surv(data2.my11.X$os_months, data2.my11.X$os_status),
                   family="cox", alpha=1, lambda = cv.fit$lambda.1se)
gl.fit$beta[,1]
coef(gl.fit, s="lambda.min")


x_train <- build.x(modelFormula, data=data, contrasts=FALSE, sparse=TRUE)
y_train <- build.y(modelFormula, data=data) %>% as.interger()- 1

x_val <- build.x(modelFormula, data=data, contrasts=FALSE, sparse=TRUE)
y_val <- build.y(modelFormula, data=data) %>% as.interger()- 1

mod_glmnet <- cv.glmnet(x=x_train, y=y_train, family="binomial", nfolds=10)
corefpath(mod_glmne) 
coefplot(mod_glmnet,lambda='lambda.min', sort='magnitude')





##############
###boosting###
##############

library("gbm")
my.survival.objectB=Surv(data.us$os_months, data.us$os_status)
cox.glm=gbm(my.survival.objectB ~ sex + iss_com + WHO + pathway + ind.treat + BD.y1_n0 + ptype + lightc + Tx + HRD + CDKN2C_homo_del, data=data2.my11.Xc)

library("BhGLM")



#split for train_test##

set.seed(2256)
#ind <- sample(2,nrow(data2.my11.X), replace= T, prob= c(0.7,0.3))
#train <- data2.my11.X[ind==1,]
#test <- data2.my11.X[ind==2,]
library(caTools)
sample.split(data2.my11.X$os_status, SplitRatio = 0.7) -> split_values
train1 <- subset(data2.my11.X, split_values==T)
test1 <- subset(data2.my11.X, split_values==F)





