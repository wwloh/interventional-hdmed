rm(list=ls())
libraries_check <- c("data.table","glmnet","ranger")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)
sessionInfo()

## download data from supplemental materials of Biometrics paper:
## https://doi.org/10.1111/biom.12421
load("hdmed/biom12421-sup-0002-suppdata-code/Dataset_for_Biometrics.RData")
# id to preserve initial ordering of data
mydata1 <- data.table(cbind("id"=1:nrow(mydata1), mydata1))
setkey(mydata1)
nrow(mydata1)

# exposure for analysis ########################################################
# dichotomize exposure
mydata1[, A := as.integer(hsa.miR.223 > median(hsa.miR.223))]

# outcome for analysis ########################################################
# loss to follow up (y = NA)
mydata1[is.na(y),list(y,fut,death)]
# compare between the treatment groups
mydata1[is.na(y),table(is.na(y),A)]
# recode as 0
mydata1[is.na(y),y := 0]

# selected mediators ##########################################################
## glmnet
#### check for NAs
colnames(mydata1)[1:12] # check column names for mediators
mydata1_genes <- mydata1[,10:ncol(mydata1)] # mediators are in columns >= 10
dim(mydata1_genes)
genes.noNA <- which(apply(mydata1_genes,2,function(x) all(!is.na(x))))
mydata1_genes <- mydata1_genes[,genes.noNA,with=FALSE]
dim(mydata1_genes)
#### standardize gene expression values
mydata1_genes <- scale(mydata1_genes)
#### include confounders
C <- mydata1[,list(
  age,gender,race,
  "age_gender"=age*gender,
  "age_race"=age*race,
  "gender_race"=gender*race,
  "age_gender_race"=age*gender*race,A)]
#### penalty factors to always include confounders in LASSO
penfac <- c(rep(0,ncol(C)),rep(1,ncol(mydata1_genes)))
mydata1_genes <- cbind(C,mydata1_genes)
dim(mydata1_genes)

fit <- glmnet(x=as.matrix(mydata1_genes), y=mydata1$y, 
              family = "binomial", penalty.factor=penfac)
# plot(fit)
dim(coef(fit))
fit$df
min(fit$df)==ncol(C) # smallest fitted model has all confounders + treatment
# smallest model with at least ten mediators
(coef_names <- rownames(coef(fit))[
  which(coef(fit)[,which.max(fit$df>(fit$df[1]+10))]!=0)])
(gene_names <- coef_names[(ncol(C)+1+1):length(coef_names)])

# dataset for analysis ########################################################
Lnames <- c("age","gender","race")
data <- mydata1[, c("id",Lnames,"A",gene_names,"y"), with=FALSE]
setkey(data)
setnames(data, 
         old=c(gene_names,"y"), 
         new=c(paste0("M",1:length(gene_names)),"Y"))
setkey(data)

# outcome model
fitY.form <- as.formula(
  paste0("Y~A+",paste(Lnames,collapse="*"),"+",
         paste("M",1:length(gene_names),sep="",collapse="+")))
# propensity score model
fitA.form <- as.formula(paste0("A~",paste(Lnames,collapse="*")))

data <- as.data.frame(data)

# check for convergence and violations of positivity assumption
fitY <- glm(fitY.form, family = binomial("logit"), data = data)
fitY$converged
fitA <- glm(fitA.form, family = binomial("logit"), data = data)
min(min(fitA$fitted.values),1-max(fitA$fitted.values))>.Machine$double.eps*1e1
rm(fitY,fitA)

source("funs-dataex.R")

n_mc <- 1e3
n_boots <- 2e2

ptm <- proc.time()[3]
res_obs <- OneMCestimator_noMmodels(data=data,mc_draws=n_mc)
proc.time()[3] - ptm
# 2 mins with 1000 MC

# relabel observed data
obsdata <- data; rm(data)

# initialize for parallel MC jobs
args <- 1
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
}
(seed <- as.integer(args[1]))
set.seed(seed*n_boots)

# bootstrap 
ptm <- proc.time()[3]
res_boot <- list()
for (bb in 1:n_boots) {
  # resample <- TRUE
  # while(resample) {
    boot_data <- obsdata[sample.int(nrow(obsdata),replace=TRUE),]
    # check for convergence and violations of positivity assumption
    fitY <- glm(fitY.form, family = binomial("logit"), data = boot_data)
    fitA <- glm(fitA.form, family = binomial("logit"), data = boot_data)
    resample <- !(fitY$converged &
                    min(min(fitA$fitted.values),1-max(fitA$fitted.values))>
                    .Machine$double.eps*1e1)
    rm(fitY,fitA)
  # }
  boot_data[,"id"] <- 1:nrow(boot_data)
  boot_res <- OneMCestimator_noMmodels(data=boot_data,mc_draws=n_mc)
  res_boot[[bb]] <- c("resample"=resample,"boot"=bb,unlist(boot_res))
  # (1:nrow(data)) %in% sort(unique(boot_ids))
}
proc.time()[3] - ptm
res_boot.dt <- data.table(do.call(rbind,res_boot))
setkey(res_boot.dt)

save(res_boot,file=paste0("dataex-bootres-",seed,".Rdata"))
q()

# combine results =============================================================
setwd("dataex/")
myfiles <- list.files()
myfiles <- myfiles[grep(pattern=".Rdata",myfiles)]
res_boot_df <- NULL
for (ll in myfiles) {
  load(ll)
  res_boot_df <- rbind(res_boot_df,cbind("boot"=1:length(res_boot),
                                         do.call(rbind,res_boot)))
  rm(res_boot)
  cat(which(myfiles %in% ll),"\n")
}
resample_boot_df <- res_boot_df[,grep("resample",colnames(res_boot_df))]
wts_boot_df <- res_boot_df[,grep("we.raw",colnames(res_boot_df))]
res_boot_df <- res_boot_df[,grep("resample",colnames(res_boot_df),invert=TRUE)]
res_boot_df <- res_boot_df[,grep("we",colnames(res_boot_df),invert=TRUE)]
## check that bootstrap samples (across parallel jobs) are unique
by(res_boot_df,res_boot_df[,"boot"], function(x) nrow(unique(x))) # must be >1!

res_boot_df <- res_boot_df[,grep("boot",colnames(res_boot_df),invert=TRUE)]
res_boot_df <- cbind(res_boot_df,"sum"=rowSums(
  res_boot_df[,grep("theta",colnames(res_boot_df))]))

res_boot.dt <- data.table(res_boot_df)
setkey(res_boot.dt)

res_obs <- unlist(res_obs)
wts_obs <- res_obs[grep("we",names(res_obs))]
res_obs <- res_obs[grep("we",names(res_obs),invert=TRUE)]
res_obs <- c(res_obs,"sum"=sum(res_obs[grep("theta",names(res_obs))]))

res <- rbind("est"=res_obs,
             "se"=apply(res_boot.dt, 2, sd),
             apply(res_boot.dt, 2, quantile, probs=c(.025,.975)))
# direct, joint/mutual, significant indirect effects
res_summary <- res[,unique(c(grep("theta",colnames(res),invert=TRUE),
                             grep("theta",colnames(res))))]
                             # which(sign(res[3,]) == sign(res[4,]))))]
M.signif <- lapply(strsplit(colnames(res_summary),"[.]t"),function(x) x[2])
colnames(res_summary) <- c("Direct effect",
                           "Joint indirect effect ",
                           "Indirect effect via mutual dependence",
                           "Sum of mediator-specific indirect effects",
                           paste0("quad ",
                                  gene_names[na.omit(as.integer(M.signif))]))
res_summary <- t(res_summary)

library("xtable")
xtable(res_summary)

xtable(cbind(round(res_summary[,1:2],2),apply(res_summary, 1, function(x) {
  x_ <- round(x,2)
  paste0("(",x_[3],", ",x_[4],")")
})))

# bootstrap sampling dist
unique(res_boot.dt)
pdf("plot-dataex-deie.pdf",width=9,height=3)
par(mfrow=c(1,3))
plot(res_boot.dt$gamma.g1,res_boot.dt$gamma.g0, pch=".",
     xlab="Joint indirect effect",
     ylab="Direct effect")
abline(v=res_obs["gamma.g1"],h=res_obs["gamma.g0"])
plot(res_boot.dt$gamma.g1,res_boot.dt$sum, pch=".",
     xlab="Joint indirect effect",
     ylab="Sum of mediator-specific indirect effects")
abline(v=res_obs["gamma.g1"],h=res_obs["sum"])
plot(res_boot.dt$gamma.g1,res_boot.dt$mu, pch=".",
     xlab="Joint indirect effect",
     ylab="Indirect effect via mutual dependence")
abline(v=res_obs["gamma.g1"],h=res_obs["mu"])
dev.off()

# weights
pdf("plot-dataex-weights-minmax.pdf",width=6,height=6)
par(mfrow=c(1,1))
plot(wts_boot_df[,c(1,5)], pch=".", 
     xlim=range(wts_boot_df[,1]),ylim=range(wts_boot_df[,5]),
     main="Reciprocal of predicted probability of observed treatment",
     xlab="Minimum",ylab="Maximum")
abline(v=wts_obs[1],h=wts_obs[5])
abline(a=0,b=1,lty=2)
dev.off()