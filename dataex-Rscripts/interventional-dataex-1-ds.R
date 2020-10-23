rm(list=ls())
libraries_check <- c("data.table","glmnet")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)
sessionInfo()

source("hdmed-dataprep.R")

# initialize for parallel MC jobs
args <- 1
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
  source("funs-hd-Yonly.R")
} else {
  source("../sims/funs-hd-Yonly.R")
}
(seed <- as.integer(args[1]))
set.seed(seed)

if (seed==1) {
  boot_id <- 1:nrow(mydata1)
} else {
  boot_id <- sort(sample.int(nrow(mydata1),replace=TRUE)) # bootstrap sample
}
if (identical(mydata1[boot_id,id],boot_id)) {
  Data <- mydata1[boot_id,list(id,y,age,gender,race,A)] # drop all mediators
  Data[,"id" := 1:nrow(Data)]
  setnames(Data,old="y",new="Y")
  Data_Monly <- mydata1_genes[boot_id,] # mediators only
  rm(mydata1,mydata1_genes)
    
  # check for convergence and violations of positivity assumption
  fitA <- glm(fitA.form, family = binomial("logit"), data = Data)
  fitA.OK <- fitA$converged & !fitA$boundary &
    min(min(fitA$fitted.values),1-max(fitA$fitted.values))>.Machine$double.eps*1e1
  if (!fitA.OK) q()
}

# generate model matrix of baseline confounders (without an intercept) for LASSO
Lmodelmat <- model.matrix(fitA)[,-1]
rm(fitA)
# save model matrix **with all mediators** for predicting potential outcomes
Data_M <- as.data.frame(cbind(
  "id"=Data$id,"(Intercept)"=1,Lmodelmat,"A"=Data[,A],Data_Monly))

# confounder selection ########################################################
# default lambda sequence in glmnet for cross-validation
# https://stackoverflow.com/questions/23686067/default-lambda-sequence-in-glmnet-for-cross-validation
## Standardize variables: (need to use n instead of (n-1) as denominator)
mysd <- function(y) sqrt(sum((y-mean(y))^2)/length(y))
x <- cbind(Lmodelmat,Data[,A],Data_Monly)
sx <- scale(x, scale = apply(x, 2, mysd))
sx <- as.matrix(sx, ncol = 20, nrow = 100)
sy <- as.vector(scale(Data[,Y], scale = mysd(Data[,Y])))

## Calculate lambda path (first get lambda_max):
lambda_max <- max(abs(colSums(sx*sy)))/nrow(Data)
epsilon <- .0001
K <- 100
lambdapath <- exp(seq(log(lambda_max), log(lambda_max*epsilon),length.out = K))
summary(lambdapath)
rm(mysd,x,sx,sy,epsilon,lambdapath)

lambda.seq <- exp(seq(log(lambda_max), log(ncol(Data_Monly)/nrow(Data)*1e-6),
                      length.out = K*2))
summary(lambda.seq)

filename <- paste0("dataex-boot_ds-",seed,".Rdata")
if (!file.exists(file=filename)) {
  # select confounders for outcome model ######################################
  res <- ConfounderSelection(L=Lmodelmat,
                             A=Data[,A],
                             M=Data_Monly,
                             Y=Data[,Y],
                             fitY.family="binomial",
                             addL=TRUE,
                             lambda.seq=lambda.seq)
  if (seed==1) {
    save(lambda.seq,res,file=filename)  
  } else {
    save(boot_id,Data_M,lambda.seq,res,file=filename)
  }
} else {
  load(file=filename)
}

n_mc <- 50
# estimate indirect effects via each mediator using selected predictors #####
res.lasso <- NULL
ptm=proc.time()[3]
for (s in 1:length(res)) {
  res.lasso[[s]] <- NA
  if (all(!is.na(res[[s]]))) {
    res.lasso[[s]] <- OneMCestimator_noMmodels_oneIEonly(
      data=Data_M,
      mc_draws=n_mc,
      fitY.coef=res[[s]],
      Ms.i=paste0("M",s))  
  }
  
  cat("M",s,":",res.lasso[[s]], "\n")
  cat("time taken (mins) =", round((proc.time()[3]-ptm)/60), 
      "; approx. left (mins) =", round((proc.time()[3]-ptm)*(length(res)/s-1)/60),
      "\n")
}

if (seed==1) {
  save(res,res.lasso,file=filename)  
} else {
  save(boot_id,res,res.lasso,file=filename)
}

q()
