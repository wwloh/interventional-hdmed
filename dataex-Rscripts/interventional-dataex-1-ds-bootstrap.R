rm(list=ls())
libraries_check <- c("data.table","glmnet")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)
sessionInfo()

source("hdmed-dataprep.R")
source("funs-double_selection.R")

# double selection ############################################################
# initialize for parallel MC jobs
args <- 102
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
}
(seed <- as.integer(args[1]))
set.seed(NULL)

## bootstrapping double selection
fitA.OK <- FALSE
while (!fitA.OK) {
  boot_id <- sort(sample.int(nrow(mydata1),replace=TRUE)) # bootstrap sample
  Data <- mydata1[boot_id,list(id,y,age,gender,race,A)] # drop all mediators
  Data[,"id" := 1:nrow(Data)]
  setnames(Data,old="y",new="Y")
  Data_Monly <- mydata1_genes[boot_id,] # mediators only
  rm(mydata1,mydata1_genes)
    
  # check for convergence and violations of positivity assumption
  fitA <- glm(fitA.form, family = binomial("logit"), data = Data)
  fitA.OK <- fitA$converged & !fitA$boundary &
    min(min(fitA$fitted.values),1-max(fitA$fitted.values))>.Machine$double.eps*1e1
}

# generate model matrix of baseline confounders (without an intercept) for LASSO
Lmodelmat <- model.matrix(fitA)[,-1]
rm(fitA)
# save model matrix **with all mediators** for predicting potential outcomes
Data_M <- as.data.frame(cbind(
  "id"=Data$id,"(Intercept)"=1,Lmodelmat,"A"=Data[,A],Data_Monly))

res <- DoubleSelection(L=Lmodelmat,
                       A=Data[,A],
                       M=Data_Monly,
                       Y=Data[,Y],
                       fitY.family="binomial",
                       fitM.family="gaussian")
# about 100 mins for all mediators

if (exists("res")) {
  filename <- paste0("dataex-boot_ds-",seed,".Rdata")
  save(boot_id,res,file=filename)  
} else {
  q()
}


# estimate indirect effects via each mediator using selected predictors #######
res.ds <- NULL
ptm=proc.time()[3]
for (s in 1:length(res)) {
  # unique predictors selected in either mediator or outcome model
  dat.Ms.allnames <- unique(unlist(lapply(res[[s]], function(x)
    names(x[abs(x)>.Machine$double.eps]))))
  dat.Ms <- Data_M[,grep("(Intercept)", dat.Ms.allnames,invert=TRUE,value=TRUE),
                   drop=FALSE]
  dat.Ms <- cbind(dat.Ms, "Y"=Data$Y)
  # check for non-convergence or boundary solutions
  fitY <- glm(Y~., family = binomial("logit"), data = dat.Ms)
  fitY.OK <- fitY$converged & !fitY$boundary & (abs(min(fitY$coefficients))<1e6)
  rm(dat.Ms)
  
  if (fitY.OK==TRUE) {
    # save estimated coefficients
    fitY.coef <- coef(fitY)
    ## remove quotes around interactions
    names(fitY.coef) <- gsub(pattern="[`]",replacement="",names(fitY.coef))
    # same order of predictors as in double selection procedure
    fitY.coef.Ms <- fitY.coef[names(res[[s]]$Y)]
    rm(fitY.coef)
    names(fitY.coef.Ms) <- names(res[[s]]$Y)
    # fill in missing coefficients with zeroes
    fitY.coef.Ms[is.na(fitY.coef.Ms)] <- 0
  } else {
    # use (penalized) coefficients from outcome model directly
    fitY.coef.Ms <- res[[s]]$Y
  }
  
  res.ds <- c(res.ds, 
              OneMCestimator_noMmodels_oneIEonly(
                data=Data_M,
                mc_draws=10,
                fitY.coef=fitY.coef.Ms,
                Ms.i=paste0("M",s)))
  
  cat("M",s,":",res.ds[s],"\n")
  cat("time taken (mins) =", round((proc.time()[3]-ptm)/60), 
      "; approx. left (mins) =", round((proc.time()[3]-ptm)*(length(res)/s-1)/60),
      "\n")
}

save(boot_id,res,res.ds,file=filename)
q()
