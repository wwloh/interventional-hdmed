rm(list=ls())
libraries_check <- c("data.table","glmnet")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)
sessionInfo()

# simulation settings
simsettings <- expand.grid("t"=c(1,2)*50, # number of possible mediators
                           "n"=c(1,2)*50, # sample size
                           "bsig"=c(0,1,5), # coef of non-mediators
                           "ybin"=0:1 # binary or continuous outcome
                           )
simsettings <- simsettings[simsettings$t==simsettings$n,]
simsettings <- simsettings[-c(1:4),] # completed jobs
row.names(simsettings) <- NULL
simsettings

# initialize for parallel cluster jobs
args <- 6
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
  # additional seed is to calculate true population effects
  simsettings <- simsettings[rep(1:nrow(simsettings),each=1001),]
  nrow(simsettings)
}
(seed <- as.integer(args[1]))
rm(args)

filename <- paste(c(rbind(names(simsettings),
                          unlist(simsettings[seed,]))),collapse="_")
(filename <- gsub(pattern="[.]",replacement="",filename))
(t <- as.integer(simsettings[seed,"t"]))
(n <- as.integer(simsettings[seed,"n"]))
(bsig <- simsettings[seed,"bsig"])
(y_binary <- as.logical(simsettings[seed,"ybin"]))
N <- 1e4 # population size
n_mc <- 100 # number of MC draws

set.seed(9000) # fix seed for population data
source("interventional-randomA-sim3-gen_data.R")
# source("funs-hd-Yonly.R")
source("funs-hd-Yonly-morelambdas.R") # try other options for binary Y

fitA.form <- as.formula("A~1") # randomly assigned treatment

# population value of indirect effects ========================================
fitY.form <- as.formula(paste0("Y~A+L+",paste("M",1:t,sep="",collapse="+")))
if (y_binary) {
  fitY.pop <- glm(formula=fitY.form,family=binomial("logit"),data=popdat)  
} else {
  fitY.pop <- lm(formula=fitY.form,data=popdat)  
}
if (!grepl("[.]",row.names(simsettings)[seed])) {
  popeff <- list()
  ptm=proc.time()[3]
  for (s in 1:t) {
    popeff[[s]] <- OneMCestimator_noMmodels_oneIEonly(
      data=popdat, 
      mc_draws=ifelse(t<=100,n_mc,10),
      fitY.coef=fitY.pop$coefficients,
      Ms.i=paste0("M",s))
    
    cat("M",s,":",popeff[[s]],"\n")
    cat("time taken (mins) =", round((proc.time()[3]-ptm)/60,1), 
        "; approx. left (mins) =", round((proc.time()[3]-ptm)*(t/s-1)/60),
        "\n")
  }
  # 70 mins for 50 mediators and 1e4 observations
  pop.setting <- simsettings[seed,]
  save(pop.setting,popdat,popeff,fitY.pop,
       file=paste0("sim-hd-popeffs-",filename,".Rdata"))
  q()
}
rm(fitY.form)

# estimates for a given sample ================================================
OneEst <- function(Data, verbose=FALSE) {
  
  # select confounders for outcome model ######################################
  res <- ConfounderSelection(L=Data[,"L",drop=FALSE],
                             A=Data[,"A"],
                             M=Data[,grep(pattern="M",names(Data),value=TRUE)],
                             Y=Data[,"Y"],
                             fitY.family=ifelse(y_binary,"binomial","gaussian"),
                             addL=TRUE)
  
  # estimate indirect effects via each mediator using selected predictors #####
  res.lasso <- NULL
  ptm=proc.time()[3]
  for (s in 1:length(res)) {
    res.lasso[[s]] <- NA
    if (all(!is.na(res[[s]]))) {
      res.lasso[[s]] <- OneMCestimator_noMmodels_oneIEonly(
        data=Data,
        mc_draws=ifelse(n<=100,n_mc,10),
        fitY.coef=res[[s]],
        Ms.i=paste0("M",s))  
    }
    
    if (verbose) {
      cat("M",s,":",res.lasso[[s]], "\n")
      cat("time taken (mins) =", round((proc.time()[3]-ptm)/60), 
          "; approx. left (mins) =", round((proc.time()[3]-ptm)*(t/s-1)/60),
          "\n")
    }
  }
  return( unlist(list("est"=res.lasso)) )
}

set.seed(seed=NULL) # re-initializes
(nsims <- ceiling(100/t))
stm <- proc.time()[3]
simres <- lapply(1:nsims, function(sim) {
  # sample observations
  Data <- popdat[sample(nrow(popdat),n,replace=FALSE),]
  Data$id <- 1:nrow(Data)
  row.names(Data) <- NULL
  
  # point estimates
  res.pt <- OneEst(Data,verbose=TRUE)
  
  if (t<=100) {
    # ordinary bootstrap
    n.boots <- n_mc
    res.bca <- lapply(1:n.boots, function(i) {
      boot_idx <- sort(unique(sample(n,replace=TRUE)))
      boot.stm <- tryCatch(system.time(
        res.boot <- OneEst(Data[boot_idx,])
      )[3], error=function(cond) return(NA))
      if (all(!is.na(boot.stm))) {
        return(res.boot)
      } else {
        return(NULL)
      }
    })
    res.bca <- do.call(rbind,res.bca[!sapply(res.bca, is.null)])
  } else {
    res.bca <- NULL
  }
  
  res.popYcoef <- NULL
  for (s in 1:t) {
    # true values of coefficients
    res.popYcoef[[s]] <- OneMCestimator_noMmodels_oneIEonly(
      data=Data,
      mc_draws=ifelse(n<=100,n_mc,10),
      fitY.coef=fitY.pop$coefficients,
      Ms.i=paste0("M",s))
  }
  res.pt <- c(res.pt,unlist(list("true"=res.popYcoef)))
  
  return( list("settings"=simsettings[seed,],"res"=res.pt,"boot"=res.bca) )
})
proc.time()[3]-stm
# 100 min per sim for 50 mediators with n=50
# simres <- do.call(rbind,lapply(simres, unlist))

save(simres,file=paste0("interventional-randomA-sim3-hd-ds-",
                        filename,"-",seed,".Rdata"))
q()
