rm(list=ls())
libraries_check <- c("data.table")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)
sessionInfo()

# initialize for parallel cluster jobs
args <- 4
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
}
(seed <- as.integer(args[1]))
rm(args)
set.seed(seed+9000)

# simulation settings
simsettings <- expand.grid("ym"=c("noncontinuous","continuous"), # scale
                           "t"=c(10,50)) # number of mediators

# specific setting
ss <- rep(1:nrow(simsettings),each=50)[seed]
(ymscale <- simsettings[ss,"ym"])
(t <- simsettings[ss,"t"])
# read in population dataset with counterfactual mediators and potential outcomes
load(file = paste0("sim2-potential_outcomes-",ymscale,"-M",t,".Rdata"))

fitY.form <- as.formula(paste0("Y~A+L+",paste("M",1:t,sep="",collapse="+")))

source("funs-noMmodels-Ymodel.R")

nsims <- 40
n_mc <- 5e2
ptm=proc.time()[3]
simres <- lapply(1:nsims, function(sim) {
  
  res_list <- list()
  
  # generate sample and observed values
  n <- 500 # sample size
  ii <- sort(sample.int(max(dat$id),n,replace=FALSE)) # sample individuals
  Ai <- rbinom(n, 1, 0.5) # assign treatment
  Data <- list()
  for (i in 1:n) {
    # revealed outcomes based on assigned treatment
    Data[[i]] <- dat[id==ii[i] & a0==Ai[i] & a1==Ai[i]]
  }
  Data <- rbindlist(Data)
  Data[, c("a.i","a0") := NULL]
  Data[, "id" := 1:n]
  setnames(Data,c("a1","Y.a"),c("A","Y"))
  setkey(Data)
  Data <- data.frame(Data)
  
  # (non-)linear Y model
  res <- OneMCestimator_noMmodels(Data, mc_draws = n_mc,
                                  y_cont = (ymscale=="continuous"))
  
  names(res) <- paste0(names(res),".true")
  res_list <- c(res_list,res)
  rm(res)
  
  # POC estimator
  res <- OnePOCestimator(Data)
  
  names(res) <- paste0(names(res),".poc")
  res_list <- c(res_list,list(res))
  rm(res)
  
  return(c("ym"=ymscale,"t"=t,unlist(res_list)))
})
proc.time()[3]-ptm
# 1 min per sim

warnings() # print warnings

save(simres,file=paste0("interventional-randomA-sim2-hd-",seed,".Rdata"))
q()