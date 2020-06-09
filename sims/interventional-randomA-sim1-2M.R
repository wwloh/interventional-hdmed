rm(list=ls())
libraries_check <- c("data.table")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)
sessionInfo()

# simulation settings =========================================================
simsettings <- expand.grid(
  "Um"=c(FALSE,TRUE), # unobserved confounder of mediators
  "corrM"=c(FALSE,TRUE), # correlated mediators
  "contY"=c(TRUE,FALSE), # continuous or binary outcome
  ## settings for observed data analysis
  "n"=c(50,500), # sample size
  "mc"=100 # number of MC draws
)
(simsettings <- simsettings[simsettings$Um >= simsettings$corrM,])

# initialize for parallel cluster jobs
args <- 4
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
  simsettings <- simsettings[rep(1:nrow(simsettings),each=10),]
  nrow(simsettings)
}
(seed <- as.integer(args[1]))
rm(args)

# specific setting
ss <- seed
(nu <- as.numeric(simsettings[ss,"Um"]))
(corrM <- as.numeric(simsettings[ss,"corrM"]))
(contY <- simsettings[ss,"contY"])
(n <- as.integer(simsettings[ss,"n"])) 
(n_mc <- as.integer(simsettings[ss,"mc"]))

# data-generating process =====================================================
t <- 2 # number of mediators

muL <- 1 # mean of observed mediator-outcome confounder
muU <- 0 # mean of unobserved mediator-mediator confounder 

# parameters for mediator models
alpha <- list()
alpha[[1]] <- c("a0"=0,"aA"=0.6,"al"=0.2,"au"=-0.8*nu)
alpha[[2]] <- c("a0"=0,"aA"=0,"a1"=1,"al"=0.2,"au"=-0.8*nu)
alpha

# parameters for outcome model
beta <- c("b0"=0,"bA"=0,"b1"=0.2,"b2"=-0.2,"b12"=0.4*corrM,"bl"=0.2)
beta

if (contY==TRUE) {
  alpha <- lapply(alpha, "*", 2)
}

# covariance of mediator errors under treatment
cov_M1M2_A1 <- 0.8*corrM

# closed-form solutions allowing for M1-M2 interactions
E_M10 <- alpha[[1]]["a0"]+alpha[[1]]["al"]*muL+alpha[[1]]["au"]*muU
## effect of treatment on M2
a_2A <- alpha[[2]]["a1"]*alpha[[1]]["aA"]
E_M20 <- alpha[[2]]["a0"]+alpha[[2]]["al"]*muL+alpha[[2]]["au"]*muU+
  alpha[[2]]["a1"]*E_M10

true_effects <- NULL
## direct effect
true_effects["de"] <- beta["bA"] 
## indirect effect via M1
true_effects["ie1"] <- (beta["b1"]+beta["b12"]*E_M20)*alpha[[1]]["aA"]
## indirect effect via M2
true_effects["ie2"] <- (beta["b2"]+beta["b12"]*E_M10)*a_2A
## indirect effect via mutual dependence
true_effects["ie12"] <- beta["b12"]*cov_M1M2_A1
## difference between joint-sum and mutual
true_effects["ie_diff"] <- beta["b12"]*a_2A*alpha[[1]]["aA"]
true_effects

fitY.form <- as.formula(paste0("Y~A+L+",paste("M",1:t,sep="",collapse="+"),
                               "+",paste("M",1:t,sep="",collapse=":")))

source("funs-sim1-gen_data.R")
source("funs-noMmodels-Ymodel.R")
source("funs-NEmodel.R")

# generate fixed population and true effects
set.seed(9000)
popdat <- OneData(N=1e5,contY=contY)
pop_effects <- list()
pop_effects[["ie"]] <- OneMCestimator_noMmodels(data = popdat, 
                                                mc_draws = n_mc, 
                                                y_cont = contY, 
                                                pt_est = TRUE)$pt
## true parameter values from data-generating mechanism
oracle_M1_par <- list()
oracle_M1_par[["alpha"]] <- alpha[[1]]
oracle_M1_par[["sigma"]] <- 1
pop_effects[["ne"]] <- OneSteenEstimator(data=popdat,y_cont=contY,oracle_M1_par)
## remove oracle estimators
pop_effects[["ne"]][["ne.o"]] <- NULL
pop_effects <- unlist(pop_effects)
names(pop_effects) <- paste0(names(pop_effects),".pop")
pop_effects
true_effects; sum(true_effects[-1])

# estimates for sampled data ==================================================
set.seed(seed)
nsims <- 100
ptm=proc.time()[3]
simres <- lapply(1:nsims, function(sim) {
  
  res_list <- list()
  
  # generate sample and observed values
  fitY.OK <- FALSE
  while (!fitY.OK) {
    Data <- popdat[sort(sample(nrow(popdat),n,replace=FALSE)),]
    if (contY==FALSE) {
      # check if outcome model can be fitted
      fitY.check <- glm(fitY.form, family = binomial("logit"), data = Data)
      fitY.OK <- fitY.check$converged & !fitY.check$boundary & 
        min(fitY.check$fitted.values) > .Machine$double.eps*1e1 &
        max(fitY.check$fitted.values) < (1-.Machine$double.eps*1e1)
      rm(fitY.check)
    } else {
      fitY.OK <- TRUE
    }
  }
  rownames(Data) <- NULL
  Data[,"id"] <- 1:n
  
  # (non-)linear Y model
  res <- OneMCestimator_noMmodels(data = Data, 
                                  mc_draws = n_mc, 
                                  y_cont = contY, 
                                  pt_est = TRUE)$pt
  
  names(res) <- paste0(names(res),".intv")
  res_list <- c(res_list,res)
  rm(res)
  
  # NE estimator
  res <- OneSteenEstimator(data=Data,y_cont=contY,oracle_M1_par)
  
  names(res) <- paste0(names(res),".natu")
  res_list <- c(res_list,list(res))
  rm(res)
  
  # misspecify causal ordering by swapping indices of M1 and M2
  correct_M1 <- Data[,"M1"]
  correct_M2 <- Data[,"M2"]
  Data[,"M1"] <- correct_M2
  Data[,"M2"] <- correct_M1
  res <- OneMCestimator_noMmodels(data = Data, 
                                  mc_draws = n_mc, 
                                  y_cont = contY, 
                                  pt_est = TRUE)$pt
  
  names(res) <- paste0(names(res),".intv.m2m1")
  res_list <- c(res_list,res)
  rm(res)
  
  # NE estimator
  res <- OneSteenEstimator(data=Data,y_cont=contY,oracle_M1_par)
  ## remove oracle estimators
  res$ne.o <- NULL
  names(res) <- paste0(names(res),".natu.m2m1")
  res_list <- c(res_list,list(res))
  rm(res)
  
  return( c(simsettings[ss,],pop_effects,unlist(res_list)) )
})
proc.time()[3]-ptm
# 9 sec per sim for n=500 and 100 MC draws

warnings() # print warnings

filename <- paste(c(rbind(names(simsettings),
                          unlist(simsettings[ss,]))),collapse="_")
save(simres,file=paste0("interventional-randomA-sim1-2M-",
                        filename,"-",seed,".Rdata"))
q()
