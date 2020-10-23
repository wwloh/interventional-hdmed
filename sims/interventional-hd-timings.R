rm(list=ls())
libraries_check <- c("data.table","glmnet")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

source("funs-hd-Yonly.R")

alltimes <- NULL
N <- 1e4 # population size
n_mc <- 100 # number of MC draws
y_binary <- FALSE
bsig <- 0
fitA.form <- as.formula("A~1") # randomly assigned treatment
for (t in c(50,100)) {
  source("interventional-randomA-sim3-gen_data.R")
  fitY.form <- as.formula(paste0("Y~A+L+",paste("M",1:t,sep="",collapse="+")))
  if (y_binary) {
    fitY.pop <- glm(formula=fitY.form,family=binomial("logit"),data=popdat)  
  } else {
    fitY.pop <- lm(formula=fitY.form,data=popdat)  
  }
  mytimes <- NULL
  for (i in 1:20) {
    # sample observations
    Data <- popdat[sample(nrow(popdat),1000,replace=FALSE),]
    Data$id <- 1:nrow(Data)
    row.names(Data) <- NULL
    
    mytimes_i <- NULL
    ptm=proc.time()[3]
    cat (i, " ############################################################# \n")
    for (s in 1:t) {
      OneMCestimator_noMmodels_oneIEonly(
        data=Data, 
        mc_draws=n_mc, 
        fitY.coef=fitY.pop$coefficients,
        Ms.i=paste0("M",s))
      
      mytimes_i[[s]] <- proc.time()[3]-ptm
      cat(s, "; time taken (mins) =", round((proc.time()[3]-ptm)/60,1), 
          "; approx. left (mins) =", round((proc.time()[3]-ptm)*(t/s-1)/60),
          "\n")
    }
    mytimes[[i]] <- mytimes_i
  }
  alltimes <- c(alltimes,mytimes)
}
names(alltimes) <- paste0("p",c(50,100))
save(alltimes,file="interventional-hd-timings.Rdata")
q()

load(file="interventional-hd-timings.Rdata")
mean(unlist(alltimes[1:20])/50)/60
# [1] 0.06272628
mean(unlist(alltimes[21:40])/100)/60
# [1] 0.1167481