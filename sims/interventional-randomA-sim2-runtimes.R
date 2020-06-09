rm(list=ls())
libraries_check <- c("data.table")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)
sessionInfo()

# simulation settings
simsettings <- expand.grid("t"=c(9,50,100), # number of mediators
                           "n"=c(1e3,1e4))

# initialize for parallel cluster jobs
args <- 1
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
  simsettings <- simsettings[
    unlist(sapply(1:nrow(simsettings), function(i) 
      rep(i,ifelse(simsettings$n[i]==1e4,1,20))))
    ,]
  nrow(simsettings)
}
(seed <- as.integer(args[1]))
rm(args)

(t <- as.integer(simsettings[seed,"t"])) # number of possible mediators
(n <- as.integer(simsettings[seed,"n"])) # sample size
N <- 1e4 # population size

set.seed(9000) # fix seed for population data
source("interventional-randomA-sim3-gen_data.R")

set.seed(seed=NULL) # re-initializes
Data <- popdat[sample(nrow(popdat),n,replace=FALSE),]
Data$id <- 1:nrow(Data)

# oracle knowledge of the true mediators ##################################
fitY.form <- as.formula(paste0("Y~A+L+",paste("M",1:k,sep="",collapse="+")))
source("funs-noMmodels-Ymodel.R")

ptm=proc.time()[3]
Data_M1k <- Data[,c("id","A","L",paste0("M",1:k),"Y")]
popeffects <- OneMCestimator_noMmodels(data=Data_M1k,
                                       mc_draws=100,
                                       y_cont=FALSE, 
                                       pt_est=TRUE)$pt
(ptm <- proc.time()[3]-ptm)

res <- c("oracle"=popeffects[paste0("t",1:k)],"oracle"=ptm) # save runtime

# all candidate mediators #####################################################
rm(OneMCestimator_noMmodels,fitY.form,Data_M1k,popeffects)
fitY.form <- as.formula(paste0("Y~A+L+",paste("M",1:t,sep="",collapse="+")))
source("funs-noMmodels-Ymodel.R")

ptm=proc.time()[3]
popeffects <- OneMCestimator_noMmodels(data=Data,
                                       mc_draws=100,
                                       y_cont=FALSE, 
                                       pt_est=TRUE)$pt
(ptm <- proc.time()[3]-ptm)

res <- c(res,"all"=popeffects[paste0("t",1:t)],"all"=ptm) # save runtime

warnings() # print warnings

ss <- seed
filename <- paste(c(rbind(names(simsettings),
                          unlist(simsettings[ss,]))),collapse="_")
save(res,file=paste0("interventional-randomA-sim3-hd-",
                     filename,"-",seed,".Rdata"))
q()

# results =====================================================================
setwd("sim3-runtimes/")
myfiles <- list.files()
myfiles <- myfiles[grep(pattern=".Rdata",myfiles)]
sim_results_list <- NULL
for (ll in myfiles) {
  load(ll)
  seed <- as.integer(rev(strsplit(strsplit(ll,"-")[[1]],".Rdata"))[[1]])
  tM <- max(as.integer(unlist(lapply(strsplit(names(res),split="t"),"[",2))),
            na.rm=TRUE)
  sim_results_list[[seed]] <- unlist(list(
    "t"=tM,
    res[grep("elapsed",names(res),value=TRUE)]))
  rm(res)
}
res <- data.table(do.call(rbind,sim_results_list))
setnames(res,"all.elapsed","elapsed")
setkey(res)
res

res[,.N,by=t]
res[,elapsed.perM := elapsed/t]
setkey(res)

ci.emp <- res[,as.list(range(elapsed)/60),by=t]
plot(res[,mean(elapsed/60),by=t],
     xlab="Number of candidate mediators",ylab="Time (minutes)",
     main="Average time (n = 1000)",
     pch=20,ylim=range(ci.emp[,-1]))
for (tt in 1:nrow(ci.emp)) {
  lines(ci.emp[tt,rep(t,2)],ci.emp[tt,-1,with=FALSE],lwd=2)
}

library("xtable")
xtable(merge(res[,mean(elapsed/60),by=t],ci.emp))

ci.emp <- res[,as.list(range(elapsed.perM)/60),by=t]
plot(res[,mean(elapsed.perM/60),by=t],
     xlab="Number of candidate mediators",ylab="Time (minutes)",
     main="Average time per mediator (n = 1000)",
     pch=20,ylim=range(ci.emp[,-1]))
for (tt in 1:nrow(ci.emp)) {
  lines(ci.emp[tt,rep(t,2)],ci.emp[tt,-1,with=FALSE],lwd=2)
}

library("xtable")
xtable(merge(res[,mean(elapsed.perM/60),by=t],ci.emp))
