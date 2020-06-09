rm(list=ls())
libraries_check <- c("data.table","glmnet")#,"coxed")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)
sessionInfo()

# simulation settings
simsettings <- expand.grid("t"=c(9,50,100), # number of mediators
                           "n"=c(1,2,10,40)*50)
simsettings <- simsettings[simsettings$t==100,]

# initialize for parallel cluster jobs
args <- 12
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
  simsettings <- simsettings[rep(1:nrow(simsettings),each=1000),]
  nrow(simsettings)
}
(seed <- as.integer(args[1]))
rm(args)

(t <- as.integer(simsettings[seed,"t"])) # number of possible mediators
(n <- as.integer(simsettings[seed,"n"])) # sample size
N <- 1e4 # population size

set.seed(9000) # fix seed for population data
source("interventional-randomA-sim3-gen_data.R")

fitA.form <- as.formula("A~1") # randomly assigned treatment
source("dataex-Rscripts/funs-double_selection.R")

# estimates for a given sample ================================================
OneEst <- function(Data, n_mc) {
  
  # save model matrix **with all mediators** for predicting potential outcomes
  Data_M <- Data
  Data_M$id <- NULL
  Data_M$Y <- NULL
  Data_M <- cbind("id"=1:nrow(Data_M),"(Intercept)"=1,Data_M)
  
  # double selection ############################################################
  res <- DoubleSelection(L=Data[,"L",drop=FALSE],
                         A=Data[,"A"],
                         M=Data[,grep(pattern="M",names(Data),value=TRUE)],
                         Y=Data[,"Y"],
                         fitY.family="binomial",
                         fitM.family="gaussian")
  
  # estimate indirect effects via each mediator using selected predictors #######
  res.ds <- NULL
  ptm=proc.time()[3]
  for (s in 1:length(res)) {
    # unique predictors selected in either mediator or outcome model
    dat.Ms.allnames <- unique(unlist(lapply(res[[s]], function(x)
      names(x[abs(x)>.Machine$double.eps]))))
    dat.Ms.allnames <- c(dat.Ms.allnames,"Y")
    dat.Ms <- Data[,grep("(Intercept)", dat.Ms.allnames,invert=TRUE,value=TRUE),
                   drop=FALSE]
    # check for non-convergence or boundary solutions
    fitY <- glm(Y~., family = binomial("logit"), data = dat.Ms)
    fitY.OK <- fitY$converged & !fitY$boundary
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
                  mc_draws=n_mc,
                  fitY.coef=fitY.coef.Ms,
                  Ms.i=paste0("M",s)))
    
    cat("M",s,":",res.ds[s],"\n")
    cat("time taken (mins) =", round((proc.time()[3]-ptm)/60), 
        "; approx. left (mins) =", round((proc.time()[3]-ptm)*(length(res)/s-1)/60),
        "\n")
  }
  return(res.ds)
}

set.seed(seed=NULL) # re-initializes
nsims <- 1
stm <- proc.time()[3]
simres <- lapply(1:nsims, function(sim) {
  # sample observations
  Data <- popdat[sample(nrow(popdat),n,replace=FALSE),]
  Data$id <- 1:nrow(Data)
  row.names(Data) <- NULL
  
  # point estimates
  res.ds <- OneEst(Data,n_mc=100)
  
  if (t<1000) {
    n.boots <- 50
    # ordinary bootstrap
    res.bca <- do.call(rbind,lapply(1:n.boots, function(i) {
      OneEst(Data[sample(n,replace=TRUE),],n_mc=1)
    }))
    # res.bca <- apply(res.bca, 2, coxed::bca)
  } else {
    res.bca <- matrix(NA,nrow=2,ncol=length(res.ds))
  }
  
  return( list("t"=t,"n"=n,"res"=res.ds,
               "res.bca"=res.bca) )
})
proc.time()[3]-stm
# 25 min per sim for 200 mediators with n=200
# simres <- do.call(rbind,lapply(simres, unlist))

filename <- paste(c(rbind(names(simsettings),
                          unlist(simsettings[seed,]))),collapse="_")
save(simres,file=paste0("interventional-randomA-sim3-hd-ds-",
                        filename,"-",seed,".Rdata"))
q()

# results =====================================================================
k <- 5  # number of mediators with non-zero indirect effects
m_i.ie <- 1:k # indirect effects
## population values
subfolder <- "sim3-popvals/"
myfiles <- list.files(subfolder)
myfiles <- grep(pattern=".Rdata",myfiles,value=TRUE)
res_list <- list()
for (ll in myfiles) {
  # results
  load(paste0(subfolder,ll))
  tM <- max(as.integer(unlist(lapply(strsplit(names(res),split="t"),"[",2))),
            na.rm=TRUE)
  oracle.ie <- res[paste0("oracle.",paste0("t",m_i.ie))]
  all.ie <- res[paste0("all.",paste0("t",m_i.ie))]
  names(oracle.ie) <- names(all.ie) <- paste0("t",m_i.ie) 
  res_list <- c(res_list, list(
    # as.list(c("M"=tM,"est"=0,oracle.ie)),
    as.list(c("M"=tM,"est"=1,all.ie))))
  rm(res)
  cat(ll,"\n")
}
popeff <- rbindlist(res_list,fill=TRUE)
setkey(popeff)
popeff

library("coxed")
## sample estimates
subfolder <- "sim3-ds/"
myfiles <- list.files(subfolder)
myfiles <- grep(pattern=".Rdata",myfiles,value=TRUE)
res_list <- list()
res_list.ci <- list()
for (ll in myfiles) {
  # results
  load(paste0(subfolder,ll))
  simres.ci <- lapply(simres, "[[", "res.bca")
  simres <- do.call(rbind,lapply(simres, function(x) unlist(x[c("t","n","res")])))
  simsetting <- unique(simres[,c("t","n")])
  simsetting <- paste(c(rbind(c("t","n"),simsetting)),collapse="_")
  res_list[[simsetting]] <- c(res_list[[simsetting]],list(data.frame(simres)))
  if (!all(unlist(lapply(simres.ci, function(x) apply(x,2,is.na))))) {
    simres.ci <- lapply(simres.ci, function(x) 
      cbind("t"=unique(simres[,"t"]),"n"=unique(simres[,"n"]),
            # bca-transformed CIs
            "M"=1:ncol(x),t(apply(x,2,coxed::bca))))
    row.names(simres.ci) <- NULL
    res_list.ci[[simsetting]] <- c(res_list.ci[[simsetting]],
                                   list(data.frame(do.call(rbind,simres.ci))))
  }
  rm(simres,simres.ci)
  cat(ll,"\n")
}

lapply(res_list,function(x) sum(unlist(lapply(x,nrow))))

res.summ <- lapply(res_list, function(tt) {
  tt.df <- do.call(rbind,tt)
  # sample estimates (for all mediators) ####
  esteffects <- tt.df[,grep("res.t",names(tt.df),value=TRUE)]
  # ranks for all possible mediators
  tt.ranks <- t(apply(esteffects, 1, function(tt.x) 
    rank(-abs(tt.x),na.last="keep")))
  tt.ranks.k <- t(apply(tt.ranks, 1, "%in%", m_i.ie)) # ranked in top 5
  colnames(tt.ranks.k) <- colnames(tt.ranks)
  
  # population values
  popeffects <- unlist(popeff[M==ncol(esteffects) & est==1,
                              paste0("t",m_i.ie), with=FALSE])
  names(popeffects) <- names(esteffects)[m_i.ie]
  
  # results for true mediators
  tt.res <- rbind(
    "nsims"=rep(nrow(esteffects),length(m_i.ie)),
    "true.hd"=popeffects[m_i.ie],
    "est.hd"=colMeans(esteffects[,m_i.ie]),
    "se.hd"=apply(esteffects[,m_i.ie],2,sd),
    "rank.av"=colMeans(tt.ranks[,m_i.ie]),
    "rank.sd"=apply(tt.ranks[,m_i.ie],2,sd),
    "rank.k.av"=colMeans(tt.ranks.k[,m_i.ie])
    )
  
  # other mediators for comparison
  m_compare <- names(c(
    which.max(abs(colMeans(esteffects[,-(m_i.ie)]))) # largest abs. magnitude
  ))
  
  tt.res_compare <- rbind(
    "nsims"=rep(nrow(esteffects),length(m_compare)),
    "true.hd"=rep(0,length(m_compare)),
    "est.hd"=colMeans(esteffects[,m_compare, drop=FALSE]),
    "se.hd"=apply(esteffects[,m_compare, drop=FALSE],2,sd),
    "rank.av"=colMeans(tt.ranks[,m_compare, drop=FALSE]),
    "rank.sd"=apply(tt.ranks[,m_compare, drop=FALSE],2,sd),
    "rank.k.av"=colMeans(tt.ranks.k[,m_compare, drop=FALSE])
  )
  
  return( cbind(tt.res,tt.res_compare) )
})
res.summ
library("xtable")
lapply(res.summ[grep("t_9",names(res.summ))], xtable)
lapply(res.summ[grep("t_50",names(res.summ))], xtable)
lapply(res.summ[grep("t_100",names(res.summ))], xtable)

# plots 
for (tt in c(50,100)) {
  res <- rbindlist(lapply(res_list, rbindlist),fill=TRUE)
  res <- res[t==tt]
  setkey(res)
  
  pdf(paste0("plot-sim3-est-",tt,".pdf"),width=8,height=6)
  par(mfrow=c(2,3))
  for (s in 1:k) {
    est.density.s <- NULL
    for (ni in unique(res$n)) {
      est.density.s <- c(
        est.density.s,
        list(density(unlist(res[n==ni,paste0("res.t",s),with=FALSE]))))
    }
    names(est.density.s) <- unique(res$n)
    # my_xlim <- range(est.density.s[[which.max(unique(res$n))]]$x)
    my_xlim <- range(est.density.s[[which(names(est.density.s)==tt)]]$x)
    
    plot(est.density.s[[1]],
         # main=paste0("M",s ," (ds)"),
         main=paste0("M",s),
         xlab="Indirect effect estimates",
         xlim=my_xlim,
         ylim=range(unlist(lapply(est.density.s, function(x) x$y))),
         col=0)
    for (i in 1:length(est.density.s)) {
      lines(est.density.s[[i]], col=i)
    }
    legend("topright",legend=paste0("n=",unique(res$n)),
           col=1:length(est.density.s),lty=1,bty="n",cex=.8)
    abline(v=popeff[M==unique(res$t) & est==1,paste0("t",s), with=FALSE],lty=2)
  }
  dev.off()
  rm(est.density.s)
}

res_list.ci <- do.call(rbind,lapply(res_list.ci, function(x) do.call(rbind,x)))
row.names(res_list.ci) <- NULL
res_list.ci <- data.table(res_list.ci)
setkey(res_list.ci)
res_list.ci[,.N,by=list(t,n,M)][,unique(N),by=list(t,n)]

pM <- 100
res_list.ci_t <- res_list.ci[t==pM]
maxM <- as.integer(strsplit(colnames(
  res.summ[grep(paste0("t_",pM),names(res.summ))][[1]])[k+1],"t")[[1]][2])
for (s in 1:pM) {
  if (s <= k) {
    ie.true <- as.numeric(popeff[est==1 & M==pM,paste0("t",s),with=FALSE])
  } else {
    ie.true <- 0
  }
  res_list.ci_t[M==s,"true" := ie.true]
}
res_list.ci_t[, c("cover","halfwidth") := list(V4 <= true & V5 >= true, (V5-V4)/2)]
for (nn in unique(res_list.ci_t$n)) {
  print(xtable(t(res_list.ci_t[n==nn & (M<=k | M==maxM),lapply(.SD, mean),
                             by=list(n,M), .SDcols=c("cover","halfwidth")])))
}

