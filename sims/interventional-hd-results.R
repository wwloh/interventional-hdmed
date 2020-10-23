rm(list=ls())
libraries_check <- c("data.table","xtable")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

# results =====================================================================
nM <- 100 # number of mediators
# nM <- 50 # number of mediators
subfolder <- paste0("sim3-yonly-",nM,"M/")
myfiles <- list.files(subfolder)
myfiles <- grep(pattern=".Rdata",myfiles,value=TRUE)
res_list <- list()
for (ll in myfiles) {
  # results
  load(paste0(subfolder,ll))
  res_list <- c(res_list,simres)
  rm(simres)
}
rm(myfiles)

# true values
subfolder <- paste0("popeffs-",nM,"M/")
myfiles <- list.files(subfolder)
myfiles <- grep(pattern=".Rdata",myfiles,value=TRUE)
popeffs <- list()
for (ll in myfiles) {
  # results
  load(paste0(subfolder,ll))
  popeff <- unlist(popeff)
  names(popeff) <- paste0("m",1:length(popeff))
  popeffs <- c(popeffs,list(c(unlist(pop.setting),popeff)))
  rm(popeff)
}
ie_true <- data.table(do.call(rbind,popeffs))
setkey(ie_true)

# calculate CIs
simres <- lapply(res_list, function(x) {
  # meths <- unique(unlist(lapply(strsplit(names(x$res),split="[.]"),"[",1)))
  meths <- c("est","true")
  x_ci <- apply(x$boot,2,quantile, probs=c(.025,.975), na.rm=TRUE)
  x_true <- ie_true[x$settings]
  x_summ_meths <- NULL
  for (mm in meths) {
    if (mm != "true") {
      x_ci_mm <- x_ci[,grep(mm,x=colnames(x_ci),value=TRUE)]
      # half-width
      x_ci_mm.hw <- apply(x_ci_mm, 2, diff)
      # coverage
      x_ci_mm.cover <- x_ci_mm <- sapply(1:ncol(x_ci_mm), function(j) {
        (x_ci_mm[1,j] <= x_true[,paste0("m",j),with=FALSE]) &
          (x_ci_mm[2,j] >= x_true[,paste0("m",j),with=FALSE])
      })  
      x_summ_mm <- rbind("pt"=x$res[grep(mm,x=names(x$res),value=TRUE)],
                         "cover"=x_ci_mm.cover,"halfwidth"=x_ci_mm.hw)
    } else {
      x_summ_mm <- rbind("pt"=x$res[grep(mm,x=names(x$res),value=TRUE)])
    }
    colnames(x_summ_mm) <- paste0("m",1:ncol(x_summ_mm))
    x_summ_mm
    x_summ_meths <- rbind(x_summ_meths,
                          data.frame("meth"=mm,"est"=row.names(x_summ_mm),
                                     x_summ_mm))
  }
  return( data.frame(x$settings,x_summ_meths,row.names = NULL) )
})
simres <- rbindlist(simres,fill=TRUE)
setkey(simres)

k <- 5  # number of mediators with non-zero indirect effects
m_i.ie <- 1:k # indirect effects
simsettings <- c("t","n","bsig","ybin","meth","est")
(simres_unique <- unique(simres[,list(t,n,bsig,ybin)]))

# number of sims
simres[, .N, by=simsettings]
# average estimate
simres_est <- simres[est=="pt", lapply(.SD,mean,na.rm=TRUE), by=simsettings]
setkey(simres_est)
# empirical SE
simres_ese <- simres[est=="pt", lapply(.SD, sd, na.rm=TRUE), by=simsettings]
setkey(simres_ese)
# ranks (of absolute magnitudes)
simres_ptest <- simres[est=="pt"]
setkey(simres_ptest)
M_Ranks <- function(onerow) {
  m_t <- as.integer(onerow[,"t"])
  m_names <- paste0("m",1:m_t)
  m_ranks <- rank(-abs(as.numeric(onerow[,..m_names])),na.last="keep")
  names(m_ranks) <- paste0("rank",1:m_t)
  m_topk <- m_ranks %in% m_i.ie # ranked in top 5
  names(m_topk) <- paste0("topk",1:m_t)
  as.list(c(onerow[,..simsettings],m_ranks,m_topk))
}
simres_ranks <- simres_ptest[,M_Ranks(.SD), by=row.names(simres_ptest)]
simres_ranks[,row.names := NULL]
setkey(simres_ranks)
simres_ranks <- simres_ranks[, lapply(.SD,mean,na.rm=TRUE), by=simsettings]
# coverage
simres_cover <- simres[est=="cover", lapply(.SD,mean,na.rm=TRUE), by=simsettings]
setkey(simres_cover)

for (ss in 1:nrow(simres_unique)) {
  for (mymeth in c("est","true")) {
    ss_list <- list()
    ss_list[["true"]] <- ie_true[simres_unique[ss]]
    ss_list[["est"]] <- simres_est[meth==mymeth][simres_unique[ss]]
    ss_list[["ese"]] <- simres_ese[meth==mymeth][simres_unique[ss]]
    ss_list[["cover"]] <- simres_cover[meth==mymeth][simres_unique[ss]]
    ss_list[["ranks"]] <- simres_ranks[meth==mymeth][simres_unique[ss]]
    # highest ranked
    nonm <- as.integer(apply(ss_list$ranks, 1, function(m) {
      m_t <- as.integer(m["t"])
      which.max(m[paste0("rank",(k+1):m_t)])+k
    }))
    # (two-sided) tail probability that biases were non-zero
    biaszero.pv <- matrix(sapply(1:(k+1), function(j) {
      mu <- ss_list$true[,paste0("m",c(1:k,nonm))[j],with=FALSE]
      mu_hat <- ss_list$est[,paste0("m",c(1:k,nonm))[j],with=FALSE]
      mu_se <- ss_list$ese[,paste0("m",c(1:k,nonm))[j],with=FALSE]
      pnorm(abs(unlist(mu_hat)-unlist(mu))/unlist(mu_se),lower.tail=FALSE)*2
    }),nrow=1)
    
    simres_ss <- rbind(
      ss_list$true[,paste0("m",c(1:k,nonm)),with=FALSE],
      ss_list$est[,paste0("m",c(1:k,nonm)),with=FALSE],
      ss_list$ese[,paste0("m",c(1:k,nonm)),with=FALSE],
      ss_list$cover[,paste0("m",c(1:k,nonm)),with=FALSE],
      ss_list$ranks[,paste0("rank",c(1:k,nonm)),with=FALSE],
      ss_list$ranks[,paste0("topk",c(1:k,nonm)),with=FALSE],
      biaszero.pv,
      use.names=FALSE
    )
    print(xtable(
      cbind(simres_unique[ss],
            "summ"=c(names(ss_list),"topk","biaszero.pv"),simres_ss)),
      include.rownames=FALSE)
    
    if (mymeth=="est") {
      # point estimates
      pt_est <- simres[meth==mymeth & est=="pt"][simres_unique[ss]][
        ,paste0("m",c(1:k,nonm)),with=FALSE]
      filename <- paste(c(rbind(names(simres_unique[ss]),
                                unlist(simres_unique[ss]))),collapse="_")
      filename <- gsub(pattern="[.]",replacement="",filename)
      pdf(paste0("plot-sim3-est-",filename,".pdf"),width=8,height=6)
      par(mfrow=c(2,3))
      for (s in 1:ncol(pt_est)) {
        est.density.s <- density(unlist(pt_est[,s,with=FALSE]),na.rm=TRUE)
        my_xlim <- range(est.density.s$x)
        plot(est.density.s,
             main=gsub("m","M",colnames(pt_est)[s]),
             xlab="Indirect effect estimates",
             xlim=my_xlim,
             ylim=range(est.density.s$y),
             col=0)
        lines(est.density.s, col=1)
        abline(v=ss_list$true[,paste0("m",c(1:k,nonm)),with=FALSE][,s,with=FALSE],
               lty=2)
      }
      dev.off()
      rm(est.density.s)
    }
  }
}
