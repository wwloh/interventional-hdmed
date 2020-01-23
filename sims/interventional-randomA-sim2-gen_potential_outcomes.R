rm(list=ls())
libraries_check <- c("data.table")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)
sessionInfo()
set.seed(9000)

# simulation settings
simsettings <- expand.grid("ym"=c("noncontinuous","continuous"), # scale
                           "t"=c(10,50)) # number of mediators

for (ss in 1:nrow(simsettings)) {
  ymscale <- simsettings[ss,"ym"]
  t <- simsettings[ss,"t"]
  
  N <- 1e5 # population size
  a0 <- 0; al <- 0.25; nu <- 0.25 # parameters for mediator models
  b0 <- -1; ba <- 0.5; bl <- 1 # parameters for outcome model
  
  alpha <- rep(0,t)
  alpha[1] <- 1
  alpha[(t/2+1):t] <- 1.5
  
  beta <- rep(0,t)
  beta[1] <- 1
  beta[2:(t/2)] <- 1.5
  
  Bmat <- matrix(0, nrow=t, ncol=t)
  for (s in 2:t) {
    for (j in 1:(s-1)) {
      if (beta[s] != 0) {
        Bmat[s,j] <- (alpha[j]==0)*exp(j-s)
      } else if (beta[s] == 0) {
        Bmat[s,j] <- 0.5*(j==1)
      }
    }
  }
  IBmat.inv <- solve(diag(1,nrow=t,ncol=t)-Bmat)
  
  L <- rnorm(N)
  U <- rnorm(N)

  # generate counterfactual mediator values under treatment and control
  Mast <- list()
  for (aa in 0:1) {
    A <- rep(aa, N)
    Mast_a <- IBmat.inv %*% 
      ((alpha %o% A) + matrix(rep(a0+al*L+nu*U,each=t), 
                              nrow=t,ncol=N,byrow=FALSE))
    for (ss in 1:t) {
      if (ymscale=="noncontinuous") {
        Mast_a[ss,] <- rbinom(n=ncol(Mast_a), 1, 
                              exp(Mast_a[ss,])/(1+exp(Mast_a[ss,])))
      } else {
        Mast_a[ss,] <- Mast_a[ss,] + rnorm(n=ncol(Mast_a))
      }
    }
    M_a <- data.frame(t(Mast_a))
    colnames(M_a) <- paste0("M",1:t)
    Mast[[aa+1]] <- M_a
  }
  names(Mast) <- paste0("A",0:1)
  # lapply(Mast, colMeans)
  
  # duplicated data for each individual for direct and (joint) indirect effects
  dat <- data.table(id = 1:N, L)
  setkey(dat)
  ## hypothetical treatment levels
  a0.vals <- c(0,0,1)
  a1.vals <- c(0,1,1)
  ## set counterfactual mediator values
  dat.list <- list()
  for (ai in 1:3) {
    dat.ai <- cbind(dat,"a.i"=ai,"a0"=a0.vals[ai],"a1"=a1.vals[ai])
    if (a1.vals[ai]==0) {
      dat.ai <- cbind(dat.ai, Mast$A0)
    } else {
      dat.ai <- cbind(dat.ai, Mast$A1) 
    }
    dat.list[[ai]] <- dat.ai
  }
  dat <- rbindlist(dat.list)
  setkey(dat)
  
  # set potential outcomes
  Ya.ast <- b0 + ba*dat$a0 + bl*dat$L + 
    tcrossprod(as.matrix(dat[, grepl("M",key(dat)), with=FALSE]),t(beta))[,1]
  if (ymscale=="noncontinuous") {
    Y.a <- rbinom(n=length(Ya.ast),1,prob=exp(Ya.ast)/(1+exp(Ya.ast)))
  } else {
    Y.a <- Ya.ast + rnorm(n=length(Ya.ast), sd=sqrt(t))
  }
  dat <- cbind(dat,Y.a)
  setkey(dat)
  
  # save population values
  save(dat,file = paste0("sim2-potential_outcomes-",ymscale,"-M",t,".Rdata"))
  
  # average over all individuals
  dat <- dat[, lapply(.SD,mean), by=list(a.i,a0,a1), .SDcols="Y.a"]
  # true values of direct and joint indirect effects
  if (ymscale=="noncontinuous") {
    logit <- function(x) log(x/(1-x))
    gamma <- c(logit(dat[3, Y.a])-logit(dat[2, Y.a]),
               logit(dat[2, Y.a])-logit(dat[1, Y.a]))
  } else {
    gamma <- c(dat[3, Y.a]-dat[2, Y.a],dat[2, Y.a]-dat[1, Y.a])
  }
  names(gamma) <- paste0("g",0:1)
  cat(as.character(ymscale),t,"\n")
  print(gamma)
}
# noncontinuous 10 
# g0        g1 
# 0.3351156 0.1438164 
# continuous 10 
# g0        g1 
# 0.4787982 1.0002125 
# noncontinuous 50 
# g0        g1 
# 0.2928249 0.1630121 
# continuous 50 
# g0        g1 
# 0.4746873 1.0826835 