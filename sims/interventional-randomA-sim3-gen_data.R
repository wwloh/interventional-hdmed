# generate data for entire population
k <- 5  # number of mediators with non-zero indirect effects
L <- rnorm(N) # observed mediator-outcome confounder 

# split possible mediators into three groups
m_i.ie <- 1:k # indirect effects
m_i.noY <- (k+1):ceiling((t+k)/2) # no effect on Y
m_i.noA <- (ceiling((t+k)/2)+1):t # not affected by A
identical(c(m_i.ie,m_i.noY,m_i.noA),1:t)

fitY.form <- as.formula(paste0("Y~A+L+",paste("M",1:t,sep="",collapse="+")))

fitY.OK <- FALSE
while (!fitY.OK) {
  # randomly assigned treatment
  A <- rbinom(n=N,1,0.5)
  M <- matrix(NA, nrow=N, ncol=t)
  for (s in 1:t) {
    # parameters for mediator models
    alpha_s <- (s==1)*1.0
    a0 <- rnorm(1,sd=0.1)
    al <- rnorm(1,sd=0.1)
    if (s>1 && (s %in% m_i.ie)) {
      bs_lag1 <- 1
    } else if (!(s %in% m_i.ie) & (s != min(m_i.noA))) {
      bs_lag1 <- rnorm(1,sd=0.1)
    } else {
      bs_lag1 <- 0
    }
    if ((s > (k+1)) && (s %in% m_i.noY)) {
      cs_k <- rnorm(1,sd=0.1)
    } else {
      cs_k <- 0
    }
    # observed mediators
    M[,s] <- a0 + alpha_s*A + al*L + rnorm(N,sd=0.1)
    if (s>1) {
      M[,s] <- M[,s] + bs_lag1*M[,s-1]
    }
    if (s>k) {
      M[,s] <- M[,s] + cs_k*M[,k]
    }
  }
  colnames(M) <- paste0("M",1:t)
  # parameters for outcome model
  beta <- rep(0,t)
  beta[m_i.ie] <- 0.2
  beta[m_i.noA] <- rnorm(length(m_i.noA),sd=0.1)
  b0 <- 0; ba <- 0; bl <- rnorm(1,sd=0.1)
  # observed outcome
  Y.ast <- b0 + ba*A + bl*L + (M %*% beta)[,1]
  Y <- rbinom(n=N,size=1,prob=exp(Y.ast)/(1+exp(Y.ast)))
  
  popdat <- data.frame(cbind("id"=1:N,L,A,M,Y))
  
  # check if outcome model can be fitted
  fitY.check <- glm(fitY.form, family = binomial("logit"), data = popdat)
  fitY.OK <- fitY.check$converged & !fitY.check$boundary & 
    min(fitY.check$fitted.values) > .Machine$double.eps*1e1 &
    max(fitY.check$fitted.values) < (1-.Machine$double.eps*1e1)
}

rm(fitY.form)