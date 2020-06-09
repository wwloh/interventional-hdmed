library("mvtnorm")
OneData <- function(N,contY) {
  
  t <- 2 # number of mediators
  
  L <- rnorm(N,mean=muL) # observed mediator-outcome confounder 
  U <- rnorm(N,mean=muU) # unobserved mediator-mediator confounder 
  A <- rbinom(N, 1, 0.5) # randomly assigned treatment
  M <- matrix(NA,nrow=N,ncol=t) # generate mediator values
  ## correlated error terms for mediators that may depend on treatment
  sigmaM_A0 <- diag(rep(1,t))
  sigmaM_A1 <- matrix(cov_M1M2_A1,nrow=t,ncol=t)
  diag(sigmaM_A1) <- 1
  M_eps <- matrix(NA,nrow=N,ncol=t)
  M_eps[A==0,] <- rmvnorm(n=N-sum(A),mean=rep(0,t),sigma=sigmaM_A0)
  M_eps[A==1,] <- rmvnorm(n=sum(A),mean=rep(0,t),sigma=sigmaM_A1)
  for (pp in 1:t) {
    M[,pp] <- alpha[[pp]]["a0"] + alpha[[pp]]["aA"]*A + 
      alpha[[pp]]["al"]*L + alpha[[pp]]["au"]*U
    if (pp>1) {
      # A --> M1 --> M2
      M[,pp] <- M[,pp] + alpha[[pp]]["a1"]*M[,1]
    }
    M[,pp] <- M[,pp] + M_eps[,pp]
  }
  colnames(M) <- paste0("M",1:t)
  
  # outcomes
  Y.ast <- beta["b0"] + beta["bA"]*A + beta["bl"]*L + 
    (M %*% beta[paste0("b",1:t)])[,1] + beta["b12"]*M[,"M1"]*M[,"M2"]
  if (contY==TRUE) {
    Y <- Y.ast + rnorm(N)
  } else {
    Y <- rbinom(n=N,size=1,prob=exp(Y.ast)/(1+exp(Y.ast)))
  }
  
  Data <- data.frame(cbind("id"=1:N,L,A,M,Y))
  return(Data)
}
