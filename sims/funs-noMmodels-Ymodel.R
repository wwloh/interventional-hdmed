# helper function to create duplicated data for one individual
Dupdata <- function(t) {
  out <- diag(1,nrow=t,ncol=t) # for indirect effects via each mediator
  out <- rbind(0,out,1)  # first row all 0s, last row all 1s
  out <- cbind(0,out) # first column all zeroes for a0
  out <- rbind(0,out[nrow(out),],1,out) # for direct, joint, and mutual effects
  colnames(out) <- paste0("a",0:t)
  out <- cbind("Marg"=c(rep(0,3),rep(1,nrow(out)-3)),"a.i"=1:nrow(out),out)
  return(out)
}

OneMCestimator_noMmodels <- function(data,mc_draws=1e2,y_cont,pt_est=TRUE) {
  res <- list()

  # fit outcome model =========================================================
  if (y_cont==TRUE) {
    fitY <- glm(fitY.form, family = gaussian("identity"), data = data)  
  } else {
    fitY <- glm(fitY.form, family = binomial("logit"), data = data)  
  }
  
  # helper function to predict Y for different outcome models
  PredictY <- function(onedat) {
    ## first column is always treatment
    colnames(onedat)[1] <- "A"
    Ya <- predict.glm(fitY, type="response", newdata=onedat)
    return(Ya)
  }
  
  # mediator column names
  Mnames <- colnames(data)[grep("M",colnames(data))]
  t <- length(Mnames) # number of distinct mediators
  ## observed mediator values in each treatment group
  Mobs <- list(data[data$A==0, Mnames],data[data$A==1, Mnames])
  n_A0 <- sum(data$A==0)
  n_A1 <- sum(data$A==1)
  
  # duplicated data for each individual =======================================
  dat <- data.table(data[,c("id","L")]) #### keep only id, L to reduce memory
  setkey(dat)
  alevels <- dat[,as.data.table(Dupdata(t)),by=id]
  setkey(alevels)
  dat <- merge(dat,alevels,all.x=TRUE)
  setkey(dat)
  rm(alevels)
  
  SampleMs <- function(mydt,av_mc) {
    # mydt is a data.table for each id containing the relevant columns
    # e.g., mydt <- dat[id==1]
    mydt_mc <- mydt[rep(1:nrow(mydt),each=mc_draws)]
    mydt_mc <- cbind(mydt_mc,"mc"=rep(1:mc_draws,times=nrow(mydt)))
    setkey(mydt_mc)
    
    Mtilde <- data.frame(matrix(NA,nrow=nrow(mydt_mc),ncol=t))
    colnames(Mtilde) <- Mnames
    
    ## sample *joint* mediator values
    #### a1==0
    a1_0 <- mydt_mc$Marg==0 & mydt_mc$a1==0
    Mtilde[a1_0,] <- Mobs[[1]][sample(n_A0, size=sum(a1_0), replace=TRUE),]
    rm(a1_0)
    #### a1==1
    a1_1 <- mydt_mc$Marg==0 & mydt_mc$a1==1
    Mtilde[a1_1,] <- Mobs[[2]][sample(n_A1, size=sum(a1_1), replace=TRUE),]
    rm(a1_1)
    
    ## sample *marginal* mediator values
    M.as <- lapply(1:t, function(s) {
      a_s <- unlist(mydt_mc[Marg==1,paste0("a",s),with=FALSE])
      Ms.as <- rep(NA,length(a_s))
      Ms.as[a_s==0] <- Mobs[[1]][sample(n_A0, size=sum(a_s==0), replace=TRUE),s]
      Ms.as[a_s==1] <- Mobs[[2]][sample(n_A1, size=sum(a_s==1), replace=TRUE),s]
      return(Ms.as)
    })
    Mtilde[mydt_mc$Marg==1,] <- do.call(cbind,M.as); rm(M.as)
    mydt_mc <- cbind(mydt_mc,Mtilde); rm(Mtilde)
    setkey(mydt_mc)
    
    ## predict potential outcomes using sampled mediator values
    mydt_mc[, "Y.a" := PredictY(onedat=mydt_mc[, c("a0","L",Mnames), with=FALSE])]
    setkey(mydt_mc)
    if (av_mc==TRUE) {
      ## average over MC draws
      mydt_mc_means <- mydt_mc[,lapply(.SD, mean),
                               by=c("Marg","a.i",paste0("a",0:t)),.SDcols="Y.a"]
      setkey(mydt_mc_means)
      return(mydt_mc_means[,Y.a])
    } else {
      return(mydt_mc)  
    }
  }
  
  if (pt_est==TRUE) {
    # average potential outcomes for each individual only
    dat[, "Y.a" := SampleMs(.SD,av_mc=TRUE), by=id]
    setkey(dat)
    ## average over all individuals
    mu_hat <- dat[, lapply(.SD,mean), 
                  by=c("Marg","a.i",paste0("a",0:t)),.SDcols="Y.a"]
    setnames(mu_hat,old="Y.a",new="Y")
    setkey(mu_hat)
  }
  
  # Delta transformation matrix for effects
  eff_mat <- matrix(0,nrow=t+4,ncol=nrow(mu_hat))
  eff_mat[1,2:3] <- c(-1,1) # direct effect
  eff_mat[2,1:2] <- c(-1,1) # joint indirect effect
  for (s in 1:(t+1)) {
    eff_mat[s+2,c(4,s+4)] <- c(-1,1) # indirect effect via Ms
  }
  # indirect effect via mutual
  eff_mat[nrow(eff_mat),] <- eff_mat[2,] - eff_mat[nrow(eff_mat)-1,]
  
  if (fitY$family$family=="binomial") {
    mu_hat[, Y := log(Y/(1-Y))] # logit function
  }
  est <- (eff_mat %*% mu_hat[,Y])[,1]
  names(est) <- c(paste0("g",0:1),paste0("t",1:t),"t_margsum","mu")
  
  if (pt_est==TRUE) {
    var_est <- NULL
  }
  return( list("pt"=est, "vcov"=var_est) )
}
