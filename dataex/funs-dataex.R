# helper function to create duplicated data for one individual
Dupdata <- function(t) {
  out <- diag(1,nrow=t+1,ncol=t+1)
  out[1,] <- 0 # first row all 0s
  out <- rbind(out, 1) # last row all 1s except first entry
  out[nrow(out),1] <- 0
  colnames(out) <- paste0("a",0:t)
  out <- cbind("a.i"=1:nrow(out),out)
  return(out)
}

OneMCestimator_noMmodels <- function(data,mc_draws=2) {
  res <- list()

  ## fit outcome model --------------------------------------------------------
  fitY <- glm(fitY.form, family = binomial("logit"), data = data)
  # helper function to predict Y for different outcome models
  PredictY <- function(onedat) {
    ## first column is always treatment
    colnames(onedat)[1] <- "A"
    Ya <- predict.glm(fitY, type="response", newdata=onedat)
    return(Ya)
  }
  
  # sampling probabilities for each individual's mediator values
  fitA <- glm(fitA.form, family = binomial("logit"), data = data)
  pA1hat <- predict(fitA,type="response") # predicted prob. of treatment
  pA0hat <- 1-pA1hat # predicted prob. of control
  wt_A1 <- 1/pA1hat # inverse prob. of treatment weight
  wt_A0 <- 1/pA0hat # inverse prob. of control weight
  ind_A1 <- which(data$A==1) # indices for individuals in treatment group
  ind_A0 <- which(data$A==0) # indices for individuals in control group
  wt_A1.raw <- wt_A1[ind_A1]
  wt_A1 <- wt_A1.raw/sum(wt_A1.raw)
  wt_A0.raw <- wt_A0[ind_A0]
  wt_A0 <- wt_A0.raw/sum(wt_A0.raw)
  
  ## mediator column names
  Mnames <- colnames(data)[grep("M",colnames(data))]
  t <- length(Mnames) # number of distinct mediators
  ## observed mediator values
  Mobs <- data[, Mnames]
  
  ## direct and joint indirect effects ----------------------------------------
  #### joint mediators as a single mediator
  dat <- data.table(data[,c("id",Lnames)])
  setkey(dat)
  #### duplicate rows
  alevels <- dat[,as.data.table(cbind(
    "a.i"=1:3,
    "a0"=c(0,0,1),
    "a1"=c(0,1,1))),
    by=id]
  setkey(alevels)
  dat <- merge(dat,alevels,all.x=TRUE)
  setkey(dat)
  rm(alevels)
  
  #### sample *joint* mediator values
  SampleMs_joint <- function(mydt) {
    # mydt is a data.table for each id containing the relevant columns
    # e.g., mydt <- dat[id==1]
    mydt_mc <- mydt[rep(1:nrow(mydt),each=mc_draws)]
    setkey(mydt_mc)
    
    M.a1 <- data.frame(matrix(NA,nrow=nrow(mydt_mc),ncol=t))
    colnames(M.a1) <- Mnames
    M.a1[mydt_mc$a1==0,] <- Mobs[
      sample(x=ind_A0,size=sum(mydt_mc$a1==0),replace=TRUE,prob=wt_A0),]
    M.a1[mydt_mc$a1==1,] <- Mobs[
      sample(x=ind_A1,size=sum(mydt_mc$a1==1),replace=TRUE,prob=wt_A1),]
    
    mydt_mc <- cbind(mydt_mc,M.a1)
    setkey(mydt_mc)
    mydt_mc[, "Y.a" := PredictY(onedat=mydt_mc[, c("a0",Lnames,Mnames),
                                               with=FALSE])]
    setkey(mydt_mc)
    
    ## average over MC draws
    mydt_mc_means <- mydt_mc[,lapply(.SD, mean),by=a.i, .SDcols="Y.a"]
    setkey(mydt_mc_means)
    return(mydt_mc_means[,Y.a])
  }
  dat[, "Y.a" := SampleMs_joint(.SD), by=id]
  setkey(dat)
  
  # average over all individuals
  dat <- dat[, lapply(.SD,mean), by=list(a.i,a0,a1), .SDcols="Y.a"]
  logit <- function(x) log(x/(1-x))
  gamma <- c(logit(dat[3, Y.a])-logit(dat[2, Y.a]),
             logit(dat[2, Y.a])-logit(dat[1, Y.a]))

  names(gamma) <- paste0("g",0:1)
  res[["gamma"]] <- gamma
  rm(gamma,dat)

  ## mediator-specific indirect effects ---------------------------------------
  dat <- data.table(data[,c("id",Lnames)])
  setkey(dat)
  alevels <- dat[,as.data.table(Dupdata(t)),by=id]
  setkey(alevels)
  dat <- merge(dat,alevels,all.x=TRUE)
  setkey(dat)
  rm(alevels)
  
  #### sample *marginal* mediator values
  SampleMs_marginal <- function(mydt) {
    # mydt is a data.table for each id containing the relevant columns
    # e.g., mydt <- dat[id==1]
    mydt_mc <- mydt[rep(1:nrow(mydt),each=mc_draws)]
    setkey(mydt_mc)
    
    M.as <- lapply(1:t, function(s) {
      a_s <- unlist(mydt_mc[,paste0("a",s),with=FALSE])
      Ms.as <- rep(NA,length(a_s))
      Ms.as[a_s==0] <- Mobs[
        sample(x=ind_A0,size=sum(a_s==0),replace=TRUE,prob=wt_A0),s]
      Ms.as[a_s==1] <- Mobs[
        sample(x=ind_A1,size=sum(a_s==1),replace=TRUE,prob=wt_A1),s]
      return(Ms.as)
    })
    Ma <- data.table(do.call(cbind,M.as))
    setnames(Ma,Mnames) # do not order by setting key!!
    mydt_mc <- cbind(mydt_mc,Ma)
    setkey(mydt_mc)
    mydt_mc[, "Y.a" := PredictY(onedat=mydt_mc[, c("a0",Lnames,Mnames), 
                                               with=FALSE])]
    setkey(mydt_mc)
    
    ## average over MC draws
    mydt_mc_means <- mydt_mc[,lapply(.SD, mean),by=a.i, .SDcols="Y.a"]
    setkey(mydt_mc_means)
    return(mydt_mc_means[,Y.a])
  }
  dat[, "Y.a" := SampleMs_marginal(.SD), by=id]
  setkey(dat)
  
  # average over all individuals
  dat <- dat[, lapply(.SD,mean), by=c("a.i",paste0("a",0:t)), .SDcols="Y.a"]
  theta <- rep(NA,t)
  for (s in 1:t) {
    theta[s] <- logit(dat[s+1, Y.a])-logit(dat[1, Y.a])
  }
  names(theta) <- paste0("t",1:t)
  res[["theta"]] <- theta
  res[["mu"]] <- logit(dat[nrow(dat), Y.a])-logit(dat[1, Y.a])
  res$mu <- as.numeric(res$gamma["g1"])-res$mu
  rm(theta,dat)
  
  # weights
  res[["we.raw"]] <- c(summary(wt_A0.raw),summary(wt_A1.raw))
  res[["we.scaled"]] <- c(summary(wt_A0),summary(wt_A1))
  return(res)
}