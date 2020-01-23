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

OneMCestimator_noMmodels <- function(data,mc_draws=1e3, y_cont) {
  # data=OneData();mc_draws=2;y_cont=TRUE
  res <- list()

  ## fit outcome model --------------------------------------------------------
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
  
  ## mediator column names
  Mnames <- colnames(data)[grep("M",colnames(data))]
  t <- length(Mnames) # number of distinct mediators
  ## observed mediator values in each treatment group
  Mobs <- list(data[data$A==0, Mnames],data[data$A==1, Mnames])
  n_A0 <- sum(data$A==0)
  n_A1 <- sum(data$A==1)
  
  ## direct and joint indirect effects ----------------------------------------
  #### joint mediators as a single mediator
  dat <- data.table(data[,c("id","L")]) #### keep only id, L to reduce memory
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
    M.a1[mydt_mc$a1==0,] <- Mobs[[1]][
      sample(1:n_A0, size=sum(mydt_mc$a1==0), replace=TRUE),]
    M.a1[mydt_mc$a1==1,] <- Mobs[[2]][
      sample(1:n_A1, size=sum(mydt_mc$a1==1), replace=TRUE),]
    
    mydt_mc <- cbind(mydt_mc,M.a1)
    setkey(mydt_mc)
    mydt_mc[, "Y.a" := PredictY(onedat=mydt_mc[, c("a0","L",Mnames), with=FALSE])]
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
  if (y_cont==TRUE) {
    gamma <- c(dat[3, Y.a]-dat[2, Y.a],dat[2, Y.a]-dat[1, Y.a])
  } else {
    logit <- function(x) log(x/(1-x))
    gamma <- c(logit(dat[3, Y.a])-logit(dat[2, Y.a]),
               logit(dat[2, Y.a])-logit(dat[1, Y.a]))
  }
  names(gamma) <- paste0("g",0:1)
  res[["gamma"]] <- gamma
  rm(gamma,dat)

  ## mediator-specific indirect effects ---------------------------------------
  dat <- data.table(data[,c("id","L")]) #### keep only id, L to reduce memory
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
      Ms.as[a_s==0] <- 
        Mobs[[1]][sample(1:n_A0, size=sum(a_s==0), replace=TRUE),s]
      Ms.as[a_s==1] <- 
        Mobs[[2]][sample(1:n_A1, size=sum(a_s==1), replace=TRUE),s]
      return(Ms.as)
    })
    Ma <- data.table(do.call(cbind,M.as))
    setnames(Ma,Mnames) # do not order by setting key!!
    mydt_mc <- cbind(mydt_mc,Ma)
    setkey(mydt_mc)
    mydt_mc[, "Y.a" := PredictY(onedat=mydt_mc[, c("a0","L",Mnames), with=FALSE])]
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
    if (y_cont==TRUE) {
      theta[s] <- dat[s+1, Y.a]-dat[1, Y.a] 
    } else {
      theta[s] <- logit(dat[s+1, Y.a])-logit(dat[1, Y.a])
    }
  }
  names(theta) <- paste0("t",1:t)
  res[["theta"]] <- theta
  if (y_cont==TRUE) {
    res[["mu"]] <- dat[nrow(dat), Y.a]-dat[1, Y.a]
  } else {
    res[["mu"]] <- logit(dat[nrow(dat), Y.a])-logit(dat[1, Y.a])
  }
  res$mu <- as.numeric(res$gamma["g1"])-res$mu
  
  rm(theta,dat)
  
  return(res)
}

OnePOCestimator <- function(data) {
  fitY <- glm(fitY.form, family = gaussian("identity"), data = data)
  ## mediator column names
  Mnames <- colnames(data)[grep("M",colnames(data))]
  t <- length(Mnames) # number of distinct mediators
  bpaths <- rep(NA,t)
  for (s in 1:t) {
    fitMs <- glm(as.formula(paste0("M",s,"~A+L")),
                 family = gaussian("identity"), data = data)
    bpaths[s] <- coef(fitMs)["A"]
  }
  return(c(coef(fitY)["A"],bpaths*coef(fitY)[Mnames]))
}