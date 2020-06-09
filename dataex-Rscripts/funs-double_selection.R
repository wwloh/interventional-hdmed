DoubleSelection <- function(L,A,M,Y,fitY.family,fitM.family) {
  t <- ncol(M) # number of mediators
  L <- as.matrix(L)
  A <- as.numeric(A)
  M <- as.matrix(M)
  Y <- as.numeric(Y)
  
  # glmnet: Regression with LASSO to select predictors ==========================
  ## https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html#log
  ## continuous: use deviance/MSE as the criterion for 10-fold cross-validation
  ## binary: use misclassification error
  fitY.type.measure <- ifelse(fitY.family=="gaussian","deviance","class")
  fitM.type.measure <- ifelse(fitM.family=="gaussian","deviance","class")
  glmnet.control(mxit = 2e3) # increase max. iterations from default of 100
  res <- list()
  ptm=proc.time()[3]
  for (s in 1:t) {
    ## outcome model
    #### penalty factors to always include mediator Ms
    penfac <- rep(1,t)
    penfac[s] <- 0
    penfac <- c(rep(1,ncol(L)),0,penfac) # ordered as L, A, M
    refit <- TRUE
    while (refit) {
      fitY.stm <- tryCatch(system.time(
        fit.Y_Ms <- glmnet::cv.glmnet(x=as.matrix(cbind(L,A,M)),
                                      y=Y,
                                      family=fitY.family,
                                      penalty.factor=penfac,
                                      type.measure = fitY.type.measure)
      )[3], error=function(cond) return(NA))
      refit <- is.na(fitY.stm)
      if (refit==TRUE) {
        # penalize without cross-validation
        fit.Y_Ms <- glmnet::glmnet(x=as.matrix(cbind(L,A,M)),
                                   y=Y,
                                   family=fitY.family,
                                   penalty.factor=penfac)
        if (max(fit.Y_Ms$lambda) < Inf) {
          refit <- FALSE
        } else {
          # ignore all other mediators
          fit.Y_Ms <- glm(Y~.,family=fitY.family,
                          data=data.frame(cbind(L,A,M[,s])))
          refit <- FALSE
        }
      }
    }
    if (any(class(fit.Y_Ms)=="cv.glmnet")) {
      fit.Y_Ms.coef <- coef.glmnet(fit.Y_Ms, s = "lambda.min")[,1]
    } else if (any(class(fit.Y_Ms)=="glmnet")) {
      fit.Y_Ms.coef <- coef.glmnet(fit.Y_Ms, s = min(fit.Y_Ms$lambda))[,1]
    } else {
      fit.Y_Ms.coef <- coef(fit.Y_Ms)
      fit.Y_Ms.coef.Ms <- as.numeric(fit.Y_Ms.coef[length(fit.Y_Ms.coef)])
      fit.Y_Ms.coef.noM <- fit.Y_Ms.coef[-length(fit.Y_Ms.coef)]
      ## fill in zeroes for all other mediators
      fit.Y_Ms.coef.allM <- rep(0,ncol(M))
      fit.Y_Ms.coef.allM[s] <- fit.Y_Ms.coef.Ms
      names(fit.Y_Ms.coef.allM) <- paste0("M",1:ncol(M))
      fit.Y_Ms.coef <- c(fit.Y_Ms.coef.noM,fit.Y_Ms.coef.allM)
    }
    rm(fit.Y_Ms,penfac,fitY.stm)
    
    ## mediator model
    refit <- TRUE
    while (refit) {
      fitM.stm <- tryCatch(system.time(
        fit.Ms <- glmnet::cv.glmnet(x=as.matrix(cbind(L,A,M[,-s])),
                                    y=M[,s],
                                    family=fitM.family,
                                    type.measure = fitM.type.measure)
      )[3], error=function(cond) return(NA))
      refit <- is.na(fitM.stm)
    }
    fit.Ms.coef <- coef.glmnet(fit.Ms, s = "lambda.min")[,1]
    rm(fit.Ms,fitM.stm)
    
    # selected predictors in either outcome or mediator model ===================
    res[[s]] <- list("Y"=fit.Y_Ms.coef,"M"=fit.Ms.coef)
    cat(s, "; time taken (mins) =", round((proc.time()[3]-ptm)/60), 
        "; approx. left (mins) =", round((proc.time()[3]-ptm)*(t/s-1)/60),
        "\n")
  }
  return(res)
}


OneMCestimator_noMmodels_oneIEonly <- function(data,mc_draws,fitY.coef,Ms.i) {
  
  # sampling probabilities for each individual's mediator values ==============
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
  
  # mediator column names
  Mnames <- grep("M",colnames(data),value=TRUE)
  t <- length(Mnames) # number of distinct mediators
  ## observed mediator values
  Mobs <- data[, Mnames]
  
  # mediator-specific indirect effect via one mediator only ===================
  ss <- which(Ms.i==Mnames)
  dat <- data.table(data[,1:(which(names(data)=="A")-1)])
  setkey(dat)
  #### duplicate rows
  dupas <- matrix(0,nrow=2,ncol=t) # hypothetical treatment levels for Ms
  dupas[2,ss] <- 1 # only treatment for Ms
  dupas <- cbind(0,dupas) # a0
  colnames(dupas) <- paste0("a",0:t)
  dupas <- cbind("a.i"=1:2,dupas)
  
  alevels <- dat[,as.data.table(dupas),by=id]
  setkey(alevels)
  dat <- merge(dat,alevels,all.x=TRUE)
  setkey(dat)
  rm(alevels)
  
  #### sample *marginal* mediator values
  mydt_mc <- dat[rep(1:nrow(dat),each=mc_draws)]
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
  Mtilde <- do.call(cbind,M.as); rm(M.as)
  colnames(Mtilde) <- Mnames
  mydt_mc <- cbind(mydt_mc,Mtilde); rm(Mtilde)
  setkey(mydt_mc)
  
  ## predict potential outcomes using sampled mediator values =================
  setnames(mydt_mc,"a0","A")
  # estimated coefficients in outcome model
  Ya <- (as.matrix(mydt_mc[,names(fitY.coef),with=FALSE]) %*% fitY.coef)[,1]
  Expit <- function(y.logit) {
    if (exp(y.logit)==Inf) {
      return( 1.0 )
    } else {
      return( exp(y.logit)/(1+exp(y.logit)) )
    }
  }
  Ya <- sapply(Ya, Expit)
  mydt_mc[, "Y.a" := Ya]
  setkey(mydt_mc)
  
  ## average over MC draws and all individuals
  dat <- mydt_mc[, lapply(.SD,mean), by=c("a.i","A",paste0("a",1:t)),
                 .SDcols="Y.a"]
  logit <- function(x) {
    if (1-x<.Machine$double.eps) {
      log(1/.Machine$double.eps)
    } else if (x<.Machine$double.eps) {
      log(.Machine$double.eps)
    } else {
      log(x/(1-x))
    }
  }
  theta <- logit(dat[2, Y.a])-logit(dat[1, Y.a])
  names(theta) <- paste0("t",ss)
  
  return(theta)
}
