ConfounderSelection <- function(L,A,M,Y,fitY.family,addL,lambda.seq=NULL) {
  predY_names <- c("A",colnames(L),colnames(M))
  t <- ncol(M) # number of mediators
  L <- as.matrix(L)
  A <- as.numeric(A)
  M <- as.matrix(M)
  Y <- as.numeric(Y)
  X <- as.matrix(cbind(A,L,M))
  
  # glmnet: Regression with LASSO to select predictors ==========================
  ## https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html#log
  ## https://glmnet.stanford.edu/articles/glmnet.html#logistic-regression
  ## continuous: use deviance/MSE as the criterion for cross-validation
  ## binary: use misclassification error
  fitY.type.measure <- ifelse(fitY.family=="gaussian","deviance","class")
  glmnet.control(mxit = 2e3) # increase max. iterations from default of 100
  res <- list()
  ptm=proc.time()[3]
  for (s in 1:t) {
    ## outcome model
    #### penalty factors to always include mediator Ms
    penfac <- rep(1,t)
    penfac[s] <- 0
    if (addL) {
      penfac.L <- rep(0,ncol(L))
    } else {
      penfac.L <- rep(1,ncol(L))
    }
    penfac <- c(0,penfac.L,penfac) # ordered as A, L, M
    # https://stackoverflow.com/questions/23686067/default-lambda-sequence-in-glmnet-for-cross-validation
    if (is.null(lambda.seq)) {
      if (fitY.family=="binomial") {
        lambda.seq <- seq(from=1e-6,to=1,length.out=200)*log(t)/n    
      }
    }
    fitY.stm <- tryCatch(system.time(
      fit.Y_Ms <- glmnet::cv.glmnet(x=X,
                                    y=Y,
                                    lambda=lambda.seq,
                                    family=fitY.family,
                                    penalty.factor=penfac,
                                    type.measure = fitY.type.measure,
                                    standardize=FALSE,
                                    nfolds=min(length(Y),100))
    )[3], error=function(cond) return(NA))
    if (all(!is.na(fitY.stm))) {
      # relaxed LASSO: refit lasso using only selected predictors
      # https://stats.stackexchange.com/questions/82466/why-use-lasso-estimates-over-ols-estimates-on-the-lasso-identified-subset-of-var
      fit.Y_Ms.coef <- coef(fit.Y_Ms, s = "lambda.min")[,1]
      fit.Y_Ms.selected <- (fit.Y_Ms.coef != 0)[-1] # remove intercept
      penfac.selected <- penfac[fit.Y_Ms.selected]
      fit.Y_Ms.selected <- names(fit.Y_Ms.selected)[fit.Y_Ms.selected]
      fitY_refit <- any(penfac.selected!=0)
      if (fitY_refit) {
        ## any predictors that need not be forced into model
        fitY_refit.stm <- tryCatch(system.time(
          fitY_refit_Ms <- glmnet::cv.glmnet(x=X[,fit.Y_Ms.selected],
                                             y=Y,
                                             lambda=lambda.seq,
                                             family=fitY.family,
                                             penalty.factor=penfac.selected,
                                             type.measure = fitY.type.measure,
                                             standardize=FALSE,
                                             nfolds=min(length(Y),100))
        )[3], error=function(cond) return(NA))
        if (all(!is.na(fitY_refit.stm))) {
          rm(fit.Y_Ms.coef)
          fit.Y_Ms.coef <- coef(fitY_refit_Ms, s = "lambda.min")[,1]
        } else {
          fitY_refit <- FALSE
        }
      }
      if (!fitY_refit) {
        ## only predictors that need to be forced into model
        fit.Y_Ms.selected <- paste0("Y~",paste(fit.Y_Ms.selected,collapse="+"))
        if (fitY.family=="gaussian") {
          fit.Y_Ms.reest <- lm(fit.Y_Ms.selected,data=data.frame(X,Y))
        } else {
          fit.Y_Ms.reest <- glm(fit.Y_Ms.selected,data=data.frame(X,Y),
                                family=binomial("logit"))
        }
        fit.Y_Ms.coef <- fit.Y_Ms.reest$coef
      }
    } else {
      fit.Y_Ms.coef <- rep(NA,length(predY_names))
      names(fit.Y_Ms.coef) <- predY_names
    }
    # selected predictors
    res[[s]] <- fit.Y_Ms.coef
    cat("M",s,":",fit.Y_Ms.coef, "\n")
    cat("time taken (mins) =", round((proc.time()[3]-ptm)/60), 
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
  mtilde_mat <- mydt_mc[,grep(pattern="(Intercept)",x=names(fitY.coef),
                              invert=TRUE,value=TRUE), with=FALSE]
  mtilde_mat <- cbind("(Intercept)"=1,mtilde_mat)
  mtilde_mat <- mtilde_mat[,names(fitY.coef),with=FALSE]
  Ya <- (as.matrix(mtilde_mat) %*% fitY.coef)[,1]
  y_binary <- all(sort(unique(data$Y))==(0:1))
  if (y_binary) {
    Logit <- function(x) {
      if (1-x<.Machine$double.eps) {
        log(1/.Machine$double.eps)
      } else if (x<.Machine$double.eps) {
        log(.Machine$double.eps)
      } else {
        log(x/(1-x))
      }
    }
    Expit <- function(y.logit) {
      if (exp(y.logit)==Inf) {
        return( 1.0 )
      } else {
        return( exp(y.logit)/(1+exp(y.logit)) )
      }
    }
    Ya <- sapply(Ya, Expit)  
  }
  mydt_mc[, "Y.a" := Ya]
  setkey(mydt_mc)
  
  ## average over MC draws and all individuals
  dat <- mydt_mc[, lapply(.SD,mean), by=c("a.i","A",paste0("a",1:t)),
                 .SDcols="Y.a"]
  if (y_binary) {
    theta <- Logit(dat[2, Y.a])-Logit(dat[1, Y.a])
  } else {
    theta <- dat[2, Y.a]-dat[1, Y.a]
  }
  names(theta) <- paste0("t",ss)
  
  return(theta)
}
