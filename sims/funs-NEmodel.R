OneSteenEstimator <- function(data,y_cont,oracle_M1_par) {
  
  fitM1 <- glm(M1 ~ A + L, family = gaussian("identity"), data = data)
  if (y_cont) {
    fitY <- glm(fitY.form, family = gaussian("identity"), data = data)  
  } else {
    fitY <- glm(fitY.form, family = binomial("logit"), data = data)
  }
  
  # expanded dataset
  n <- nrow(data)
  a_aux <- lapply(1:n, function(i) {
    A <- data$A[i]
    cbind(data[rep(i,4),],
          "replicate"=1:4,
          "a0"=rep(c(A,1-A),times=2),
          "a1"=rep(c(A,1-A),each=2),
          "a2"=rep(A,4)) # observed exposure A
  })
  dat <- data.frame(do.call(rbind,a_aux))
  
  # estimated weights
  W1.numer <- dnorm(
    x=dat$M1,
    mean=predict.glm(fitM1,newdata=data.frame("A"=dat$a1,"L"=dat$L),
                     type="response"),
    sd=sqrt(summary(fitM1)$dispersion))
  W1.denom <- dnorm(
    x=dat$M1,
    mean=predict.glm(fitM1,newdata=data.frame("A"=dat$a2,"L"=dat$L),
                     type="response"),
    sd=sqrt(summary(fitM1)$dispersion))
  W <- W1.numer/W1.denom
  dat <- cbind(dat,W); rm(W)
  
  # oracle weights
  alpha_coef <- oracle_M1_par$alpha
  W1.numer <- dnorm(
    x=dat$M1, 
    mean=alpha_coef["a0"] + alpha_coef["aA"]*dat$a1 + alpha_coef["al"]*dat$L,
    sd=oracle_M1_par$sigma)
  W1.denom <- dnorm(
    x=dat$M1, 
    mean=alpha_coef["a0"] + alpha_coef["aA"]*dat$a2 + alpha_coef["al"]*dat$L,
    sd=oracle_M1_par$sigma)
  W.oracle <- W1.numer/W1.denom
  dat <- cbind(dat,W.oracle); rm(W.oracle)
  
  dat[, "Yhat"] <- predict.glm(fitY,newdata=data.frame(
    "A"=dat$a0,"L"=dat$L,"M1"=dat$M1,"M2"=dat$M2),
    type="response")
  
  # fit NE models using either estimated or oracle weights
  if (y_cont) {
    fit_mod2 <- glm(Yhat ~ a0 + a1 * a2,
                    family = gaussian("identity"), 
                    data = dat, weights = W)
    fit_mod2.oracle <- glm(Yhat ~ a0 + a1 * a2,
                           family = gaussian("identity"),
                           data = dat, weights = W.oracle)
  } else {
    fit_mod2 <- glm(Yhat ~ a0 + a1 * a2,
                    family = binomial("logit"), 
                    data = dat, weights = W)
    fit_mod2.oracle <- glm(Yhat ~ a0 + a1 * a2,
                           family = binomial("logit"), 
                           data = dat, weights = W.oracle)
  }
  theta1 <- coef(fit_mod2)
  names(theta1) <- paste0("t",lapply(
    strsplit(names(theta1),"a"),paste0,collapse=""))
  theta1.oracle <- coef(fit_mod2.oracle)
  names(theta1.oracle) <- paste0("t",lapply(
    strsplit(names(theta1.oracle),"a"),paste0,collapse=""))
  
  return(list("ne"=theta1,"ne.o"=theta1.oracle))
}