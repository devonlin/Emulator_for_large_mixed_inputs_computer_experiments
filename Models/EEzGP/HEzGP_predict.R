HEzGP_predict=function (X_new, model, MSE_on = 0) 
{
  if (!is(model, "HEzGP model")) {
    stop("    The 2nd input should be a model of class \"EzGP model\".")
  }
  if (length(MSE_on) != 1) {
    stop("    MSE_on should be a scalar flag. Set it to 1 to turn it \"on\".")
  }
  X <- model$data$X
  Y <- model$data$Y
  p <- model$data$p
  q <- model$data$q
  m <- model$data$m
  n <- nrow(model$data$X)
  tau <- model$data$tau
  parv <- model$param
  covm <- cov_Hm(X, p, q, m, n, parv, tau,models = 2)
  y = as.matrix(Y)
  Tm = try(chol(covm), silent = TRUE)
  if ("try-error" %in% class(Tm)) {
    return(NULL)
  }
  m1 = as.matrix(c(rep(1, n)))
  invT = backsolve(Tm, diag(dim(Tm)[1]))
  invc = invT %*% t(invT)
  mu = as.numeric(1/sum(invc) * (t(m1) %*% invc %*% y))
  npar = 1 + q + p + p 
  psum <- function(x1, x2, par2) {
    return(sum(-par2 * (x1 - x2)^2))
  }
  covx <- function(w1, w2, parv) {
    par1 = parv[1:(q + 1)]
    par2 = parv[(q + 2):(q + 1 + p)]
    par3 = parv[(q + 2 + p):npar]
    x1 = w1[1:p]
    z1 = w1[(p + 1):(p + q)]
    x2 = w2[1:p]
    z2 = w2[(p + 1):(p + q)]
    res1 = par1[1] * exp(psum(x1, x2, par2))
    # for (i in 1:q) {
    #   if (z1[i] != z2[i]) {
    #     res1 = res1 + 0
    #   }
    #   else {
    #     l = z1[i]
    
    qual=c()
    for(i in 1:q){
      if ((z1[i] != z2[i] )){
        qual_ind=0
      }else{
        qual_ind=1
      }
      qual[i]= par1[i + 1]*qual_ind
    }
        res1 = res1 + sum(qual) * exp(psum(x1, x2, par3))
    #   }
    # }
    return(res1)
  }
  prey <- function(wn) {
    covv = matrix(0, n)
    for (i in 1:n) {
      if (sum(round(wn, 5) != round(X[i, ], 5)) > 0) {
        covv[i] = covx(wn, X[i, ], parv)
      }
      else {
        covv[i] = sum(parv[1:(q + 1)]) + tau
      }
    }
    gamma = as.matrix(covv)
    Y_hat = mu + t(gamma) %*% invc %*% (y - mu * m1)
    MSE = sum(parv[1:(q + 1)]) - (t(gamma) %*% invc %*% gamma)
    +((1 - t(m1) %*% invc %*% gamma)^2)/(t(m1) %*% invc %*% 
                                           m1)
    result <- c(Y_hat, MSE)
    return(result)
  }
  prediction <- list()
  if (is.vector(X_new) == TRUE) {
    p_all <- ncol(model$data$X)
    if (length(X_new) != p_all) {
      stop("    The dimensionality of X_new is not correct!")
    }
    value <- prey(X_new)
    if (MSE_on) {
      prediction$Y_hat <- value[1]
      prediction$MSE <- value[2]
    }
    else {
      prediction$Y_hat <- value[1]
    }
  }
  else if (is.matrix(X_new) == TRUE) {
    p_all <- ncol(model$data$X)
    if (ncol(X_new) != p_all) {
      stop("    The dimensionality of X_new is not correct!")
    }
    Y_hat = c()
    MSE = c()
    nn = nrow(X_new)
    for (i in 1:nn) {
      value = prey(X_new[i, ])
      Y_hat[i] = value[1]
      MSE[i] = value[2]
    }
    if (MSE_on) {
      prediction$Y_hat <- Y_hat
      prediction$MSE <- MSE
    }
    else {
      prediction$Y_hat <- Y_hat
    }
  }
  else {
    stop("    X_new must be a matrix or a vector")
  }
  return(prediction)
}
