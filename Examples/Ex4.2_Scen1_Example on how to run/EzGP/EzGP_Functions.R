#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#EzGP Functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#covariance EzGP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cov_EzGP <- function(W1,W2, p, q, m, parv, tau = 0) {
  # W1 : matrix n1 x (p+q) 
  # W2 : matrix n2 x  (p+q) 
  npar = 1 + q + p + p*sum(m)
  ## a help function
  psum <- function(x1,x2, par2)
  {
    return(sum(-par2*(x1-x2)^2))
  }
  n1 <- nrow(W1)
  n2 <- nrow(W2)
  # calculating covariance between two inputs w1 and w2
  covx_R <- function(w1,w2, parv){
    #variance parameter sigma^2
    par1 = parv[1:(q+1)]
    #correlation parameter in G0
    par2 = parv[(q+2):(q+1+p)]
    #correlation parameter in G1 to Gq
    par3 = parv[(q+2+p): npar]
    x1 = w1[1:p]
    z1 = w1[(p+1):(p+q)]
    x2 = w2[1:p]
    z2 = w2[(p+1):(p+q)]
    res1 = par1[1]*exp(psum(x1,x2,par2))
    
    # if (models == 0){
    #   
    for (i in 1:q){
      if(z1[i] != z2[i]){
        res1 = res1+0
      }
      else{
        l = z1[i]
        res1 = res1 + par1[i+1]*exp(psum(x1,x2, par3[(sum(m[1:i])*p - m[i]*p + (l-1)*p + 1) : (sum(m[1:i])*p - m[i]*p + (l-1)*p + p)]))
      }
    }
    return(res1)
  }
  
  
  if (n1==n2 & (sum(W1)-sum(W2))==0){
    X=W1
    covx.m <- function(w.12, parv){
      return( covx_R(w1 = w.12[1:(p+q)], w2 = w.12[(p+q+1):(2*p+2*q)], parv) )
    }
    
    #expand grid to avoid for looping
    rcoord <- cbind(
      rep(seq_len(n1 - 1), times = rev(seq_len(n1 - 1))),
      unlist(lapply(
        X = rev(seq_len(n1 - 1)),
        FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = n1)))
    
    #covariance matrix
    covm = matrix(0,n1,n2)
    # first compute the vector of elements in covariance matrix
    Rtemp <- apply(cbind(X[rcoord[, 1], ], X[rcoord[, 2], ]), 1, FUN = covx.m, parv=parv)
    covm[rcoord] <- Rtemp
    covm <- covm + t(covm)
    diag(covm) <- sum(parv[1:(q+1)]) + tau
  }else{
    covm <- matrix(0, n1, n2)
    for(i in seq(1,n1)){
      for(j in seq(1,n2)){
        covm[i,j] <-covx_R(w1=W1[i,], w2=W2[j ,], parv)
      }
    }
  }
  return(covm)		
}


#2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#EzGP_fit_partioned  model to save current covariance
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##modeling
EzGP_fit_partioned <- function(X, Y, p, q, m, tau = 0, lb = "T", ub = "T",  x0 = "T",darg_est=NULL,
                         xtol_rel = 1.0e-5, maxeval = 100, algorithm = "NLOPT_LD_LBFGS"){
  
  message('     The model fitting process may take a while.')
  # function input checking
  if (missing(X)){
    stop('    X must be provided.')
  }
  if (!is.matrix(X) && !is.data.frame(X)){
    stop('    X must be a matrix or a data frame.')
  }
  
  # Number of training data
  n = nrow(X)
  
  if (missing(Y)){
    stop('    Y must be provided.')
  }
  if (!all(is.finite(Y)) || !is.numeric(Y)){
    stop('    All the elements of Y must be finite numbers.')
  }
  if (is.vector(Y) == TRUE) {
    y <- as.matrix(Y)
  } else if (is.matrix(Y) == TRUE){
    if (ncol(Y) == 1){
      y <- Y
    } else{
      stop('    The response value (i.e. y) for each observation should has only one.')
    }
  }
  
  if (n != nrow(y)){
    stop('    The number of rows (i.e., observations) in X and Y should match!')
  }
  
  if (!is.finite(p) || !is.numeric(p)){
    stop('    p must be a finite number')
  }
  
  if (!is.finite(q) || !is.numeric(q)){
    stop('    q must be a finite number')
  }
  
  if (ncol(X) != (p+q)){
    stop('    The number of quantitative and qualitative factors does not match the input training data')
  }
  
  if (is.vector(m) == FALSE){
    stop('    m must be a vetcor')
  }
  
  if (length(m) != q){
    stop('    m must contain the number of levels of each of the qualitative factors, i.e. m should be a vector with length q')
  }
  
  
  if (is.matrix(lb) == TRUE){
    stop('    lb must be "T" (default setting) or a vector of finite numbers.')
  } else if (is.vector(lb) == TRUE){
    if ((length(lb) == 1) && (lb != "T")){
      stop('    lb must be "T" (default setting) or a vector of finite numbers')
    } else if (length(lb) != 1){
      if (is.numeric(lb) == FALSE){
        stop('    lb must be "T" (default setting) or a vector of finite numbers')
      } else if (length(lb) != npar){
        stop('    The length of lb must match the number of the parameters in the model')
      }
    }
  }
  
  if (is.matrix(ub) == TRUE){
    stop('    ub must be "T" (default setting) or a vector of finite numbers.')
  } else if (is.vector(ub) == TRUE){
    if ((length(ub) == 1) && (ub != "T")){
      stop('    ub must be "T" (default setting) or a vector of finite numbers')
    } else if (length(ub) != 1){
      if (is.numeric(ub) == FALSE){
        stop('    ub must be "T" (default setting) or a vector of finite numbers')
      } else if (length(ub) != npar){
        stop('    The length of ub must match the number of the parameters in the model')
      }
    }
  }
  
  if (is.matrix(x0) == TRUE){
    stop('    x0 must be "T" (default setting) or a vector of finite numbers.')
  } else if (is.vector(x0) == TRUE){
    if ((length(x0) == 1) && (x0 != "T")){
      stop('    x0 must be "T" (default setting) or a vector of finite numbers')
    } else if (length(x0) != 1){
      if (is.numeric(x0) == FALSE){
        stop('    x0 must be "T" (default setting) or a vector of finite numbers')
      } else if (length(x0) != npar){
        stop('    The length of x0 must match the number of the parameters in the model')
      }
    }
  }
  
  # total number of parameters in the model
  npar = 1 + q + p + p*sum(m)
  
  eval_f_list <- function(parv) {
    ## a help function
    psum <- function(x1,x2, par2)
    {
      return(sum(-par2*(x1-x2)^2))
    }
    
    
    #check if all parameters are all positive
    if (min(parv) < 0)
    {
      return(list( "objective" = 10000000000,
                   "gradient" = rep(NA, npar)))
    }
    
    #expand grid to avoid for looping
    rcoord <- cbind(
      rep(seq_len(n - 1), times = rev(seq_len(n - 1))),
      unlist(lapply(
        X = rev(seq_len(n - 1)),
        FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = n)))
    
    #covariance matrix
    covm = cov_EzGP(W1=X,W2=X, p=p, q=q, m=m, parv=parv, tau = tau)
    
    Tm = try(chol(covm), silent=TRUE)
    #round(t(Tm)%*%Tm,2) == round(covm,2)
    if ('try-error' %in% class(Tm)) {
      return(list( "objective" = 10000000000,
                   "gradient" = rep(NA, npar)))
    }
    m1 = as.matrix(c(rep(1,n)))
    invT = backsolve(Tm, diag(dim(Tm)[1]))
    invc = invT%*%t(invT)
    
    #likelihood function value
    MLE_result = 2*sum(log(diag(Tm))) + (t(y) %*% invc %*% y) - 1/sum(invc)*(t(m1) %*% invc %*% y)^2
    
    #calcaulte analytical gradients
    # trend estimator
    mu = as.numeric(1/sum(invc)*(t(m1) %*% invc %*% y))
    
    ### derivative for variance parameter sigma^2_0
    grad_var0 <- function()
    {
      ### derative function of sigma_0 for wi, wj
      gradf_var <- function(w12)
      {
        x1 = w12[1:p]
        x2 = w12[(p+q+1):(2*p+q)]
        #correlation parameter in G0
        par2 = parv[(q+2):(q+1+p)]
        gx = exp(psum(x1,x2,par2))
        return(as.numeric(gx))
      }
      der_m = matrix(0,n,n)
      Rtemp_var1 <- apply(cbind(X[rcoord[, 1], ], X[rcoord[, 2], ]), 1, FUN = gradf_var)
      der_m[rcoord] <- Rtemp_var1
      der_m <- der_m + t(der_m)
      diag(der_m) = 1
      result = sum(diag(invc %*% der_m)) - t(y-mu) %*% invc %*% der_m %*% invc %*% (y-mu)
      return(result)
    }
    #grad_var0()
    
    ### derivative for variance parameter sigma^2_h h =1, ..., q
    grad_var <- function(h)
    {
      ### derative function of sigma_h for wi, wj
      gradf_var <- function(w12,h)
      {
        x1 = w12[1:p]
        z1 = w12[(p+1):(p+q)]
        x2 = w12[(p+q+1):(p+q+p)]
        z2 = w12[(p+q+p+1):(2*p+2*q)]
        if(z1[h] != z2[h]){
          return(0)
        }
        else{
          l = as.numeric(z1[h])
          par3 = parv[(q+2+p): npar]
          gx = exp(psum(x1,x2, par3[(sum(m[1:h])*p - m[h]*p + (l-1)*p + 1) : (sum(m[1:h])*p - m[h]*p + (l-1)*p + p)]))
          return(as.numeric(gx))
        }
      }
      der_m = matrix(0,n,n)
      Rtemp_var1 <- apply(cbind(X[rcoord[, 1], ], X[rcoord[, 2], ]), 1, FUN = gradf_var, h=h)
      der_m[rcoord] <- Rtemp_var1
      der_m <- der_m + t(der_m)
      diag(der_m) = 1
      result = sum(diag(invc %*% der_m)) - t(y-mu) %*% invc %*% der_m %*% invc %*% (y-mu)
      return(result)
    }
    
    ### derivative for correlation parameter theta_(0)_s for s=1,...,p
    grad_cor0 <- function(s)
    {
      ### derative function of theta_0_s for wi, wj
      gradf_cor0 <- function(w12,s)
      {
        x1 = w12[1:p]
        #z1 = w12[(p+1):(p+q)]
        x2 = w12[(p+q+1):(p+q+p)]
        #z2 = w12[(p+q+p+1):(2*p+2*q)]
        #correlation parameter in G0
        par2 = parv[(q+2):(q+1+p)]
        gx = -parv[1] * (x1[s] - x2[s])^2 * exp(psum(x1,x2,par2))
        return(as.numeric(gx))
      }
      der_m = matrix(0,n,n)
      Rtemp_cor <- apply(cbind(X[rcoord[, 1], ], X[rcoord[, 2], ]), 1, FUN = gradf_cor0, s=s)
      der_m[rcoord] <- Rtemp_cor
      der_m <- der_m + t(der_m)
      diag(der_m) = 0
      result = sum(diag(invc %*% der_m)) - t(y-mu) %*% invc %*% der_m %*% invc %*% (y-mu)
      return(result)
    }
    
    ### derivative for correlation parameter theta_(h)_l_s for h=1,...,q, l=1,...,m_h, s=1,...,p
    grad_cor <- function(h,l,s)
    {
      ### derative function of theta_0_s for wi, wj
      gradf_corhs <- function(w12,h,l,s)
      {
        x1 = w12[1:p]
        z1 = w12[(p+1):(p+q)]
        x2 = w12[(p+q+1):(p+q+p)]
        z2 = w12[(p+q+p+1):(2*p+2*q)]
        if((z1[h] != l) | (z2[h] != l)){
          return(0)
        }
        else{
          #variance parameter sigma^2
          par1 = parv[1:(q+1)]
          #correlation parameter in G1 to Gq
          par3 = parv[(q+2+p): npar]
          gx = -par1[h+1] * (x1[s] - x2[s])^2 * exp(psum(x1,x2, par3[(sum(m[1:h])*p - m[h]*p + (l-1)*p + 1) : (sum(m[1:h])*p - m[h]*p + (l-1)*p + p)]))
          return(as.numeric(gx))
        }
      }
      der_m = matrix(0,n,n)
      Rtemp_cor <- apply(cbind(X[rcoord[, 1], ], X[rcoord[, 2], ]), 1, FUN = gradf_corhs, h=h, l=l, s=s)
      der_m[rcoord] <- Rtemp_cor
      der_m <- der_m + t(der_m)
      diag(der_m) = 0
      result = sum(diag(invc %*% der_m)) - t(y-mu) %*% invc %*% der_m %*% invc %*% (y-mu)
      return(result)
    }
    ###analytical gradients results
    grad_res = grad_var0()
    for (i in 1:q) {
      grad_res = c(grad_res, grad_var(i))
    }
    for (i in 1:p)
    {
      grad_res = c(grad_res, grad_cor0(i))
    }
    for (h0 in 1:q)
    {
      for (l0 in 1:m[h0])
      {
        for (s0 in 1:p)
        {
          grad_res = c(grad_res, grad_cor(h0,l0,s0))
          #print(grad_res)
        }
      }
    }
    return( list( "objective" = MLE_result,
                  "gradient" = grad_res))
  }
  if(length(darg_est)>1){
    # Set lower and upper bounds of control, and initial settings
    if (lb == "T") {
      lb_new <- rep(0.1, npar)
      lb_new=darg_est$min
    } else {
      lb_new <- lb
    }
    
    if (ub == "T") {
      
      ub_new <- c(rep(100, (q+1)), rep(10, (npar-q-1)))
      ub_new=darg_est$max
    } else {
      ub_new <- ub
    }
    
    if (x0 == "T") {
      x0_new = as.vector((lb_new + ub_new)/2)
      x0_new=darg_est$start
    } else {
      x0_new <- x0
    }
  }else{
    # Set lower and upper bounds of control, and initial settings
    if (lb == "T") {
      lb_new <- rep(0.1, npar)
    } else {
      lb_new <- lb
    }
    
    if (ub == "T") {
      
      ub_new <- c(rep(100, (q+1)), rep(10, (npar-q-1)))
    } else {
      ub_new <- ub
    }
    
    if (x0 == "T") {
      x0_new = as.vector((lb_new + ub_new)/2)
    } else {
      x0_new <- x0
    }
  }
  opts <- list("algorithm" = algorithm,
               "xtol_rel" = xtol_rel,
               "maxeval" = maxeval)
  
  res <- nloptr::nloptr( x0=x0_new,
                         eval_f=eval_f_list,
                         lb=lb_new,
                         ub=ub_new,
                         opts=opts)
  
  ##mle of parameters
  solpar = as.vector(res$solution)
  #Phi and Phi_inv
  covm <- cov_m(X, p, q, m, n, solpar, tau,models =0)
  y = as.matrix(Y)
  Tm = try(chol(covm), silent=TRUE)
  #round(t(Tm)%*%Tm,2) == round(covm,2)
  if ('try-error' %in% class(Tm)) {
    return(NULL)
  }
  invT = backsolve(Tm, diag(dim(Tm)[1]))
  invc = invT%*%t(invT)
  
  
  ## Save the fitted model
  model <- NULL
  model$param <- solpar
  model$data <- list('X' = X, 'Y' = Y, "p" = p, "q" = q, "m" = m, "tau" = tau,
                     "x0" = x0_new, "lb" = lb_new, "ub" = ub_new,
                     "xtol_rel" = xtol_rel, "maxeval" = maxeval)
  model$cov=list("Phi"=covm, "Phi_inv"=invc)
  class(model) <- 'EzGP model'
  
  return(model)
}


#3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#find ALC 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
update_EzGP_alc<- function(xm1, Xn, X,Y,n,parv,covm,invc,m){
  npar = 1 + q + p + p*sum(m)
  if(is.null(dim(xm1))) xm1 <- matrix(xm1, nrow = 1)
  #X_n+1 X_n+1
  Km1m1=cov_EzGP(W1=xm1,W2=xm1, p=p, q=q, m=m, parv=parv, tau = tau)
  
  #X_n+1 and Xn
  K_m.m1=cov_EzGP(W1=X,W2=xm1, p=p, q=q, m=m, parv=parv, tau = tau)
  
  #Xref and  X
  Kmn=cov_EzGP(W1=X,W2=Xn, p=p, q=q, m=m, parv=parv, tau = tau)
  
  #Xref X_n+1
  Km1n=cov_EzGP(W1=xm1,W2=Xn, p=p, q=q, m=m, parv=parv, tau = tau)
  
  gam <- invc %*% K_m.m1
  #vj+1
  rho <- Km1m1 - (t(K_m.m1) %*% gam) 
  
  GG=gam %*% t(gam) / as.numeric(rho)
  
  ALC_eq8=t(Kmn)%*%GG%*%Kmn+2*t(Kmn)%*%(gam/ as.numeric(rho))%*%Km1n+(Km1n)^2/as.numeric(rho)
  return(ALC_eq8)
}

#4
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Update alc adaptively  same parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
update_EzGP_adaptive=function(xm1,La_train,X,Y,parv,covm,invc,n,m){
  npar = 1 + q + p + p*sum(m)
  if(is.null(dim(xm1))) xm1 <- matrix(xm1, nrow = 1)
  #X_n+1 X_n+1
  Km1m1=cov_EzGP(W1=xm1,W2=xm1, p=p, q=q, m=m, parv=parv, tau = tau)
  
  #X_n+1 and Xn
  K_m.m1 =cov_EzGP(W1=X,W2=xm1, p=p, q=q, m=m, parv=parv, tau = tau)
  
  gam <- invc %*% K_m.m1
  #vj+1
  rho <- Km1m1 - (t(K_m.m1) %*% gam) 
  #K(x) partioned
  Km1 <- rbind(cbind(covm, K_m.m1),cbind(t(K_m.m1), Km1m1))
  # #K(x)^-1
  Km1i <- rbind(cbind(invc + gam %*% t(gam) / as.numeric(rho), -gam/as.numeric(rho)), cbind(-t(gam)/as.numeric(rho), 1/as.numeric(rho)))
  X=La_train[,1:(p+q)]
  Y=La_train[,1+p+q]
  n=nrow(La_train)
  
  ## Save the fitted model
  model <- NULL
  model$param <- parv
  model$data <- list('X' = X, 'Y' =Y, "p" = p, "q" = q, "m" = m, "tau" = tau,'n'=n)
  model$cov=list("Phi"=Km1, "Phi_inv"=Km1i)
  
  class(model) <- 'EzGP model'
  return(model)
  
}


#5
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##initial paramters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getDs_mixed <- function(X, p=0.1, samp.size=1000)
{
  if(nrow(X) > samp.size) {
    i <- sample(1:nrow(X), samp.size)
    X <- X[i,]
  }
  D <- distance(X)
  D <- D[upper.tri(D)]
  D <- D[D > 0]
  dstart <- as.numeric(quantile(D, p=p))
  dmax <- max(D, na.rm=TRUE)
  dmin <- min(D, na.rm=TRUE)
  return(list(start=dstart, min=dmin, max=dmax))
}

#6
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#update loglikelihood 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
update_EzGP_loglike<- function(data, model, tau = 0, lb = "T", ub = "T",  x0 = "T",
                         xtol_rel = 1.0e-5, maxeval = 100, algorithm = "NLOPT_LD_LBFGS"){
  
  
  x_ll=data[,1:(p+q)]
  y_ll=data[,1+p+q]
  X <- model$data$X
  Y <- model$data$Y
  p <- model$data$p
  q <- model$data$q
  m <- model$data$m
  n <- nrow(model$data$X)
  tau <- model$data$tau
  parv <- model$param
  Phi_gg= model$cov$Phi
  Phi_gg_inv= model$cov$Phi_inv
  m1= model$cov$m1
  mu=model$cov$mu
  MLE_m=model$loglike
  ## a help function
  psum <- function(x1,x2, par2)
  {
    return(sum(-par2*(x1-x2)^2))
  }
  npar = 1 + q + p + p*sum(m)
  
  
  
  if(is.null(dim(x_ll))) x_ll <- matrix(x_ll, nrow = 1)
  ###add m+1 data
  Xm1 <- rbind(model$data$X,x_ll)
  Ym1 <- c(model$data$Y,y_ll)
  
  message('     The model fitting process may take a while.')
  # function input checking
  if (missing(Xm1)){
    stop('    X must be provided.')
  }
  if (!is.matrix(Xm1) && !is.data.frame(Xm1)){
    stop('    X must be a matrix or a data frame.')
  }
  
  # Number of training data
  n = nrow(Xm1)
  
  if (missing(Ym1)){
    stop('    Y must be provided.')
  }
  if (!all(is.finite(Ym1)) || !is.numeric(Ym1)){
    stop('    All the elements of Y must be finite numbers.')
  }
  if (is.vector(Ym1) == TRUE) {
    y <- as.matrix(Ym1)
  } else if (is.matrix(Ym1) == TRUE){
    if (ncol(Ym1) == 1){
      y <- Ym1
    } else{
      stop('    The response value (i.e. y) for each observation should has only one.')
    }
  }
  
  if (n != nrow(y)){
    stop('    The number of rows (i.e., observations) in X and Y should match!')
  }
  
  if (!is.finite(p) || !is.numeric(p)){
    stop('    p must be a finite number')
  }
  
  if (!is.finite(q) || !is.numeric(q)){
    stop('    q must be a finite number')
  }
  
  if (ncol(Xm1) != (p+q)){
    stop('    The number of quantitative and qualitative factors does not match the input training data')
  }
  
  if (is.vector(m) == FALSE){
    stop('    m must be a vetcor')
  }
  
  if (length(m) != q){
    stop('    m must contain the number of levels of each of the qualitative factors, i.e. m should be a vector with length q')
  }
  
  
  if (is.matrix(lb) == TRUE){
    stop('    lb must be "T" (default setting) or a vector of finite numbers.')
  } else if (is.vector(lb) == TRUE){
    if ((length(lb) == 1) && (lb != "T")){
      stop('    lb must be "T" (default setting) or a vector of finite numbers')
    } else if (length(lb) != 1){
      if (is.numeric(lb) == FALSE){
        stop('    lb must be "T" (default setting) or a vector of finite numbers')
      } else if (length(lb) != npar){
        stop('    The length of lb must match the number of the parameters in the model')
      }
    }
  }
  
  if (is.matrix(ub) == TRUE){
    stop('    ub must be "T" (default setting) or a vector of finite numbers.')
  } else if (is.vector(ub) == TRUE){
    if ((length(ub) == 1) && (ub != "T")){
      stop('    ub must be "T" (default setting) or a vector of finite numbers')
    } else if (length(ub) != 1){
      if (is.numeric(ub) == FALSE){
        stop('    ub must be "T" (default setting) or a vector of finite numbers')
      } else if (length(ub) != npar){
        stop('    The length of ub must match the number of the parameters in the model')
      }
    }
  }
  
  if (is.matrix(x0) == TRUE){
    stop('    x0 must be "T" (default setting) or a vector of finite numbers.')
  } else if (is.vector(x0) == TRUE){
    if ((length(x0) == 1) && (x0 != "T")){
      stop('    x0 must be "T" (default setting) or a vector of finite numbers')
    } else if (length(x0) != 1){
      if (is.numeric(x0) == FALSE){
        stop('    x0 must be "T" (default setting) or a vector of finite numbers')
      } else if (length(x0) != npar){
        stop('    The length of x0 must match the number of the parameters in the model')
      }
    }
  }
  
  #expand grid to avoid for looping
  rcoord <- cbind(
    rep(seq_len(n - 1), times = rev(seq_len(n - 1))),
    unlist(lapply(
      X = rev(seq_len(n - 1)),
      FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = n)))
  
  # total number of parameters in the model
  
  eval_f_list <- function(parv) {
    #print(parv)
    ## a help function
    psum <- function(x1,x2, par2)
    {
      return(sum(-par2*(x1-x2)^2))
    }
    
    
    if(is.null(dim(x_ll))) x_ll <- matrix(x_ll, nrow = 1)
    #A=Phi_gg
    Phi_gg=cov_EzGP(W1=X,W2=X, p=p, q=q, m=m, parv=parv, tau = tau)
    #B=Phi_gl
    Phi_gl=cov_EzGP(W1=X,W2=x_ll, p=p, q=q, m=m, parv=parv, tau = tau)
    #D=Phi_ll
    Phi_ll=cov_EzGP(W1=x_ll,W2=x_ll, p=p, q=q, m=m, parv=parv, tau = tau)
    
    
    #Phi_mm
    Phi_mm<- rbind(cbind(Phi_gg, Phi_gl),cbind(t(Phi_gl), Phi_ll))
    #D_block
    D_block=solve(Phi_ll-(t(Phi_gl))%*%solve(Phi_gg)%*%Phi_gl)
    
    Phi_mm_inv= rbind(cbind(((solve(Phi_gg))+(solve(Phi_gg))%*%Phi_gl%*%D_block%*%t(Phi_gl)%*%(solve(Phi_gg))),
                            (-(solve(Phi_gg))%*%Phi_gl%*%D_block)), 
                      cbind((-D_block%*%t(Phi_gl)%*%(solve(Phi_gg))),D_block))

    #check if all parameters are all positive
    if (min(parv) < 0)
    {
      return(list( "objective" = 10000000000,
                   "gradient" = rep(NA, npar)))
    }
    
    Tm = try(chol(Phi_mm), silent=TRUE)
    #round(t(Tm)%*%Tm,2) == round(covm,2)
    if ('try-error' %in% class(Tm)) {
      return(list( "objective" = 10000000000,
                   "gradient" = rep(NA, npar)))
    }
    invT = backsolve(Tm, diag(dim(Tm)[1]))
    
    
    m1 = as.matrix(c(rep(1,n)))
    
    #likelihood function value
    MLE_result =2*sum(log(diag(Tm))) + (t(y) %*% Phi_mm_inv %*% y) - 1/sum(Phi_mm_inv)*(t(m1) %*% Phi_mm_inv %*% y)^2
    #print(MLE_result)
    #calcaulte analytical gradients
    # trend estimator
    mu = as.numeric(1/sum(Phi_mm_inv)*(t(m1) %*% Phi_mm_inv %*% y))
    ### derivative for variance parameter sigma^2_0
    grad_var0 <- function()
    {
      ### derative function of sigma_0 for wi, wj
      gradf_var <- function(w12)
      {
        x1 = w12[1:p]
        x2 = w12[(p+q+1):(2*p+q)]
        #correlation parameter in G0
        par2 = parv[(q+2):(q+1+p)]
        gx = exp(psum(x1,x2,par2))
        return(as.numeric(gx))
      }
      der_m = matrix(0,n,n)
      Rtemp_var1 <- apply(cbind(Xm1[rcoord[, 1], ], Xm1[rcoord[, 2], ]), 1, FUN = gradf_var)
      der_m[rcoord] <- Rtemp_var1
      der_m <- der_m + t(der_m)
      diag(der_m) = 1
      result = sum(diag(Phi_mm_inv %*% der_m)) - t(y-mu) %*% Phi_mm_inv %*% der_m %*% Phi_mm_inv %*% (y-mu)
      return(result)
    }
    #grad_var0()
    
    ### derivative for variance parameter sigma^2_h h =1, ..., q
    grad_var <- function(h)
    {
      ### derative function of sigma_h for wi, wj
      gradf_var <- function(w12,h)
      {
        x1 = w12[1:p]
        z1 = w12[(p+1):(p+q)]
        x2 = w12[(p+q+1):(p+q+p)]
        z2 = w12[(p+q+p+1):(2*p+2*q)]
        if(z1[h] != z2[h]){
          return(0)
        }
        else{
          l = as.numeric(z1[h])
          par3 = parv[(q+2+p): npar]
          gx = exp(psum(x1,x2, par3[(sum(m[1:h])*p - m[h]*p + (l-1)*p + 1) : (sum(m[1:h])*p - m[h]*p + (l-1)*p + p)]))
          return(as.numeric(gx))
        }
      }
      der_m = matrix(0,n,n)
      Rtemp_var1 <- apply(cbind(Xm1[rcoord[, 1], ], Xm1[rcoord[, 2], ]), 1, FUN = gradf_var, h=h)
      der_m[rcoord] <- Rtemp_var1
      der_m <- der_m + t(der_m)
      diag(der_m) = 1
      result = sum(diag(Phi_mm_inv %*% der_m)) - t(y-mu) %*% Phi_mm_inv %*% der_m %*% Phi_mm_inv %*% (y-mu)
      return(result)
    }
    
    ### derivative for correlation parameter theta_(0)_s for s=1,...,p
    grad_cor0 <- function(s)
    {
      ### derative function of theta_0_s for wi, wj
      gradf_cor0 <- function(w12,s)
      {
        x1 = w12[1:p]
        #z1 = w12[(p+1):(p+q)]
        x2 = w12[(p+q+1):(p+q+p)]
        #z2 = w12[(p+q+p+1):(2*p+2*q)]
        #correlation parameter in G0
        par2 = parv[(q+2):(q+1+p)]
        gx = -parv[1] * (x1[s] - x2[s])^2 * exp(psum(x1,x2,par2))
        return(as.numeric(gx))
      }
      der_m = matrix(0,n,n)
      Rtemp_cor <- apply(cbind(Xm1[rcoord[, 1], ], Xm1[rcoord[, 2], ]), 1, FUN = gradf_cor0, s=s)
      der_m[rcoord] <- Rtemp_cor
      der_m <- der_m + t(der_m)
      diag(der_m) = 0
      result = sum(diag(Phi_mm_inv %*% der_m)) - t(y-mu) %*% Phi_mm_inv %*% der_m %*% Phi_mm_inv %*% (y-mu)
      return(result)
    }
    
    ### derivative for correlation parameter theta_(h)_l_s for h=1,...,q, l=1,...,m_h, s=1,...,p
    grad_cor <- function(h,l,s)
    {
      ### derative function of theta_0_s for wi, wj
      gradf_corhs <- function(w12,h,l,s)
      {
        x1 = w12[1:p]
        z1 = w12[(p+1):(p+q)]
        x2 = w12[(p+q+1):(p+q+p)]
        z2 = w12[(p+q+p+1):(2*p+2*q)]
        if((z1[h] != l) | (z2[h] != l)){
          return(0)
        }
        else{
          #variance parameter sigma^2
          par1 = parv[1:(q+1)]
          #correlation parameter in G1 to Gq
          par3 = parv[(q+2+p): npar]
          gx = -par1[h+1] * (x1[s] - x2[s])^2 * exp(psum(x1,x2, par3[(sum(m[1:h])*p - m[h]*p + (l-1)*p + 1) : (sum(m[1:h])*p - m[h]*p + (l-1)*p + p)]))
          return(as.numeric(gx))
        }
      }
      der_m = matrix(0,n,n)
      Rtemp_cor <- apply(cbind(Xm1[rcoord[, 1], ], Xm1[rcoord[, 2], ]), 1, FUN = gradf_corhs, h=h, l=l, s=s)
      der_m[rcoord] <- Rtemp_cor
      der_m <- der_m + t(der_m)
      diag(der_m) = 0
      result = sum(diag(Phi_mm_inv %*% der_m)) - t(y-mu) %*% Phi_mm_inv %*% der_m %*% Phi_mm_inv %*% (y-mu)
      return(result)
    }
    ###analytical gradients results
    grad_res = grad_var0()
    for (i in 1:q) {
      grad_res = c(grad_res, grad_var(i))
    }
    for (i in 1:p)
    {
      grad_res = c(grad_res, grad_cor0(i))
    }
    for (h0 in 1:q)
    {
      for (l0 in 1:m[h0])
      {
        for (s0 in 1:p)
        {
          grad_res = c(grad_res, grad_cor(h0,l0,s0))
          # print(grad_res)
        }
      }
    }
    return( list( "objective" = MLE_result,
                  "gradient" = grad_res))
  }
  
  
  
  # Set lower and upper bounds of control, and initial settings
  if (lb == "T") {
    lb_new <- rep(0.1, npar)
  } else {
    lb_new <- lb
  }
  
  if (ub == "T") {
    ub_new <- c(rep(100, (q+1)), rep(10, (npar-q-1)))
  } else {
    ub_new <- ub
  }
  
  if (x0 == "T") {
    x0_new = parv
  } else {
    x0_new <- x0
  }
  
  opts <- list("algorithm" = algorithm,
               "xtol_rel" = xtol_rel,
               "maxeval" = maxeval)
  
  res <- nloptr::nloptr(x0=x0_new,
                        eval_f=eval_f_list,
                        lb=lb_new,
                        ub=ub_new,
                        opts=opts)
  
  ##mle of parameters
  solpar = as.vector(res$solution)
  
  model<- NULL
  model$param <- solpar
  model$data <- list('X' = Xm1, 'Y' = Ym1, "p" = p, "q" = q, "m" = m, "tau" = tau,
                     "x0" = x0_new, "lb" = lb_new, "ub" = ub_new,
                     "xtol_rel" = xtol_rel, "maxeval" = maxeval )

  class(model) <- 'EzGP model'
  
  return(model)
}

