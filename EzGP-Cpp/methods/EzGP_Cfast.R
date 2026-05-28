#-----------------------------------------------------------------------------
# packages
#-----------------------------------------------------------------------------
library(twinning)
library(dbscan)
library(mltools)
library(data.table)
library(Rcpp)

#-----------------------------------------------------------------------------
# Serial Job 
#-----------------------------------------------------------------------------
f=function(iteration){
#-----------------------------------------------------------------------------
# Initials
#-----------------------------------------------------------------------------
model_name="Twin-EzGP_Cfast"
print(model_name)
load("initial.RData")
#-----------------------------------------------------------------------------
#Read functions 
#-----------------------------------------------------------------------------
source("cov_CPP.R")
source("EzGP_fit_CPP.R")
source("Update_EzGP_CPP.R")
sourceCpp("grad_fn.cpp")
sourceCpp("covx.cpp")

#-----------------------------------------------------------------------------
# data
#-----------------------------------------------------------------------------
  filename= paste0("data/train_QQ","_",iteration,".csv",sep="")
  train_QQ=as.matrix(read.csv(filename))[,2:(p+q+2)]
  filename= paste0("data/test_QQ","_",iteration,".csv",sep="")
  test_QQ=as.matrix(read.csv(filename))[,2:(p+q+2)]
#------------------------------------------------------------------------------------
  
  
#------------------------------------------------------------------------------------
#EzGP codes in c++ modifications
#(As it is a serial job I added the function of EzGP that arent using the package in the function otherwise it does not work for serial jobs)
#------------------------------------------------------------------------------------
  updateloglike_CPP<- function(data, model, tau = 0, lb = "T", ub = "T",  x0 = "T",
                               xtol_rel = 1.0e-5, maxeval = 100, algorithm = "NLOPT_LD_LBFGS"){
    
    x_ll=data[,1:(p+q)]
    y_ll=data[,1+p+q]
    X <- model$data$X
    Y <- model$data$Y
    p <- model$data$p
    q <- model$data$q
    m <- model$data$m
    tau <- model$data$tau
    parv <- model$param
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
      ## a help function
      if(is.null(dim(x_ll))) x_ll <- matrix(x_ll, nrow = 1)
      #A=Phi_gg
      Phi_gg=cov_EzGP_cpp(W1=X,W2=X, p=p, q=q, m=m, parv=parv, tau = tau)
      #B=Phi_gl
      Phi_gl=cov_EzGP_cpp(W1=X,W2=x_ll, p=p, q=q, m=m, parv=parv, tau = tau)
      #D=Phi_ll
      Phi_ll=cov_EzGP_cpp(W1=x_ll,W2=x_ll, p=p, q=q, m=m, parv=parv, tau = tau)
      
      
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
      #calcaulte analytical gradients
      # trend estimator
      mu = as.numeric(1/sum(Phi_mm_inv)*(t(m1) %*% Phi_mm_inv %*% y))
      ### derivative for variance parameter sigma^2_0
      #1
      
      ### derivative for variance parameter sigma^2_0
      grad_var0_cpp_wrap <- function()
      {
        gradf_var.m <- function(w.12, parv){
          return( gradf_var_cpp(w1= w.12[1:(p+q)], w2=w.12[(p+q+1):(2*p+2*q)], parv, q=q, p=p))
        }
        
        der_m = matrix(0,n,n)
        Rtemp_var1 <- apply(cbind(Xm1[rcoord[, 1], ], Xm1[rcoord[, 2], ]), 1, FUN = gradf_var.m ,parv=parv)
        der_m[rcoord] <- Rtemp_var1
        der_m <- der_m + t(der_m)
        diag(der_m) = 1
        result = sum(diag(Phi_mm_inv %*% der_m)) - t(y-mu) %*% Phi_mm_inv %*% der_m %*% Phi_mm_inv %*% (y-mu)
        return(result)
      }
      
      #2
      
      ### derivative for variance parameter sigma^2_h h =1, ..., q
      grad_var_1_cpp_wrap <- function(h)
      {
        gradf_var1.m <- function(w.12, parv,h){
          return( gradf_var1_cpp(w1= w.12[1:(p+q)], w2=w.12[(p+q+1):(2*p+2*q)],  parv, q=q, p=p,h=h, m=m))
        }
        der_m = matrix(0,n,n)
        Rtemp_var1 <- apply(cbind(Xm1[rcoord[, 1], ], Xm1[rcoord[, 2], ]), 1, FUN = gradf_var1.m, parv=parv,h=h)
        der_m[rcoord] <- Rtemp_var1
        der_m <- der_m + t(der_m)
        diag(der_m) = 1
        result = sum(diag(Phi_mm_inv %*% der_m)) - t(y-mu) %*% Phi_mm_inv %*% der_m %*% Phi_mm_inv %*% (y-mu)
        return(result)
      }
      
      #3
      ### derivative for correlation parameter theta_(0)_s for s=1,...,p
      grad_cor0_cpp_wrap <- function(s)
      {
        gradf_cor0.m <- function(w.12, parv,s){
          return( gradf_cor0_cpp(w1= w.12[1:(p+q)], w2=w.12[(p+q+1):(2*p+2*q)],s=s, parv, q=q, p=p, m=m))
        }
        der_m = matrix(0,n,n)
        Rtemp_cor <- apply(cbind(Xm1[rcoord[, 1], ], Xm1[rcoord[, 2], ]), 1, FUN = gradf_cor0.m, s=s,parv=parv)
        der_m[rcoord] <- Rtemp_cor
        der_m <- der_m + t(der_m)
        diag(der_m) = 0
        result = sum(diag(Phi_mm_inv %*% der_m)) - t(y-mu) %*% Phi_mm_inv %*% der_m %*% Phi_mm_inv %*% (y-mu)
        return(result)
      }
      
      #4
      ### derivative for correlation parameter theta_(h)_l_s for h=1,...,q, l=1,...,m_h, s=1,...,p
      grad_cor_cpp_wrap <- function(h,l,s)
      {
        
        gradf_corhs.m <- function(w.12, parv,h,l,s){
          return( gradf_corhs_cpp(w1= w.12[1:(p+q)], w2=w.12[(p+q+1):(2*p+2*q)],h=h, l=l, s=s, parv, q=q, p=p, m=m))
        }
        der_m = matrix(0,n,n)
        Rtemp_cor <- apply(cbind(Xm1[rcoord[, 1], ], Xm1[rcoord[, 2], ]), 1, FUN = gradf_corhs.m, h=h, l=l, s=s,parv=parv)
        der_m[rcoord] <- Rtemp_cor
        der_m <- der_m + t(der_m)
        diag(der_m) = 0
        result = sum(diag(Phi_mm_inv %*% der_m)) - t(y-mu) %*% Phi_mm_inv %*% der_m %*% Phi_mm_inv %*% (y-mu)
        return(result)
      }
      ###analytical gradients results
      grad_res = grad_var0_cpp_wrap()
      for (i in 1:q) {
        grad_res = c(grad_res, grad_var_1_cpp_wrap(i))
      }
      for (i in 1:p)
      {
        grad_res = c(grad_res, grad_cor0_cpp_wrap(i))
      }
      for (h0 in 1:q)
      {
        for (l0 in 1:m[h0])
        {
          for (s0 in 1:p)
          {
            grad_res = c(grad_res, grad_cor_cpp_wrap(h0,l0,s0))
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
    
    
    ## Save the fitted model
    model<- NULL
    model$param <- solpar
    model$data <- list('X' = Xm1, 'Y' = Ym1, "p" = p, "q" = q, "m" = m, "tau" = tau,
                       "x0" = x0_new, "lb" = lb_new, "ub" = ub_new,
                       "xtol_rel" = xtol_rel, "maxeval" = maxeval )
    
    # model$cov=list("m1"=m1, "mu"=mu ,"Phi_mm"=Phi_mm, "Phi_mm_inv"=Phi_mm_inv)
    
    class(model) <- 'EzGP model'
    
    return(model)
  }
  
  
  
  EzGP_fit_CPP <- function(X, Y, p, q, m, tau = 0, lb = "T", ub = "T",  x0 = "T",
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
      covm =cov_EzGP_cpp(W1=X,W2=X, p=p, q=q, m=m, parv=parv, tau = tau)
      
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
      grad_var0_cpp_wrap <- function()
      {
        gradf_var.m <- function(w.12, parv){
          return( gradf_var_cpp(w1= w.12[1:(p+q)], w2=w.12[(p+q+1):(2*p+2*q)], parv, q=q, p=p))
        }
        
        der_m = matrix(0,n,n)
        Rtemp_var1 <- apply(cbind(X[rcoord[, 1], ], X[rcoord[, 2], ]), 1, FUN = gradf_var.m ,parv=parv)
        der_m[rcoord] <- Rtemp_var1
        der_m <- der_m + t(der_m)
        diag(der_m) = 1
        result = sum(diag(invc %*% der_m)) - t(y-mu) %*% invc %*% der_m %*% invc %*% (y-mu)
        return(result)
      }
      
      #2
      
      ### derivative for variance parameter sigma^2_h h =1, ..., q
      grad_var_1_cpp_wrap <- function(h)
      {
        gradf_var1.m <- function(w.12, parv,h){
          return( gradf_var1_cpp(w1= w.12[1:(p+q)], w2=w.12[(p+q+1):(2*p+2*q)],  parv, q=q, p=p,h=h, m=m))
        }
        der_m = matrix(0,n,n)
        Rtemp_var1 <- apply(cbind(X[rcoord[, 1], ], X[rcoord[, 2], ]), 1, FUN = gradf_var1.m, parv=parv,h=h)
        der_m[rcoord] <- Rtemp_var1
        der_m <- der_m + t(der_m)
        diag(der_m) = 1
        result = sum(diag(invc %*% der_m)) - t(y-mu) %*% invc %*% der_m %*% invc %*% (y-mu)
        return(result)
      }
      
      #3
      ### derivative for correlation parameter theta_(0)_s for s=1,...,p
      grad_cor0_cpp_wrap <- function(s)
      {
        gradf_cor0.m <- function(w.12, parv,s){
          return( gradf_cor0_cpp(w1= w.12[1:(p+q)], w2=w.12[(p+q+1):(2*p+2*q)],s=s, parv, q=q, p=p, m=m))
        }
        der_m = matrix(0,n,n)
        Rtemp_cor <- apply(cbind(X[rcoord[, 1], ], X[rcoord[, 2], ]), 1, FUN = gradf_cor0.m, s=s,parv=parv)
        der_m[rcoord] <- Rtemp_cor
        der_m <- der_m + t(der_m)
        diag(der_m) = 0
        
        result = sum(diag(invc %*% der_m)) - t(y-mu) %*% invc %*% der_m %*% invc %*% (y-mu)
        return(result)
      }
      
      #4
      ### derivative for correlation parameter theta_(h)_l_s for h=1,...,q, l=1,...,m_h, s=1,...,p
      grad_cor_cpp_wrap <- function(h,l,s)
      {
        
        gradf_corhs.m <- function(w.12, parv,h,l,s){
          return( gradf_corhs_cpp(w1= w.12[1:(p+q)], w2=w.12[(p+q+1):(2*p+2*q)],h=h, l=l, s=s, parv, q=q, p=p, m=m))
        }
        der_m = matrix(0,n,n)
        Rtemp_cor <- apply(cbind(X[rcoord[, 1], ], X[rcoord[, 2], ]), 1, FUN = gradf_corhs.m, h=h, l=l, s=s,parv=parv)
        der_m[rcoord] <- Rtemp_cor
        der_m <- der_m + t(der_m)
        diag(der_m) = 0
        result = sum(diag(invc %*% der_m)) - t(y-mu) %*% invc %*% der_m %*% invc %*% (y-mu)
        return(result)
      }
      ###analytical gradients results
      grad_res = grad_var0_cpp_wrap()
      for (i in 1:q) {
        grad_res = c(grad_res, grad_var_1_cpp_wrap(i))
      }
      for (i in 1:p)
      {
        grad_res = c(grad_res, grad_cor0_cpp_wrap(i))
      }
      for (h0 in 1:q)
      {
        for (l0 in 1:m[h0])
        {
          for (s0 in 1:p)
          {
            grad_res = c(grad_res, grad_cor_cpp_wrap(h0,l0,s0))
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
      x0_new = as.vector((lb_new + ub_new)/2)
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
    
    
    
    ## Save the fitted model
    model <- NULL
    model$param <- solpar
    model$data <- list('X' = X, 'Y' = Y, "p" = p, "q" = q, "m" = m, "tau" = tau,
                       "x0" = x0_new, "lb" = lb_new, "ub" = ub_new,
                       "xtol_rel" = xtol_rel, "maxeval" = maxeval )
    
    
    class(model) <- 'EzGP model'
    
    return(model)
  }
  
  
  
  
  
  #pred
  cov_m <- function(X, p, q, m, n, parv, tau = 0, models = 0){
    
    ##total number of parameters in the model
    if (models == 0){
      npar = 1 + q + p + p*sum(m)
    } else if (models == 1){
      npar = 1 + q + p + (sum(m)-q)
    }
    
    ## a help function
    psum <- function(x1,x2, par2)
    {
      return(sum(-par2*(x1-x2)^2))
    }
    
    # calculating covariance between two inputs w1 and w2
    covx <- function(w1,w2, parv){
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
      
      if (models == 0){
        
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
        
      } else if (models == 1){
        for (i in 1:q){
          if(z1[i] != z2[i]){
            res1 = res1+0
          }
          else{
            l = z1[i]
            if (l==1)
            {
              res1 = res1 + par1[i+1]*exp(psum(x1,x2,1))
            }
            else {
              res1 = res1 + par1[i+1]*exp(psum(x1,x2,par3[sum((m-1)[1:i]) - (m-1)[i] + l-1]))
            }
          }
        }
        return(res1)
      }
    }
    
    ### a modified version of covx where w.12 = (w1, w2)
    covx.m <- function(w.12, parv){
      return( covx(w1 = w.12[1:(p+q)], w2 = w.12[(p+q+1):(2*p+2*q)], parv) )
    }
    
    #expand grid to avoid for looping
    rcoord <- cbind(
      rep(seq_len(n - 1), times = rev(seq_len(n - 1))),
      unlist(lapply(
        X = rev(seq_len(n - 1)),
        FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = n)))
    
    #covariance matrix
    covm = matrix(0,n,n)
    # first compute the vector of elements in covariance matrix
    Rtemp <- apply(cbind(X[rcoord[, 1], ], X[rcoord[, 2], ]), 1, FUN = covx.m, parv=parv)
    covm[rcoord] <- Rtemp
    covm <- covm + t(covm)
    diag(covm) <- sum(parv[1:(q+1)]) + tau
    
    return(covm)
  }

  
  EzGP_predict <- function(X_new, model, MSE_on = 0){
    
    ## import model and check inputs
    if (!is(model,"EzGP model")){
      stop('    The 2nd input should be a model of class "EzGP model".')
    }
    if (length(MSE_on)!=1){
      stop('    MSE_on should be a scalar flag. Set it to 1 to turn it "on".')
    }
    
    X <- model$data$X
    Y <- model$data$Y
    p <- model$data$p
    q <- model$data$q
    m <- model$data$m
    n <- nrow(model$data$X)
    tau <- model$data$tau
    parv <- model$param
    covm <- cov_m(X, p, q, m, n, parv, tau)
    
    y = as.matrix(Y)
    Tm = try(chol(covm), silent=TRUE)
    #round(t(Tm)%*%Tm,2) == round(covm,2)
    if ('try-error' %in% class(Tm)) {
      return(NULL)
    }
    m1 = as.matrix(c(rep(1,n)))
    invT = backsolve(Tm, diag(dim(Tm)[1]))
    invc = invT%*%t(invT)
    mu = as.numeric(1/sum(invc)*(t(m1) %*% invc %*% y))
    
    
    ##total number of parameters in the model
    npar = 1 + q + p + p*sum(m)
    
    ## a help function
    psum <- function(x1,x2, par2)
    {
      return(sum(-par2*(x1-x2)^2))
    }
    
    # inner function
    covx <- function(w1,w2, parv){
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
    
    prey <- function(wn){
      covv = matrix(0,n)
      for(i in 1:n){
        if (sum(round(wn,5)!=round(X[i,],5)) > 0)
        {
          covv[i] = covx(wn,X[i,],parv)
        } else {
          covv[i] = sum(parv[1:(q+1)]) + tau
        }
      }
      gamma = as.matrix(covv)
      
      # Y_hat
      Y_hat = mu  + t(gamma) %*% invc %*% (y - mu * m1)
      # MSE
      MSE = sum(parv[1:(q+1)]) - (t(gamma) %*% invc %*% gamma)
      + ((1 - t(m1) %*% invc %*% gamma)^2)/(t(m1) %*% invc %*% m1)
      result <- c(Y_hat, MSE)
      
      return(result)
    }
    
    prediction <- list()
    if (is.vector(X_new) == TRUE){
      
      p_all <- ncol(model$data$X)
      if (length(X_new) != p_all){
        stop('    The dimensionality of X_new is not correct!')
      }
      
      value <- prey(X_new)
      if (MSE_on){
        prediction$Y_hat <- value[1]
        prediction$MSE <- value[2]
      } else{
        prediction$Y_hat <- value[1]
      }
      
    } else if (is.matrix(X_new) == TRUE){
      
      p_all <- ncol(model$data$X)
      if (ncol(X_new) != p_all){
        stop('    The dimensionality of X_new is not correct!')
      }
      
      Y_hat = c()
      MSE = c()
      nn = nrow(X_new)
      for (i in 1:nn){
        value = prey(X_new[i,])
        Y_hat[i] = value[1]
        MSE[i] = value[2]
      }
      
      if (MSE_on){
        prediction$Y_hat <- Y_hat
        prediction$MSE <- MSE
      } else{
        prediction$Y_hat <- Y_hat
      }
      
    } else{
      stop('    X_new must be a matrix or a vector')
    }
    
    return(prediction)
  }
  
  
  
  
  
  
#------------------------------------------------------------------------------------
#One-hot encoding
#------------------------------------------------------------------------------------
  train_QQ_1=data.frame(train_QQ)
  #
  if(q==1){
    names <- (p+1):(p+q)
    train_QQ_1[,names] <- sapply(train_QQ_1[,names], factor)
  }else{
    
    names <- (p+1):(p+q)
    train_QQ_1[,names] <- lapply(train_QQ_1[,names], factor)
  }
  
  train_QQ_encode=one_hot(as.data.table(train_QQ_1))
  train_QQ_encode=sapply(train_QQ_encode, as.numeric)
  
  id_test=sample_id[iteration]

test_QQ_encode_qual=apply(matrix(z_true_encode[id_test,],nrow=1),2,rep,each=test_size)
test_QQ_encode=cbind(test_QQ[,1:p],test_QQ_encode_qual,test_QQ[,p+q+1])
#twinning
set.seed(1)
twin1_indices= twin(train_QQ_encode, r=r)

twin1=train_QQ[twin1_indices,]
twin2=train_QQ[-twin1_indices,]
start_time <- Sys.time()
model_gg = EzGP_fit_CPP(twin1[,1:(p+q)],twin1[,(p+q+1)],p=p,q=q,m=m,tau=tau) 
end_time <- Sys.time()
time_global=as.numeric(end_time - start_time,unit = "mins")

  
Yhat=c()
time=c()
#------------------------------------------------------------------------------------
#Prediction with each test location
#------------------------------------------------------------------------------------
  
     for(iu in 1:nr){ 
       print(iu)
      #test location
      w_star=test_QQ_encode[iu,1:(ncol(test_QQ_encode))]
      #NN
      twin22=rbind((train_QQ_encode[-twin1_indices,]),matrix(w_star,nrow=1))
      nn <- kNN(twin22[,1:(ncol(train_QQ_encode)-1)], k = local_size)
      id <- nrow(twin22)
      twin_kd_tree=twin2[nn$id[id,],]
      if(is.matrix(twin_kd_tree)!=TRUE) twin_kd_tree= matrix(twin_kd_tree,1,ncol=length(twin_kd_tree))
      
      #local-global
      #------------------------------------------------------------------------------------
      # update paramters
      start_time <- Sys.time()
      model_loglike_update_CPP=updateloglike_CPP(data=twin_kd_tree, model=model_gg,tau = tau)
      end_time <- Sys.time()
      pred= EzGP_predict(test_QQ[iu,1:(p+q)], model_loglike_update_CPP, MSE_on = 0)
      Yhat_ind= pred$Y_hat#save
      Yhat=c(Yhat,Yhat_ind)
      time_ind=as.numeric(end_time - start_time,unit = "mins")
      time=c(time,time_ind)
     }
  
  #RMSE
  RMSE=sqrt(mean(sum((Yhat-test_QQ[1:nr,(p+q+1)])^2)))
#------------------------------------------------------------------------------------
#save Results
#------------------------------------------------------------------------------------
  
  filename= paste0("RMSE","_",model_name,"_",iteration,".rds",sep="")
  saveRDS(RMSE,file=filename)
  filename= paste0("Yhat","_",model_name,"_",iteration,".rds",sep="")
  saveRDS(Yhat,file=filename)
  filename= paste0("time","_",model_name,"_",iteration,".rds",sep="")
  saveRDS(time,file=filename)
  filename= paste0("time_global","_",model_name,"_",iteration,".rds",sep="")
  saveRDS(time_global,file=filename)
  
}
iteration=as.integer(commandArgs(trailingOnly = TRUE))
f(iteration)


  
  
  
  
