#' @title The Fitting Function of \code{EzGP} Model
#'
#' @description Fits an Easy-to-Interpret Gaussian process (EzGP) model to a dataset as described in \code{reference 1}.
#'     The input variables are mixed (with both quantitative and qualitative inputs)
#'     The output variable is quantitative and scalar.
#'
#' @param X Matrix or data frame containing the inputs of training data. Each row represents the input setting of a data point and the columns are values of quantitative variables and qualitative variables.
#' @param Y Vector containing the outputs of training data points.
#' @param p Number of quantitative factors in the given dataset \code{X}.
#' @param q Number of qualitative factors in the given dataset \code{X}.
#' @param m A vector containing numbers of levels in the qualitative factors.
#' @param tau Nugget if needed. The default nugget is 0, otherwise it has to be a non-negative real value.
#' @param lb Vector with lower bounds of the parameter estimation. "T" for applying the default setting of lb (a vector of length number of parameters whose elements are all 0.1), otherwise one must provide a vector with length of the number of parameters.
#' @param ub Vector with upper bounds of the parameter estimation. "T" for applying the default setting of ub (a vector of length number of parameters whose first \code{q+1} elements are 100 and the rest \code{number of parameters - q - 1} elements are 10), otherwise one must provide a vector with the length being the number of parameters.
#' @param x0 Vector with starting values for the optimization. "T" for applying the default setting of x0 (a vector made by \code{(lb + ub)/2}), otherwise one must provide a vector with the length being the number of parameters.
#' @param xtol_rel Stopping criterion for relative change reached.
#' @param maxeval Termination condition by specifying a maximum number of function.
#' @param algorithm Optimization algorithm. See \href{https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/}{NLopt Algorithms} for more availiable algorithms.
#'
#' @import nloptr
#'
#' @return A model of class "EzGP model" list of the following items:
#' \itemize{
#' \item{\code{param}} {A list containing the estimated parameters}
#' \item{\code{data}} {A list containing the dataset and the information for fitting}
#' }
#'
#' @references
#' \enumerate{
#' \item "EzGP: Easy-to-Interpret Gaussian Process Models for Computer Experiments with Both Quantitative and Qualitative Factors", Qian Xiao, Abhyuday Mandal, C. Devon Lin, and Xinwei Deng (\doi{10.1137/19M1288462})
#' }
#'
#' @export
#'
#' @seealso
#' \code{\link[EzGP]{EzGP_predict}} to use the fitted EzGP model for prediction.\cr
#'
#' @examples
#' # Example with 3 quantitative and 3 qualitative variables (dataset included in the package):
#' #     Fit an EzGP model (with default settings), and then perform the prediction.
#' #     This example may run for a while.

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
  
  
  #Phi and Phi_inv
  #   covm =cov_EzGP_cpp(W1=X,W2=X, p=p, q=q, m=m, parv=solpar, tau = tau)
  #   y = as.matrix(Y)
  #   Tm = try(chol(covm), silent=TRUE)
  #   #round(t(Tm)%*%Tm,2) == round(covm,2)
  #   if ('try-error' %in% class(Tm)) {
  #     return(NULL)
  #   }
  #   m1 = as.matrix(c(rep(1,n)))
  #   invT = backsolve(Tm, diag(dim(Tm)[1]))
  #   invc = invT%*%t(invT)
  #   #invc=solve(covm)
  #   #print(invc)
  #   mu = as.numeric(1/sum(invc)*(t(m1) %*% invc %*% y))
  # #print(res$objective)
  #   Phi_log=res$objective-((t(y) %*% invc %*% y) -1/sum(invc)*(t(m1) %*% invc %*% y)^2)
  #   
  #   
  
  ## Save the fitted model
  model <- NULL
  model$param <- solpar
  model$data <- list('X' = X, 'Y' = Y, "p" = p, "q" = q, "m" = m, "tau" = tau,
                     "x0" = x0_new, "lb" = lb_new, "ub" = ub_new,
                     "xtol_rel" = xtol_rel, "maxeval" = maxeval )
  
  # model$cov=list("m1"=m1, "mu"=mu ,"Phi"=covm, "Phi_inv"=invc)
  # model$loglike=Phi_log
  
  class(model) <- 'EzGP model'
  
  return(model)
}

