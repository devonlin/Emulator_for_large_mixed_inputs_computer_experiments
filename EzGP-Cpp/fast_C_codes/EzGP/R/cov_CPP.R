cov_EzGP_cpp <- function(W1,W2, p, q, m, parv, tau = 0, models = 0) {
  
  # X1 : matrix n1 x d - containing training points
  # X2 : matrix n2 x d - containing test points
  
  n1 <- nrow(W1)
  n2 <- nrow(W2)

   if (n1==n2){
    X=W1
    covx.m <- function(w.12, parv){
      return( covx_wrapper(w1= w.12[1:(p+q)], w2=w.12[(p+q+1):(2*p+2*q)], parv, q=q, p=p, m=m))
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
    Rtemp <- apply(cbind(matrix(X[rcoord[, 1], ],ncol=(p+q)), matrix(X[rcoord[, 2], ],ncol=(p+q))), 1, FUN = covx.m, parv=parv)
    covm[rcoord] <- Rtemp
    covm <- covm + t(covm)
    diag(covm) <- sum(parv[1:(q+1)]) + tau
  }else{
    covm <- matrix(0, n1, n2)
    for(i in seq(1,n1)){
      for(j in seq(1,n2)){
        covm[i,j] <-covx_wrapper(w1=W1[i,], w2=W2[j ,], parv, q=q, p=p, m=m)
      }
    }
  }
  return(covm)		
}
