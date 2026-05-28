



cov_Hm <- function(X, p, q, m, n, parv, tau = 0, models = 0){

  ##total number of parameters in the model
  if (models == 0){
    npar = 1 + q + p + p*sum(m)
  } else if (models == 1){
    npar = 1 + q + p + (sum(m)-q)
  }else if (models == 2){
    npar = 1 + q + 2*p
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
    }else if (models == 2){
      
      # for (i in 1:q){
      #   if(z1[i] != z2[i]){
      #     res1 = res1+0
      #   }
      #   else{
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
      
      
      
          res1 = res1 + sum(qual)*exp(psum(x1,x2, par3))
        #}
      #}
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
#ccov_m_org=cov_m(X=X,p=p,q=q,m=m,n=nrow(X),parv=parv,tau=tau,models = 2)
