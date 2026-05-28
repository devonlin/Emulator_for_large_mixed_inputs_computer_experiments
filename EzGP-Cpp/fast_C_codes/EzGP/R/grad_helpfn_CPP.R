#1

### derivative for variance parameter sigma^2_0
grad_var0 <- function()
{
  gradf_var.m <- function(w.12, parv){
    return( gradf_var(w1= w.12[1:(p+q)], w2=w.12[(p+q+1):(2*p+2*q)], parv, q=q, p=p))
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
grad_var_1 <- function(h)
{
  gradf_var1.m <- function(w.12, parv,h){
    return( gradf_var1(w1= w.12[1:(p+q)], w2=w.12[(p+q+1):(2*p+2*q)],  parv, q=q, p=p,h=h, m=m))
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
grad_cor0 <- function(s)
{
  gradf_cor0.m <- function(w.12, parv,s){
    return( gradf_cor0(w1= w.12[1:(p+q)], w2=w.12[(p+q+1):(2*p+2*q)],s=s, parv, q=q, p=p, m=m))
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
grad_cor <- function(h,l,s)
{
  
  gradf_corhs.m <- function(w.12, parv,h,l,s){
    return( gradf_corhs(w1= w.12[1:(p+q)], w2=w.12[(p+q+1):(2*p+2*q)],h=h, l=l, s=s, parv, q=q, p=p, m=m))
  }
  der_m = matrix(0,n,n)
  Rtemp_cor <- apply(cbind(Xm1[rcoord[, 1], ], Xm1[rcoord[, 2], ]), 1, FUN = gradf_corhs.m, h=h, l=l, s=s,parv=parv)
  der_m[rcoord] <- Rtemp_cor
  der_m <- der_m + t(der_m)
  diag(der_m) = 0
  result = sum(diag(Phi_mm_inv %*% der_m)) - t(y-mu) %*% Phi_mm_inv %*% der_m %*% Phi_mm_inv %*% (y-mu)
  return(result)
}

