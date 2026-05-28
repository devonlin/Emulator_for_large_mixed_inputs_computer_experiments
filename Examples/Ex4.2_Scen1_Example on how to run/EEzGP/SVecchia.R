######   scaled Vecchia approximation (Katzfuss, Guinness, Lawrence)  #######


### necessary packages

# install.packages('GpGp')  # need version >= 0.2.2
library(GpGp) 

# install.packages('GPvecchia')
library(GPvecchia)

#####   fitting function   ########

fit_scaled_HEzGP2=function(y,inputs,inputs_encode,ms=c(30),trend='intercept',X,scale='none',q_e,p_e,m_e,
                    var.ini,ranges.ini,select=Inf,print.level=0,max.it=32,tol.dec=4,
                    n.est=min(5e3,nrow(inputs)),find.vcf=FALSE,vcf.scorefun=ls) {
  
  ## dimensions
  n=nrow(inputs)
  d=ncol(inputs)
  
  ## specify trend covariates
  if(missing(X)) {
    if(trend=='zero'){
      X=as.matrix(sample(c(-1,1),n,replace=TRUE))
    } else if(trend=='intercept'){
      X=as.matrix(rep(1,n))
      
    } else if(trend=='linear'){
      X=cbind(rep(1,n),inputs)
    } else if(trend=='pre'){
      X=as.matrix(sample(c(-1,1),n,replace=TRUE))
      beta=mean(y)
      y=y-beta
    } else stop('invalid trend option specified')
  } else trend='X'
  
  ## default variance parameter
  if(missing(var.ini)) {
    #cur.var=summary(stats::lm(y~X-1))$sigma^2
    cur.var=rep(summary(stats::lm(y~X-1))$sigma^2,(q_e+1))
    } else cur.var=var.ini
  
  ## default range parameters
  input.ranges=apply(matrix(inputs[,1:p_e],ncol=p_e),2,function(x) diff(range(x)))
  if(missing(ranges.ini)) cur.ranges=rep(2*input.ranges,2) else cur.ranges=ranges.ini
  active=rep(TRUE,d)
  
  # ## fixed nugget?
  # if(is.null(nug)){
  #   fix.nug=FALSE; nug=.01*var(y)
  # } else fix.nug=TRUE
  # 
  # ## smoothness: fixed? bessel?
  # if(is.null(nu)){
  #   covfun='matern_scaledim'
  #   cur.oth=c(3.5,nug)
  #   fix.nu=FALSE
  # } else if(nu %in% (.5+(1:4))){
  #   covfun=paste0("matern",nu*10,"_scaledim")
  #   cur.oth=nug
  #   fix.nu=FALSE
  # } else {
  #   covfun='matern_scaledim'
  #   cur.oth=c(nu,nug)
  #   fix.nu=TRUE    
  # }
  
  cur.ranges_full=c()
  for(i in 1:p_e){
    cur.ranges.p=cur.ranges[1:p_e]
    cur.ranges.q=cur.ranges[(p_e+1):(2*p_e)]
    cur.ranges_full[i]=sum(cur.ranges.p[i],cur.ranges.q[i])
  }
  cur.ranges_max=c()
  for(i in 1:p_e){
    cur.ranges.p=cur.ranges[1:p_e]
    cur.ranges.q=cur.ranges[(p_e+1):(2*p_e)]
    cur.ranges_max[i]=max(cur.ranges.p[i],cur.ranges.q[i])
  }
  
  

  ## only use subsample for estimation?
  if(n.est<n){
    ind.est=sample(1:n,n.est)
    y.full=y; inputs.full=inputs; X.full=X;inputs.full_encode=inputs_encode;
    y=y[ind.est]; inputs=inputs[ind.est,,drop=FALSE]; X=X[ind.est,,drop=FALSE];inputs_encode=inputs_encode[ind.est,,drop=FALSE];
  }


  # ## only use subsample for estimation?
  # if(n.est<n){
  #   set.seed(1)
  #   #ind.est=sample(1:n,n.est)
  #   r_vecc=round(nrow(train_QQ)/twin_size_Vecc)
  #   ind.est=twin(inputs_encode,r=r_vecc)
  #   y.full=y; inputs.full=inputs; X.full=X;inputs.full_encode=inputs_encode;
  #   y=y[ind.est]; inputs=inputs[ind.est,,drop=FALSE]; X=X[ind.est,,drop=FALSE];inputs_encode=inputs_encode[ind.est,,drop=FALSE];
  # }
  # # 
  ## decrease or remove m values larger than n
  ms=unique(ifelse(ms<length(y),ms,length(y)-1))
  
  
  ### for increasing m
  for(i.m in 1:length(ms)){
    
    m=ms[i.m]
    if(i.m<length(ms)){ tol=10^(-tol.dec-2) } else {tol=10^(-tol.dec)}
    
    ### increase maxit until convergence
    conv=FALSE
    maxit=2
    while(conv==FALSE & maxit<=max.it){
      
      if(print.level>0) {
        print(paste0('m=',m,', maxit=',maxit)); print(cur.ranges)}
      
      ## check for inactive input dims (large range params)
      active=(cur.ranges<input.ranges*select)
      if(sum(active,na.rm=TRUE)==0) stop('all inputs inactive. increase select?')
      cur.ranges[!active]=Inf
      
      
     
      ## specify how to scale input dimensions
      cur.ranges[!active]=Inf
      
      ## order and condition based on current params
      if(scale=='parms'){ scales=sqrt(cur.ranges[1:p_e])
      } else if(scale=='ranges'){ scales=sqrt(input.ranges)
      } else if(scale=='none'){ scales=1
      } else if(scale=='HE'){ scales=sqrt(cur.ranges[(p_e+1):(2*p_e)])
      }else if(scale=='add'){ scales=sqrt(cur.ranges_full)
      }else if(scale=='max'){ scales=sqrt(cur.ranges_max)
      }else if(scale=='both'){  scale1=sqrt(cur.ranges[1:p_e]); scale2=sqrt(cur.ranges[(p_e+1):(2*p_e)]);
      }else stop(paste0('invalid argument scale=',scale))
      
      if(scale=="both"){
        #fit
        theta_0_k=t(inputs_encode[,1:p_e])*scale1
        theta_k=t(inputs_encode[,1:p_e])*scale2
        scaled_locs=theta_0_k+theta_k
        #scaled_locs<-apply(scaled_locs,2,function(x) (x-min(x))/max(x-min(x)) )
        
        inputs.scaled=cbind(t(scaled_locs),inputs_encode[,(p_e+1):(ncol(train_QQ_encode)-1)])
        ord=GPvecchia::order_maxmin_exact(inputs.scaled)
        inputs.ord_encode=inputs_encode[ord,,drop=FALSE]
        inputs.ord=inputs[ord,,drop=FALSE]
        y.ord=y[ord]
        X.ord=X[ord,,drop=FALSE]
        
        theta_0_k_2=t(inputs.ord_encode[,1:p_e])*scale1
        theta_k_2=t(inputs.ord_encode[,1:p_e])*scale2
        scaled_locs_2=theta_0_k_2+theta_k_2
        #scaled_locs_2<-apply(scaled_locs_2,2,function(x) (x-min(x))/max(x-min(x)) )
        inputs.ord_nn=cbind(t(scaled_locs_2),inputs.ord_encode[,(p_e+1):(ncol(train_QQ_encode)-1)])
        NNarray=GpGp::find_ordered_nn((inputs.ord_nn),m)
        
      }else{
      ## order and condition based on current params
     # ord=GPvecchia::order_maxmin_exact(t(t(inputs)*scales))
      scaled_locs=t(inputs_encode[,1:p_e])*scales
      inputs.scaled=cbind(t(scaled_locs),inputs_encode[,(p_e+1):(ncol(train_QQ_encode)-1)])
      ord=GPvecchia::order_maxmin_exact(inputs.scaled)
      inputs.ord_encode=inputs_encode[ord,,drop=FALSE]
      inputs.ord=inputs[ord,,drop=FALSE]
      y.ord=y[ord]
      X.ord=X[ord,,drop=FALSE]
      inputs.ord_nn=cbind(t(t(inputs.ord_encode[,1:p_e])*scales),inputs.ord_encode[,(p_e+1):(ncol(train_QQ_encode)-1)])
      NNarray=GpGp::find_ordered_nn((inputs.ord_nn),m)
      }
      
      
      ## starting and fixed parameters
      cur.parms=c(cur.var,cur.ranges)
      
      fixed=NULL      
      # if(fix.nu) fixed=c(fixed,length(cur.parms)-1)
      # if(fix.nug) fixed=c(fixed,length(cur.parms))
     # active=c(active,rep(TRUE,q_e))
      ## fisher scoring
      fit=fit_model_HEzGP2(y.ord,inputs.ord,X.ord,
                          NNarray=NNarray,m_seq=m,convtol=tol,
                          start_parms=cur.parms,max_iter=maxit,
                          covfun_name=covfun,silent=(print.level<2),
                          reorder=FALSE,fixed_parms=fixed,p_e=p_e,q_e=q_e,m_e=m_h)
      cur.var=fit$covparms[1:(q_e+1)]
      cur.ranges=fit$covparms[(q_e+2):(npar)]
      conv=fit$conv
      maxit=maxit*2
      
    }
  } 
  
  
  ### prepare fit object for subsequent prediction
  fit$covparms=c(cur.var,cur.ranges)
  fit$trend=trend
  fit$p_e=p_e
  fit$q_e=q_e
  fit$m_e=m_e
  fit$NNarray=NNarray
  if(n.est<n){
    fit$y=y.full
    fit$locs=inputs.full
    fit$locs_encode=inputs.full_encode
    fit$X=X.full
  } else {
    fit$locs=inputs.ord
    fit$locs_encode=inputs.ord_encode
  }
  if(trend=='zero') {
    fit$X=as.matrix(rep(0,n))
  } else if(trend=='pre') {
    fit$betahat=beta
    fit$y=fit$y+beta
    fit$trend='intercept'
    fit$X=as.matrix(rep(1,n))
  }
  
  ### find variance correction factor, if requested
  if(find.vcf){
   # fit$vcf=fit_vcf(fit,scale=scale,scorefun=vcf.scorefun,p_e,q_e,m_e)
    fit$vcf=fit_vcf(fit,m.pred=40,n.test=min(1e3,round(nrow(fit$locs)/50)),  scale='HE',scorefun =vcf.scorefun,p_e=p_e,q_e=q_e,m_e=m_h )
  } else fit$vcf=1
  
  return(fit)
  
}





#######   prediction   ########

predictions_scaled_HEzGP2 <- function(fit,locs_pred,locs_pred_encode,m=100,nsims=0,
                               predvar=FALSE,X_pred,scale='none'){
  y_obs = fit$y
  locs_obs = fit$locs
  locs_obs_encode=fit$locs_encode
  X_obs = fit$X
  beta = fit$betahat
  covparms = fit$covparms
  covfun_name = fit$covfun_name
  n_obs <- nrow(locs_obs)
  n_pred <- nrow(locs_pred)
  if(is.null(fit$vcf)) vcf=1 else vcf=fit$vcf
  p_e=fit$p_e
  q_e=fit$q_e
  m_e=fit$m_e
  npar=1 + q_e + p_e + p_e
  
  # ## add nugget for numerical stability
  # if(covparms[length(covparms)]==0) 
  #   covparms[length(covparms)]=covparms[1]*1e-12
  if(scale=="add"){
    par_full=c()
    par.p=covparms[(q_e+1)+(1:(p_e))]
    par.q=covparms[(q_e+2+p_e):npar]
    for(i in 1:p_e){
      par_full[i]=sum(par.p[i],par.q[i])
    }
  }
  
  if(scale=="max"){
    par_max=c()
    par.p=covparms[(q_e+1)+(1:(p_e))]
    par.q=covparms[(q_e+2+p_e):npar]
    for(i in 1:p_e){
      par_max[i]=max(par.p[i],par.q[i])
    }
  }
  # specify trend if missing
  if(missing(X_pred)){
    if(fit$trend=='zero'){
      X_pred=as.matrix(rep(0,n_pred))
    } else if(fit$trend=='intercept'){
      X_pred=as.matrix(rep(1,n_pred))
    } else if(fit$trend=='linear'){
      X_pred=cbind(rep(1,n_pred),locs_pred)
    } else stop('X_pred must be specified')
  }
  
  # specify how to scale input dimensions
  if(scale=='parms'){ scales=sqrt(covparms[(q_e+1)+(1:(p_e))]) 
  } else if(scale=='add'){scales=sqrt(par_full)#scales=1/par_full
  }else if(scale=='HE'){scales=sqrt(covparms[(q_e+2+p_e):npar])
  }else if(scale=='max'){ scales=sqrt(par_max)#
  }else if(scale=='none'){ scales=1
  }else if(scale=='ranges'){ scales=sqrt(apply(locs_obs[,1:p_e],2,function(x) diff(range(x))))
  } else if(scale=='both'){ scale1=sqrt(covparms[(q_e+1)+(1:(p_e))]);scale2=sqrt(covparms[(q_e+2+p_e):npar]);
  }else stop(paste0('invalid argument scale=',scale))
  
  if(scale=="both"){
    #obs
    theta_0_k_obs=t(locs_obs_encode[,1:p_e])*scale1
    theta_k_obs=t(locs_obs_encode[,1:p_e])*scale2
    scaled_locs_obs=theta_0_k_obs+theta_k_obs
    #scaled_locs_obs<-apply(scaled_locs_obs,2,function(x) (x-min(x))/max(x-min(x)) )
    
    #pred
    theta_0_k_pred=t(locs_pred_encode[,1:p_e])*scale1
    theta_k_pred=t(locs_pred_encode[,1:p_e])*scale2
    scaled_locs_pred=theta_0_k_pred+theta_k_pred
    #scaled_locs_pred<-apply(scaled_locs_pred,2,function(x) (x-min(x))/max(x-min(x)) )
    
    inputs.locs_obs=cbind(t(scaled_locs_obs),locs_obs_encode[,(p_e+1):(ncol(train_QQ_encode)-1)])
    inputs.locs_pred=cbind(t(scaled_locs_pred),locs_pred_encode[,(p_e+1):(ncol(train_QQ_encode)-1)])
    
  }else{
    inputs.locs_obs=cbind(t(t(locs_obs_encode[,1:p_e])*scales),locs_obs_encode[,(p_e+1):(ncol(train_QQ_encode)-1)])
    inputs.locs_pred=cbind(t(t(locs_pred_encode[,1:p_e])*scales),locs_pred_encode[,(p_e+1):(ncol(train_QQ_encode)-1)])
    
  }
  
    if(nsims>0) stop('cannot produce joint samples when joint=FALSE')
    
    y  = y_obs - X_obs %*% beta
    
    # find the NNs 
    m=min(m,nrow(locs_obs))
    # NNarray=FNN::get.knnx(t(t(locs_obs)*scales),
    #                       t(t(locs_pred)*scales),m)$nn.index
    # 
  
    
   
    NNarray=FNN::get.knnx(inputs.locs_obs,
                          inputs.locs_pred,m)$nn.index
    
    means=vars=numeric(length=n_pred)
    for(i in 1:n_pred){
      
      # NN conditioning sets
      NN=NNarray[i,]
      
      # (co-)variances
      K=get(covfun_name)(covparms,rbind(locs_obs[NN,],locs_pred[i,]),p_e=p_e,q_e=q_e,m_e=m_h)
      tau=sqrt(.Machine$double.eps)
      K=K+diag(tau, nrow(K))
      cl=t(chol(K))
      
      # prediction
      means[i]=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],y[NN])
      vars[i]=cl[m+1,m+1]^2*vcf
      
    }
    means=means+c(X_pred %*% beta)
    
    if(predvar==FALSE){ 
      preds=list(means=means,NN_pred=NNarray) 
    } else {
      preds=list(means=means,vars=vars,NN_pred=NNarray)
    }
    
    return(preds)
  
}




#######   obs.pred maxmin ordering   ########
order_maxmin_pred<-function(locs, locs_pred,refine=FALSE){
  
  ord<-1:nrow(locs) #GPvecchia::order_maxmin_exact(locs)
  ord_pred <-GPvecchia::order_maxmin_exact(locs_pred)
  
  if(refine){
    
    locs_all = rbind(locs, locs_pred)
    
    n <- nrow(locs)
    m <- min( round(sqrt(n)), 200 )
    
    n_pred <- nrow(locs_pred)
    # next is to find 'ord_pred', a maxmin reordering of prediction locations
    NN <- FNN::get.knn( locs_all, k = m )$nn.index
    #NN_pred <- FNN::get.knnx( locs, locs_pred, k = 1 )$nn.dist
    # use ord, then order by NN_pred
    index_in_position <- c( ord, n + ord_pred, rep(NA,n_pred) )
    position_of_index <- order(index_in_position[1:(n+n_pred)])
    
    # move an index to the end if it is a
    # near neighbor of a previous location
    curlen <- n + n_pred
    nmoved <- 0
    for(j in (n+1):(n+2*n_pred) ){
      # nneigh tells us how many neighbors to look at
      # in order to decide whether the current point
      # has a previously ordered neighbor
      nneigh <- round( min(m,1*(n+n_pred)/(j-nmoved+1)) )
      neighbors <- NN[index_in_position[j],1:nneigh]
      if( min( position_of_index[neighbors], na.rm = TRUE ) < j ){
        nmoved <- nmoved+1
        curlen <- curlen + 1
        position_of_index[ index_in_position[j] ] <- curlen
        index_in_position[curlen] <- index_in_position[j]
        index_in_position[j] <- NA
      }
    }
    
    ord_pred <- index_in_position[ !is.na( index_in_position ) ][(n+1):(n+n_pred)] - n
    
  }
  
  return(list(ord=ord, ord_pred=ord_pred))
}




#######   find NN for prediction locations   ########
find_ordered_nn_pred <- function(locs,m,fix.first=0,searchmult=2){
  
  # if locs is a vector, convert to matrix
  if( is.null(ncol(locs)) ){
    locs <- as.matrix(locs)
  }
  
  # number of locations
  n <- nrow(locs)
  m <- min(m,n-1)
  mult <- 2
  
  # FNN::get.knnx has strange behavior for exact matches
  # so add a small amount of noise to each location
  ee <- min(apply( locs, 2, stats::sd ))
  locs <- locs + matrix( ee*1e-6*stats::rnorm(n*ncol(locs)), n, ncol(locs) )    
  
  # to store the nearest neighbor indices
  NNarray <- matrix(NA,n,m+1)
  
  # to the first mult*m+1 by brute force
  maxval <- min( mult*m + 1, n )
  if(fix.first<=maxval){
    NNarray[1:maxval,] <- GpGp::find_ordered_nn_brute(locs[1:maxval,,drop=FALSE],m)
  } else {
    maxval=fix.first
    NNarray[1:(m+1),] <- GpGp::find_ordered_nn_brute(locs[1:(m+1),,drop=FALSE],m)
    NNarray[1:maxval,1]=1:maxval
    NNarray[(m+1):maxval,1+(1:m)]=matrix(rep(1:m,maxval-m),byrow=TRUE,ncol=m)
  }
  query_inds <- min( maxval+1, n):n
  data_inds <- 1:n
  
  msearch <- m
  
  while( length(query_inds) > 0 ){
    msearch <- min( max(query_inds), round(searchmult*msearch) )
    data_inds <- 1:min( max(query_inds), n )
    NN <- FNN::get.knnx( locs[data_inds,,drop=FALSE], locs[query_inds,,drop=FALSE], msearch )$nn.index
    less_than_k <- t(sapply( 1:nrow(NN), function(k) NN[k,] <= query_inds[k]  ))
    sum_less_than_k <- apply(less_than_k,1,sum)
    ind_less_than_k <- which(sum_less_than_k >= m+1)
    
    NN_m <- t(sapply(ind_less_than_k,function(k) NN[k,][less_than_k[k,]][1:(m+1)] ))
    
    NNarray[ query_inds[ind_less_than_k], ] <- NN_m
    
    query_inds <- query_inds[-ind_less_than_k]
    
  }
  
  return(NNarray)
}



##########   line search for variance correction factor   ###
fit_vcf=function(fit,m.pred=40,n.test=min(1e3,round(nrow(fit$locs)/50)),
                 scale='none',scorefun=ls,p_e,q_e,m_e){
  
  # remove test data from fit object
  fitsearch=fit
  inds.test=sample(1:nrow(fit$locs),n.test)
  fitsearch$y=fit$y[-inds.test]
  fitsearch$locs=fit$locs[-inds.test,,drop=FALSE]
  fitsearch$locs_encode=fit$locs_encode[-inds.test,,drop=FALSE]
  fitsearch$X=fit$X[-inds.test,,drop=FALSE]
  p_e=fitsearch$p_e
  q_e=fitsearch$q_e
  m_e=fitsearch$m_e
  
  # make predictions
  preds=predictions_scaled_HEzGP2(fitsearch,locs_pred=fit$locs[inds.test,,drop=FALSE],
                         locs_pred_encode=fit$locs_encode[inds.test,,drop=FALSE],
                           m=m.pred,predvar=TRUE,scale=scale,
                           X_pred=fit$X[inds.test,,drop=FALSE])
  
  # optimize correction factor
  y.test=fit$y[inds.test]
  objfun=function(vcf) scorefun(y.test,preds$means,preds$vars*vcf)
  vcf=optimize(objfun,c(1e-6,1e6))$minimum
  
  return(vcf)
  
}


### log score
ls=function(dat,mu,sig2) -mean(dnorm(dat,mu,sqrt(sig2),log=TRUE))
