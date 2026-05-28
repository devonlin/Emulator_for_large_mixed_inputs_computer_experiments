
#-----------------------------------------------------------------------------
#package needed
library(GpGp)
library(GPvecchia)
library(Rcpp)
library(laGP)
library(mltools)
library(data.table)
#sources
sourceCpp("EzGP_Vecchia.cpp")
source("vecchia_scaled_EzGP.R")
source("fisher_scoring.R")
source("fit_model_EzGP.R")



model_name="Ez-Vecchia_encode-1"
print(model_name)

#-----------------------------------------------------------------------------
#true function
#test simple example
computer_simulator<- function(X)
{
  y=c()
  for(ic in 1:nrow(X)){
    
    X[ic,6:8]=(X[ic,6:8]-min(z_true))/(max(z_true)-min(z_true))    
    theta <- X[ic,5:8]*(2*pi)
    L     <- X[ic,1:4]
    
    thetamat <- matrix(rep(theta,times=4), 4, 4, byrow=TRUE)
    thetamatlow <- thetamat
    thetamatlow[upper.tri(thetamatlow)] <- 0
    sumtheta <- rowSums(thetamatlow)
    
    u <- sum(L*cos(sumtheta))
    v <- sum(L*sin(sumtheta))
    
    y[ic] <- 10*(u^2 + v^2)^(0.5)
  }
  return(y)
}

nsim=30
covfun=paste0("EzGP_cov")
p_e=5
q_e=3
m_h=c(3,3,3)
npar=1+p_e+q_e+p_e*sum(m_h)

train_size=200
test_size=100
#dimensions

zvector1=1:m_h[1]
zvector2=1:m_h[2]
zvector3=1:m_h[3]
z_true=as.matrix(expand.grid(zvector1,zvector2,zvector3))


ms=c(5)
local_size=max(25,3*(p_e+q_e))
la_size=10
nn_size=local_size+la_size
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

RMSE=c()
time_whole=c()
Yhat_each=matrix(NA,nrow=test_size,ncol=nsim)
for(l in 1:nsim){
  print(l)
  #-----------------------------------------------------------------------------
  filename= paste0("data/train_QQ","_",l,".csv",sep="")
  train_QQ=as.matrix(read.csv(filename))[,2:(p_e+q_e+2)]
  filename= paste0("data/test_QQ","_",l,".csv",sep="")
  test_QQ=as.matrix(read.csv(filename))[,2:(p_e+q_e+2)]
  #------------------------------------------------------------------------------------
  train_QQ_1=data.frame(train_QQ)
  test_QQ_1=data.frame(test_QQ)
  #
  if(q_e==1){
    names <- (p_e+1):(p_e+q_e)
    train_QQ_1[,names] <- sapply(train_QQ_1[,names], factor)
    test_QQ_1[,names] <- sapply(test_QQ_1[,names], factor)
  }else{
    
    names <- (p_e+1):(p_e+q_e)
    train_QQ_1[,names] <- lapply(train_QQ_1[,names], factor)
    test_QQ_1[,names] <- lapply(test_QQ_1[,names], factor)
  }
  
  #
  train_QQ_encode=one_hot(as.data.table(train_QQ_1))
  train_QQ_encode=sapply(train_QQ_encode, as.numeric)
  test_QQ_encode=one_hot(as.data.table(test_QQ_1))
  test_QQ_encode=sapply(test_QQ_encode, as.numeric)
  #------------------------------------------------------------------------------------
  #test location 
  
  start_time_whole <- Sys.time()
  fit.sv=fit_scaled_EzGP(y=train_QQ[,1+p_e+q_e],inputs=train_QQ[,1:(p_e+q_e)],
                         inputs_encode=train_QQ_encode[,1:(ncol(train_QQ_encode)-1)],
                         ms=ms[1],n.est=min(3e3,train_size*prod(m_h)),p_e=p_e,q_e=q_e,m_e=m_h,scale="none")
  Yhat_each[,l]=predictions_scaled_EzGP(fit=fit.sv,locs_pred=test_QQ[,1:(p_e+q_e)],locs_pred_encode=test_QQ_encode[,1:(ncol(test_QQ_encode)-1)],m=nn_size,scale = "none")$preds
  
  end_time_whole <- Sys.time()
  time_whole[l]=as.numeric(end_time_whole - start_time_whole,unit = "mins")
  RMSE[l]=sqrt(mean((Yhat_each[,l]-test_QQ[,1+p_e+q_e])^2))
}

filename= paste0("RMSE","_",model_name,"_",".rds",sep="")
saveRDS(RMSE,file=filename)
filename= paste0("Yhat","_",model_name,"_",".rds",sep="")
saveRDS(Yhat_each,file=filename)
filename= paste0("time","_",model_name,"_",".rds",sep="")
saveRDS(time_whole,file=filename)

