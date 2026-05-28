#-------------------------------------------------------------------------------
#Packages Needed
#-------------------------------------------------------------------------------

library(twinning)
library(laGP)
library(EzGP)
library(mltools)
library(data.table)
library(doParallel)
#-------------------------------------------------------------------------------
#For fast update of EzGP 
#-------------------------------------------------------------------------------

source("EzGP_Functions.R")
#-------------------------------------------------------------------------------
#initials 
#-------------------------------------------------------------------------------
model_name="Twin-FEzGP-node"
print(model_name)
load("initial.RData")
RMSE=c()
time_whole=c()
Yhat_each=matrix(NA,nrow=nr,ncol=nsim)
#-------------------------------------------------------------------------------
#Simulation replication
#-------------------------------------------------------------------------------
for(l in nsim_start:nsim_end){
print(l)
#-------------------------------------------------------------------------------
#Data 
#-------------------------------------------------------------------------------
filename= paste0("data/train_QQ","_",l,".csv",sep="")
train_QQ=as.matrix(read.csv(filename))[,2:(p+q+2)]
filename= paste0("data/test_QQ","_",l,".csv",sep="")
test_QQ=as.matrix(read.csv(filename))[,2:(p+q+2)]
#-------------------------------------------------------------------------------
#One-hot encoding
#-------------------------------------------------------------------------------
train_QQ_1=data.frame(train_QQ)
test_QQ_1=data.frame(test_QQ)
#
if(q==1){
  names <- (p+1):(p+q)
  train_QQ_1[,names] <- sapply(train_QQ_1[,names], factor)
  test_QQ_1[,names] <- sapply(test_QQ_1[,names], factor)
}else{
  
  names <- (p+1):(p+q)
  train_QQ_1[,names] <- lapply(train_QQ_1[,names], factor)
  test_QQ_1[,names] <- lapply(test_QQ_1[,names], factor)
}

#
train_QQ_encode=one_hot(as.data.table(train_QQ_1))
train_QQ_encode=sapply(train_QQ_encode, as.numeric)
test_QQ_encode=one_hot(as.data.table(test_QQ_1))
test_QQ_encode=sapply(test_QQ_encode, as.numeric)

#-------------------------------------------------------------------------------
#twinning 
#-------------------------------------------------------------------------------
set.seed(1) 
twin1_indices= twin(train_QQ_encode, r=r)
#-------------------------------------------------------------------------------
#Build g global points and l local points
#-------------------------------------------------------------------------------
twin1=train_QQ[twin1_indices,]
twin2=train_QQ[-twin1_indices,]

#global
start_time <- Sys.time()
model_gg = EzGP_fit_partioned(twin1[,1:(p+q)],twin1[,(p+q+1)],p=p,q=q,m=m,tau=tau) 
end_time <- Sys.time()
time_global=as.numeric(end_time - start_time,unit = "mins")


#------------------------------------------------------------------------------------
#local points change using parallel
parafunc_foreach=function(iu,twin1,p,q,m,tau,local_size,twin1_indices){
  #test location
  w_star=matrix(test_QQ_encode[iu,1:(ncol(test_QQ_encode)-1)],nrow = 1)
  #NN
  dst <- drop(distance(w_star,train_QQ_encode[-twin1_indices,1:(ncol(train_QQ_encode)-1)]))
  if(is.matrix(dst)) dst <- apply(dst, 2, min)
  twin_kd_tree=twin2[order(dst)[1:local_size],]# select l local points using NN
 #------------------------------------------------------------------------------------
  #update  the effect of l local points
  model_loglike_update=update_EzGP_loglike(data=twin_kd_tree, model=model_gg, tau = tau)
  pred= EzGP_predict(test_QQ[iu,1:(p+q)], model_loglike_update, MSE_on = 0)
  Yhat= pred$Y_hat#save
  return(list(Yhat=Yhat))   
}
#parallel process  
# Create an array from the NODESLIST environment variable
nodeslist = unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))
# Create the cluster with the nodes name. One process per count of node name.
# nodeslist = node1 node1 node2 node2, means we are starting 2 processes on node1, likewise on node2.
cl = makeCluster(nodeslist, type = "PSOCK") 
registerDoParallel(cl)
start_time_whole <- Sys.time()
results <- foreach(iu=icount(nr), .combine='rbind',
                   .packages=c('nloptr', "laGP","EzGP","twinning")) %dopar% {
                     parafunc_foreach(iu=iu,twin1=twin1,p=p,q=q,m=m,tau=tau,local_size=local_size,twin1_indices=twin1_indices)
                   }
end_time_whole <- Sys.time()
stopCluster(cl) 
time_whole[l]=as.numeric(end_time_whole - start_time_whole,unit = "mins")

#get predictions
Yhat=unlist(results[,1])
Yhat_each[,l]=Yhat

#------------------------------------------------------------------------------------
#RMSE
RMSE[l]=sqrt(mean(((Yhat-test_QQ[1:nr,(p+q+1)])^2)))

}
#------------------------------------------------------------------------------------
#save data
#------------------------------------------------------------------------------------
filename= paste0("RMSE","_",model_name,"_",".rds",sep="")
saveRDS(RMSE,file=filename)
filename= paste0("Yhat","_",model_name,"_",".rds",sep="")
saveRDS(Yhat_each,file=filename)
filename= paste0("time","_",model_name,"_",".rds",sep="")
saveRDS(time_whole,file=filename)

