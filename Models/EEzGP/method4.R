#-------------------------------------------------------------------------------
#Packages Needed
#-------------------------------------------------------------------------------

library(mltools)
library(data.table)
library(doParallel)
#-------------------------------------------------------------------------------
#Functions for HEzGP modeling and fast update and prediction 
#-------------------------------------------------------------------------------

source("HEzGP_Functions.R")
source("HEzGP_fit.R")
source("HEzGP_predict.R")
source("cov_Hm.R")
#-------------------------------------------------------------------------------
#initials 
#-------------------------------------------------------------------------------

model_name="LHEzGP-ns-1-relabel"
print(model_name)
load("initial.RData")
#-----------------------------------------------------------------------------

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

#------------------------------------------------------------------------------------
#local points change using parallel
parafunc_foreach=function(iu,p,q,m,tau,ns_twin){
  #qualitative inputs
    tar_z= as.numeric((test_QQ[iu, (p+1):(p+q)]))
    num = nrow(train_QQ)
    rows = c()
    i = 1
    for (j in 1:num){
      if (sum(train_QQ[j,(p+1):(p+q)] == tar_z) >= ns[1]){#extract the matches ns is the number of the match
        rows[i] = j
        i = i + 1
      }
    }
   LEzGP_data=train_QQ[rows,]#extract the data

#------------------------------------------------------------------------------------
#model LEzGP
    model = HEzGP_fit(LEzGP_data[,1:(p+q)],LEzGP_data[,(p+q+1)],p=p,q=q,m=m,tau=tau) 
    pred= HEzGP_predict(test_QQ[iu,1:(p+ q)], model, MSE_on = 0)
    Yhat= pred$Y_hat#save
  return(list(Yhat=Yhat))  
}

#parallel process 
# Create an array from the NODESLIST environnement variable
nodeslist = unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))
# Create the cluster with the nodes name. One process per count of node name.
# nodeslist = node1 node1 node2 node2, means we are starting 2 processes on node1, likewise on node2.
cl = makeCluster(nodeslist, type = "PSOCK") 

registerDoParallel(cl)
start_time_whole <- Sys.time()
results <- foreach(iu=icount(nr), .combine='rbind',
                   .packages=c('nloptr')) %dopar% {
                     parafunc_foreach(iu=iu,p=p,q=q,m=m,tau=tau,ns_twin=ns_twin)
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
  
  