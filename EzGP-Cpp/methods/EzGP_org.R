#-----------------------------------------------------------------------------
# packages
#-----------------------------------------------------------------------------
library(twinning)
library(dbscan)
library(EzGP)
library(mltools)
library(data.table)

#-----------------------------------------------------------------------------
# Serial Job 
#-----------------------------------------------------------------------------
f=function(iteration){
#-----------------------------------------------------------------------------
# Initials
#-----------------------------------------------------------------------------
model_name="Twin-EzGP_org"
print(model_name)
load("initial.RData")
#-----------------------------------------------------------------------------
# data
#-----------------------------------------------------------------------------
  filename= paste0("data/train_QQ","_",iteration,".csv",sep="")
  train_QQ=as.matrix(read.csv(filename))[,2:(p+q+2)]
  filename= paste0("data/test_QQ","_",iteration,".csv",sep="")
  test_QQ=as.matrix(read.csv(filename))[,2:(p+q+2)]
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
  
  #
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

  Yhat=c()
  time=c()
 
#------------------------------------------------------------------------------------
#Prediction with each test location
#------------------------------------------------------------------------------------
  #test location
     for(iu in 1:nr){ 
       print(iu)
      w_star=test_QQ_encode[iu,1:(ncol(test_QQ_encode))]
      twin22=rbind((train_QQ_encode[-twin1_indices,]),matrix(w_star,nrow=1))
      nn <- kNN(twin22[,1:(ncol(train_QQ_encode)-1)], k = local_size)
      #  # explore neighborhood point local size
      id <- nrow(twin22)
      twin_kd_tree=twin2[nn$id[id,],]
      if(is.matrix(twin_kd_tree)!=TRUE) twin_kd_tree= matrix(twin_kd_tree,1,ncol=length(twin_kd_tree))
      
    #local-global
      local_global_data=rbind(twin1,twin_kd_tree)
#------------------------------------------------------------------------------------
      #add together
      start_time <- Sys.time()
      model= EzGP_fit(local_global_data[,1:(p+q)],local_global_data[,(p+q+1)],p=p,q=q,m=m,tau=tau)
      end_time <- Sys.time()
      pred= EzGP_predict(test_QQ[iu,1:(p+ q)], model, MSE_on = 0)
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

}
iteration=as.integer(commandArgs(trailingOnly = TRUE))
f(iteration)



  
  
  
  
  
