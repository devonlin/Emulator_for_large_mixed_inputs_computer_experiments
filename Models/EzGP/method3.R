#-------------------------------------------------------------------------------
#Packages Needed
#-------------------------------------------------------------------------------

library(EzGP)
library(mltools)
library(data.table)
library(doParallel)
library(laGP)
#-------------------------------------------------------------------------------
#For fast update of EzGP and ALC
#-------------------------------------------------------------------------------

source("EzGP_Functions.R")
#-------------------------------------------------------------------------------
#initials 
#-------------------------------------------------------------------------------
model_name="La-alc-EzGP-node"
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
  npar = 1 + q + p + p*sum(m)
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
#laGP adaptive design
  close=min((1000+nn_size), nrow(train_QQ))# define number of the design points points to consider around the test location
 
   parafunc_foreach=function(iu,p,q,m,tau,local_size,la_size){
     #test location
     w_star=matrix(test_QQ_encode[iu,1:(ncol(test_QQ_encode)-1)],nrow = 1)
     #NN 
     dst <- drop(distance(w_star,train_QQ_encode[,1:(ncol(train_QQ_encode)-1)]))
      if(is.matrix(dst)) dst <- apply(dst, 2, min)
     #initials 
     La_train=train_QQ[order(dst)[1:local_size],]
     # design points to select nest points from
      test_la <- train_QQ[order(dst)[(local_size+1):close],]
      #model initial 
      model = EzGP_fit_partioned((La_train[,1:(p+q)]),La_train[,(p+q+1)],p=p,q=q,m=m,tau=tau)
      #select next points from test_la
      for(i in 1:la_size){
        ALC_total=c()
        for(al in 1:nrow(test_la)){
          #ALC maximization
          ALC_value=update_EzGP_alc(xm1 = test_la[al, 1:(p + q)], Xn = matrix(test_QQ[iu, 1:(p + q)], nrow = 1), X=model$data$X,Y= model$data$Y,
                               n=nrow(model$data$X),parv=model$param,covm= model$cov$Phi, invc= model$cov$Phi_inv,m=model$data$m)
          ALC_total=c(ALC_total,ALC_value)
        }
        ALC=which.min(-(ALC_total))
        #add next point
        new_input=matrix(test_la[ALC,],nrow=1)
        La_train=rbind(La_train,new_input)
        #update the adaptive deisgn using the same paramters
        model <- update_EzGP_adaptive(xm1=new_input,La_train=La_train,X=model$data$X, Y=model$data$Y, n=nrow(model$data$X),
                                 parv=model$param,covm= model$cov$Phi, invc= model$cov$Phi_inv,m=model$data$m)
     
        test_la=test_la[-ALC,]
      }
      #update the paramters after all added points selected
      model = EzGP_fit_partioned((La_train[,1:(p+q)]),La_train[,(p+q+1)],p=p,q=q,m=m,tau=tau)
      pred= EzGP_predict(test_QQ[iu,1:(p+ q)], model, MSE_on = 0)
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
                     .packages=c('nloptr',"EzGP","laGP")) %dopar% {
                       parafunc_foreach(iu=iu,p=p,q=q,m=m,tau=tau,local_size=local_size,la_size=la_size)
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

