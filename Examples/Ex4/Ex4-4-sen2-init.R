rm(list=ls())
#------------------------------------------------------------------------------------
#package needed
library(lhs)
library(mltools)
library(data.table)
#-----------------------------------------------------------------------------------
#dimensions
p = 1
q = 6
m=c(3,3,3,3,3,3)
train_size=5
test_size=10
tau=sqrt(.Machine$double.eps)
nsim=30
zvector1=1:m[1]
z_true=as.matrix(expand.grid(zvector1,zvector1,zvector1,zvector1,zvector1,zvector1))
z_true_1=data.frame(z_true)
#
if(q==1){
  z_true_1 <- sapply(z_true, factor)
}else{
  
  z_true_1 <- lapply(z_true_1, factor)
}

#
z_true_encode=one_hot(as.data.table(z_true_1))
z_true_encode=sapply(z_true_encode, as.numeric)
#sizes
twin1_size= round(min(50 * (p+q), max(sqrt(train_size*prod(m)), 10 * (p+q))))
twin2_size=train_size*prod(m)-twin1_size
r=round(train_size*prod(m)/twin1_size)
local_size=max(25, 3 * (p+q))
ns_twin = c(q-2)
ns =c(q-1)
la_size=10
nn_size=local_size+la_size
nr=test_size*10
#number of cores
no_cores=20
#simulations
nsim_start=21
nsim_end=30

#------------------------------------------------------------------------------------
#true function
#-----------------------------------------------------------------------------
computer_simulator=function(xx)#Example function
{
  y=c()
  for(l in 1:nrow(xx)){
    z=xx[l,2:7]
    x=xx[l,1]
    
    if(z[1]==1)
      a=cos(3*pi*x)
    else if (z[1]==2)
      a=cos(4*pi*x)
    else if (z[1]==3)
      a=cos(5*pi*x)
    
    if(z[2]==1)
      b=sin(3*pi*x)
    else if (z[2]==2)
     b=sin(4*pi*x)
    else if (z[2]==3)
      b=sin(5*pi*x)
    
    
    
    if(z[3]==1)
      c=1
    else if (z[3]==2)
      c=2
    else if (z[3]==3)
      c=3
    
    if(z[4]==1)
      d=0.1
    else if (z[4]==2)
      d=0.2
    else if (z[4]==3)
      d=0.3
    
    
 
    if(z[5]==1)
      e=x
    else if (z[5]==2)
      e=x^2
    else if (z[5]==3)
     e=x^3
    
    if(z[6]==1)
      h=-2
    else if (z[6]==2)
      h=0
    else if (z[6]==3)
      h=-2
    
    y[l] <- h+e+((a+b)*c*d)
    
  }
  return(y)
  
}


save.image("initial.RData")

#------------------------------------------------------------------------------------
#traning data
for( j in 1:nsim){
  x_M_c0=matrix(NA,nrow=train_size*prod(m),ncol=p)
  for(i in 1:prod(m)){
    x_M_c0[(((i-1)*(train_size))+1):(i*train_size),] = randomLHS(train_size,p)
    
  }
  z_QQ=apply(z_true,2,rep,each=train_size)
  train_QQ_XZ=cbind(x_M_c0,z_QQ)
  train_QQ_XZ = as.matrix(train_QQ_XZ)
  train_QQ_Y=computer_simulator(train_QQ_XZ)
  train_QQ=cbind(train_QQ_XZ,train_QQ_Y)
  filename= paste0("train_QQ","_",j,".csv",sep="")
  write.csv(train_QQ,file = filename)
  #test data
  #different level
  x_M_c0=matrix(NA,nrow=test_size*prod(m),ncol=p)
  for(i in 1:prod(m)){
    x_M_c0[(((i-1)*(test_size))+1):(i*test_size),] = randomLHS(test_size,p)
    
  }
  z_QQ=apply(z_true,2,rep,each=test_size)
  test_QQ_XZ=cbind(x_M_c0,z_QQ)
  test_QQ_XZ = as.matrix(test_QQ_XZ)
  test_QQ_Y=computer_simulator(test_QQ_XZ)
  test_QQ_full=cbind(test_QQ_XZ,test_QQ_Y,1:(prod(m)*test_size))
  id_test=sample(1:nrow(test_QQ_full),test_size*10)
  test_QQ=test_QQ_full[sort(id_test),]
  filename= paste0("test_QQ","_",j,".csv",sep="")
  write.csv(test_QQ,file = filename)
}

#-----------------------------------------------

