rm(list=ls())
#------------------------------------------------------------------------------------
#package needed
library(lhs)
library(mltools)
library(data.table)
#-----------------------------------------------------------------------------------
#dimensions
p = 3
q = 3
m=c(3,3,3)
train_size=200
test_size=100
tau=sqrt(.Machine$double.eps)
nsim=30
zvector1=1:m[1]
zvector2=1:m[2]
zvector3=1:m[3]
z_true=as.matrix(expand.grid(zvector1,zvector2,zvector3))


twin1_size= round(min(50 * (p+q), max(sqrt(train_size*prod(m)), 10 * (p+q))))
twin2_size=train_size-twin1_size
r=ceiling(train_size*prod(m)/twin1_size)
local_size=max(25,3*(p+q))
ns = q
ns_twin=q
nr=test_size
la_size=10
nn_size=local_size+la_size
#number of cpu
no_cores=20
#simulation
nsim_start=1
nsim_end=30


#-----------------------------------------------------------------------------
computer_simulator=function(xx)#Example function
{
  y=c()
  for(l in 1:nrow(xx)){
    z=xx[l,4:6]
    x=xx[l,1:3]
    
    if(z[1]==1)
      f=x[1]+x[2]^2+x[3]
    else if (z[1]==2)
      f=x[1]^2+x[2]+x[3]
    else if (z[1]==3)
      f=x[1]+x[2]^2+x[3]
    
    if(z[2]==1)
      g=cos(x[1])+cos(x[2]*2)+cos(x[3])
    else if (z[2]==2)
      g=cos(x[1])+cos(x[2]*2)+cos(x[3])
    else if (z[2]==3)
      g=cos(x[1]*2)+cos(x[2])+cos(x[3])
    
    if(z[3]==1)
      h=sin(x[1])+sin(x[2]*2)+sin(x[3])
    else if (z[3]==2)
      h=sin(x[1])+sin(x[2]*2)+sin(x[3])
    else if (z[3]==3)
      h=sin(x[1]*2)+sin(x[2])+sin(x[3])
    
    y[l] <- 100*(f+(g+h))
    
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
  id_test=sample(1:nrow(test_QQ_full),test_size)
  test_QQ=test_QQ_full[sort(id_test),]
  filename= paste0("test_QQ","_",j,".csv",sep="")
  write.csv(test_QQ,file = filename)
}

#-----------------------------------------------

