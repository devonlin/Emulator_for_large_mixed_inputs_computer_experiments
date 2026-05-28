#emulation/borehole
rm(list=ls())
#------------------------------------------------------------------------------------
#package needed
library(lhs)
library(mltools)
library(data.table)
#-----------------------------------------------------------------------------------
#dimensions
p = 4
q = 3
m=c(3,3,3)
train_size=200
test_size=500
tau=sqrt(.Machine$double.eps)
nsim=30
zvector1=1:m[1]
zvector2=1:m[2]
zvector3=1:m[3]
z_true=as.matrix(expand.grid(zvector1,zvector2,zvector3))

#sizes
twin1_size= round(min(50 * (p+q), max(sqrt(train_size*prod(m)), 10 * (p+q))))
twin2_size=train_size-twin1_size
r=ceiling(train_size*prod(m)/twin1_size)
local_size=max(25,3*(p+q))
ns = q
nr=test_size
la_size=10
nn_size=local_size+la_size


#number of the cores
no_cores=20
#simulations
nsim_start=1
nsim_end=30

#------------------------------------------------------------------------------------
#true function
#-----------------------------------------------------------------------------
computer_simulator=function(X)
{
  
  y=c()
  for(ic in 1:nrow(X)){
  X[ic,5:7]=(X[ic,5:7]-min(z_true))/(max(z_true)-min(z_true))    
  Tu <- X[ic,1] * (115600 - 63070) + 63070#x_1
  rw <- X[ic,2] * (.15 - .05) + .05 #x_2
  L  <- X[ic,3] * (1680 - 1120) + 1120 #x_3
  Tl <- X[ic,4] * (116 - 63.1) + 63.1#x_4
  Hu <- X[ic,5] * (390 - 290) + 290#x_5
  r  <- X[ic,6] * (50000 - 100) + 100 #x_6
  Kw <- X[ic,7] * (12045 - 9855) + 9855 #x_7
  
  
  frac1 <- 2 * pi * Tu * (Hu)
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  
  y[ic] <- frac1 / frac2
  }
  return((y))
  
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

