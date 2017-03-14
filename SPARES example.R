rm(list = ls())
source("SPARES.R")

n=150
p=500
p0=5

#########beta^0 and the active set
s0<-sort(sample(1:p,p0))
beta0<-rep(0,p)
beta0[s0]<-runif(p0,0.5,2)
round(cbind(s0,beta0[s0]),2)

################AR(1) correlation structure
rho=0.8
H <- abs(outer(1:p, 1:p, "-")) 
sigma <- rho^H 
dim(sigma)
round(sigma[1:10,1:10],2)

######generate data
xmat<-mvrnorm(n=n,mu=rep(0,p),Sigma=sigma)
err<-rnorm(n,1,1)
yvec<-xmat%*%beta0+err

##########oracle estimates
lm0<-lm(yvec~xmat[,s0])
summary(lm0)

###########run SPARES
sps<-SPARES(xmat,yvec,selection="lasso",nloop=500)
summary(sps)
#sm.beta: SPARES estimator
#sd: standard deviation estimate
#p: p-values of testing \beta^0_j = 0
#sel.freq: selection frequency for the predictors
#

sps$int
sps$nzero
plot(sps$sm.beta,ylim=c(-0.5,2.5))
points(s0,beta0[s0],pch=19,col="red")
plot(as.table(sps$sel.freq))

s0

coefmat<-cbind(1:p,sps$sm.beta,sps$p)
head(coefmat[order(sps$p),])
