rm(list = ls())
source("SPARES.R")

n=150
p=500
p0=5

s0<-sort(sample(1:p,p0))
beta0<-rep(0,p)
beta0[s0]<-runif(p0)+0.5
cbind(s0,beta0[s0])

################AR(1)
rho=0.8
H <- abs(outer(1:p, 1:p, "-")) 
sigma <- rho^H 
dim(sigma)
sigma[1:10,1:10]

xmat<-mvrnorm(n=n,mu=rep(0,p),Sigma=sigma)
err<-rnorm(n,1,1)
yvec<-xmat%*%beta0+err

l0<-cv.glmnet(xmat,yvec)
lam0<-l0$lambda.1se
coef0<-coef(l0,s=lam0)[-1]
set0<-(1:p)[coef0!=0]
set0

sp1<-SPARE(scale(xmat),yvec-mean(yvec),selection="lasso",n,p,lam0)
sp1$sel.set
plot(sp1$beta.hat,ylim=c(-0.5,1.5))
points(s0,beta0[s0],pch=19,col="red")

sps<-SPARES(xmat,yvec,selection="lasso",nloop=500)
summary(sps)
sps$int
plot(sps$sm.beta,ylim=c(-0.5,1.5))
points(s0,beta0[s0],pch=19,col="red")
plot(as.table(sps$sel.freq))
(1:p)[sps$p<0.05/p]

coefmat<-cbind(1:p,sps$sm.beta,sps$p)
head(coefmat[order(sps$p),])
s0
