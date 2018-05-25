#####################R main emulator codes
library(MASS)
source("abind.R",echo=FALSE)
#library(abind)
##predict nutrient concentration overtime
coef=as.matrix(read.table('coef.txt'))
new=as.matrix(read.table('input/input.txt',header=FALSE))
namm0=c("nh3","no2","no3","o2","co2","Biomass","nparticle","Height","Roughness","time")
ll=length(namm0)
time=seq(500,new[,ll],500)
ltime=length(time)
new2=matrix(NA,nrow=ltime,ncol=ll)
for(i in 1:ltime){
	new2[i,]=new
}
new2[,ll]=time
nutrient=abs(new2%*%coef)
input=new[,-ll,drop=FALSE]
input=log(input)
nutrient=log(nutrient)
###compute correlation
mycor=function(X1,X2,betai){
pdf<- diag(betai, nrow = length(betai))
tx1<- t(X1)
tx2 <- t(X2)
R1 <- X1%*%pdf%*%tx1 
R2 <- X2%*%pdf%*%tx2
S1 <- t(as.matrix(diag(R1)))%x%rep(1,nrow(X2))
S2 <- as.matrix(diag(R2))%x%t(rep(1,nrow(X1)))
a1 <- t(tx1)%*%pdf%*%tx2
a2 <- t(tx2)%*%pdf%*%tx1
return(exp(Re(t(a1)+a2-S1-S2)))
}
############
nout=4
ntime=ltime
nsim=2##number of Monte carlo simulation
X=as.matrix(read.table('X.txt',header=TRUE))##training data
Y=as.matrix(read.table('Y.txt',header=TRUE))##tr
namm=colnames(X)
tau=2
#scale2=c(0.00001000,0.00001000,0.00001000,0.00619307,0.00001000,2.00000000,2.00000000,2.00000000,2.00000000,0.00001000)
scale2=c( 0.201935506, 0.006313825, 0.000010000, 0.177466408, 0.000010000, 6.577416780, 7.721354529, 6.594165598,4.081750274, 0.004609465)
#scale2=c(1.000000e-05,3.299517e+01, 1.000000e-05 ,5.913863e-04, 1.000000e-05, 4.371568e+01, 2.431108e+01, 5.567225e+01 ,6.575740e+01, 1.000000e-05)
id2=length(scale2)
theta=scale2[-id2]	
load("Ainv")
#
pred=function(newdata){
X0=newdata#scale(newdata,center=center,scale=scale)
colnames(X0)=namm
form=~nh3+no2+no3+o2+co2+Biomass+nparticle+Height+Roughness
H<- X# model.matrix(form,data=as.data.frame(X))
## Define number of columns of output and model matrices
m<-dim(H)[2]
n<-dim(Y)[1]
n2<- dim(newdata)[1]
######model matrix 
H0<-X0#model.matrix(form,as.data.frame(X0))
A01=mycor(X0,X,betai=theta)###cross correlation
A00=mycor(X0,X0,betai=theta)##test point correlation
A00=A00+(scale2[id2]*diag(tau,dim(A00)))
iOmega=solve(t(H)%*%Ainv%*%H)
betahat=solve(t(H)%*%Ainv%*%H)%*%(t(H)%*%Ainv%*%Y)
mu_star=H0%*%betahat+t(A01)%*%Ainv%*%(Y-H%*%betahat)
r1=H0-(t(A01)%*%Ainv%*%H)
c_star=A00-(t(A01)%*%Ainv%*%A01)+(((r1)%*%(iOmega)%*%t(r1)))
Sigma=(t(Y-H%*%betahat)%*%Ainv%*%(Y-H%*%betahat))/(n-m)
svar=list()
for(i in 1:n2){
svar[[i]]=(diag(c_star)[i]*diag(Sigma))
#svar[[i]]=diag(c_star)%x%Sigma
}
K0=abind(svar,along=2)
out=list(mu=mu_star,K=t(K0))
return(out)
}
##

f=function(newdata){
npoints=1
mstar=list();result1=result2=list();mstar2=list()
m1=pred(newdata=newdata)
mu=m1$mu###predictions
varr=m1$K#variance
mstar=mvrnorm(nsim,c(mu),diag(c(varr),ncol=nout,nrow=nout))
V1=varr
V2=0
mstar2[[1]]=(mstar)
result1[[1]]=(mu)
result2[[1]]=(V1+V2)
print(paste("iteration",1,sep=" "))
##
for(j in 2:ntime){
g2a=mstar2[[j-1]]
inp=matrix(rep(nutrient[j,],nsim),nrow=nsim,byrow=TRUE)
input=cbind(inp,g2a)
m1=pred(newdata=input)
p1=array(m1$mu,c(npoints,nsim,nout))
p2=array(m1$K,c(npoints,nsim,nout))
mu=t(matrix(apply(p1,3,rowMeans)))
mu2=array(rep(mu,nsim),c(npoints,nout,nsim));mu2=aperm(mu2,c(1,3,2))
varr=apply(p2,c(1,3),mean)#expectation of variance
varr2=apply((p1-mu2)^2,3,rowMeans)#variance of expectation
V1=varr;V2=varr2
mstar=mvrnorm(nsim,c(mu),diag(c(V1),ncol=nout,nrow=nout))
mstar2[[j]]=(mstar)
result1[[j]]=(mu)
result2[[j]]=(V1)
print(paste("iteration",j,sep=" "))
}
out1=matrix(unlist(result1),nrow=ntime,ncol=nout,byrow=TRUE)
out2=matrix(unlist(result2),nrow=ntime,ncol=nout,byrow=TRUE)
upper=exp(out1+2*sqrt(out2))
lower=exp(out1-2*sqrt(out2))
VV=((upper-lower)/4)^2
out1=(out1)
out1=rbind(new[,c(6:9)],exp(out1))
drdt=(out1[-1,]-out1[-nrow(out1),])/500
colnames(drdt)=colnames(Y)
drdt=cbind(time,drdt)
colnames(VV)=colnames(out1)=colnames(out2)=colnames(Y)
#write.table(VV,file="output/variance.txt",row.names=FALSE)
write.table(drdt,file="output/Rate.txt",row.names=FALSE,col.names = TRUE)
#write.table(out1,file="output/out1.txt",row.names=FALSE)
#write.table(out2,file="output/output2.txt",row.names=FALSE)
}
tmpf=f(input)
