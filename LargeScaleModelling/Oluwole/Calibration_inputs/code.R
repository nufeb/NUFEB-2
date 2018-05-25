###########################################################################processing data
rm(list=ls())
y1=as.matrix(read.table("C:/Octave/mydata/y1.dat",header=FALSE))
y2=as.matrix(read.table("C:/Octave/mydata/y2.dat",header=FALSE))
d1=read.table("C:/Octave/mydata/design1.dat",header=FALSE)
d2=read.table("C:/Octave/mydata/design2.dat",header=FALSE)
library(abind)
library(emulator)
library(MASS)
#normalize data
min1=apply(d1,2,min)
max1=apply(d1,2,max)
min2=apply(d2,2,min)
max2=apply(d2,2,max)
X=matrix(NA,nrow=60,ncol=8)
X2=matrix(NA,nrow=60,ncol=8)
for(i in 1:8){
X[,i]=(d1[,i]-min1[i])/(max1[i]-min1[i])
X2[,i]=(d2[,i]-min1[i])/(max1[i]-min1[i])
}
#output
y1log=log(y1)##CMoutputlog
y2log=log(y2)##CMoutputUntriedTruelog
y=y1log
###########
Ntime=1827
Ptvar=p=2
m0 = matrix(1,nrow=Ptvar,ncol=1)
C0 = 10*diag(1,Ptvar)
n0 = 1
s0 = 10
dis=c(.9,.9)
# iterations: each iteration consists of 200 MH draws.
NIter = 110#0
Nburnin = 10#0
Pinput=8
r=60
Beta=array(1,c(Pinput,NIter))
Alpha=array(1.999,c(Pinput,NIter))
alphai=Alpha[,1]
betai=Beta[,1] 
design=X#d1
#load C:\Octave/share/octave/packages/statistics
################################################################BEGIN
myvar <- function (u,alpha,beta) {
n2=ncol(u)
n1=nrow(u)
 cM = matrix(0,n1,ncol = n1)
  for (i in 1:n1) {	  
for(j in 1:n1){	  
temp=0	  
for(k in 1:n2){
temp = temp + beta[k]*abs(u[i, k]- u[j, k])^alpha[k]	
}	
cM[i,j]=exp(-temp)	
cM[i,j]=cM[j,i]	  
}}
cM
}
#Sigma=myvar(design,alphai,betai)
#Sigma=corr.matrix(design,scales=betai)    
############t_var1
tvarm1=function(y,p,dis,m0,C0,s0,n0,Sigma){
arx=as.matrix(y)
r=ncol(y)
for(i in 1:r){
arx[,i]=y[,i]-rep(colMeans(y)[i],each=T)
}
T=Ntime
bv=dis[1]; bw=dis[2]	
m=list();s=list();n=list();C=list()#;e=list()
length(s)=length(n)=length(m)=length(C)=T
A=R=Q=e=d=list()
m[[1]]=m0
C[[1]]=C0
s[[1]]=s0
n[[1]]=n0
d[[1]]=s0*n0
#d0=d[[1]]
for(t in p:(T-1)){
F=arx[(t):(t-1),]
R[[t]]=C[[t-1]]/bw
Q[[t]]=(t(F)%*%R[[t]]%*%F)+ (s[[t-1]]*Sigma)
A[[t]]=R[[t]]%*%F%*%solve(Q[[t]]) 
e[[t]]=arx[t+1,]-t(F)%*%m[[t-1]]
m[[t]]=m[[t-1]]+A[[t]]%*%e[[t]]
d[[t]]=bv*d[[t-1]]+(s[[t-1]]%*%t(e[[t]])%*%solve(Q[[t]])%*%e[[t]])
n[[t]]=(bv*n[[t-1]])+r
s[[t]]=c(d[[t]]/n[[t]])
C[[t]]=(R[[t]]-A[[t]]%*%Q[[t]]%*%t(A[[t]]))*c(s[[t]]/s[[t-1]])
C[[t]]=(C[[t]]+t(C[[t]]))/2##make it a symmetric matrix
}
C[[1]]=C[[p]]#
m[[1]]=m[[p]]#
s[[1]]=s[[p]]#
n[[1]]=n[[p]]#
#now ad-hoc treatment of 1st p_tvar values 
C[p:(t+1)]=C[(p-1):t];m[p:(t+1)]=m[(p-1):t]
s[p:(t+1)]=s[(p-1):t];n[p:(t+1)]=n[(p-1):t]
CC=abind(C,along=3)
mm=abind(m,along=2)
nn=abind(n,along=1)
ss=abind(s,along=1)
list(CC=CC,nn=nn,ss=ss,mm=mm)
}

#res=tvarm1(y,p,dis,m0,C0,s0,n0,Sigma)
#####variance
myV=function(nn,ss,p,dis){
bw=dis[2]
T=length(nn)
V=rep(0,T)
nt=nn[T]
dt=nn[T]*ss[T]
VT=rgamma(1,shape=nt/2,dt/2)
V[T]=1/VT
for(t in (T-1):(p+1)){
nt=(1-bw)*nn[t]
dt=nn[t]*ss[t]
VT=(bw*VT)+rgamma(1,shape=nt/2,dt/2)	
V[t]=1/VT	
}
V[1:p]=V[p+1]
V=V
}
#res2=myV(res$nn,res$ss,p,dis)
#Vi=res2
################var2
#Vi=myV(res$nn,res$ss,p,dis)
tvarm2=function(y,p,dis,m0,C0,Vi,Sigma){
arx=as.matrix(y)
r=ncol(y)
for(i in 1:r){
arx[,i]=y[,i]-rep(colMeans(y)[i],each=T)
}
T=Ntime
bv=dis[1]; bw=dis[2]
Phi=matrix(m0,p,T)
Err=matrix(0,nrow=T,ncol=r)
###forward filtering
#m=list();C=list()#;e=list()
C=list()
m=list()
length(m)=length(C)=T
mt=m0
Ct=C0
for(t in (p+1):(T)){
F=arx[(t-1):(t-p),]
R=Ct/bv
Q=(t(F)%*%R%*%F)+ (Vi[[t-1]]*Sigma)
A=R%*%F%*%solve(Q) 
e=arx[t,]-t(F)%*%mt
mt=mt+A%*%e#G%*%mt+A%*%e
m[[t]]=mt
Ct=(R-A%*%Q%*%t(A))
Ct=(Ct+t(Ct))/2
C[[t]]=Ct
}
C[1:p]=list(C[[p+1]])#
m[1:p]=list(m[[p+1]])#
#C[p:(t+1)]=C[(p-1):t];m[p:(t+1)]=m[(p-1):t]
CC=abind(C,along=3)
mm=abind(m,along=2)
###backward sampling
h=1
Phi[,T]=mvrnorm(h,mm[,T],CC[,,T])##suppose to be T-distribution
for(t in (T-1):(p+1)){
mt=((1-bv)*mm[,t])+(bv*Phi[,(t+1)])##posterior
Err[t+1,]=arx[t+1,]-t(Phi[,(t+1)])%*%arx[t:(t-p+1),] #updated error	
Ct=(1-bw)*CC[,,t]
Phi[,t]=mvrnorm(h,mt,Ct)
}
for(t in (p+1):(p-1)){
Err[p+1,]=arx[p+1,]-t(Phi[,(p+1)])%*%arx[p:1,] #updated error
}
Phii=Phi
Erri=Err
list(Phii=Phii,Erri=Erri)
}
#res3=tvarm2(y,p,dis,m0,C0,V,Sigma)

####################
######## Metropolis algorithm ################
#logprior
MH=function(alphai,betai,design,Erri,Vi,p){
T=nrow(Erri)
r=ncol(Erri)
Sigma=corr.matrix(design,scales=betai)#myvar(design,alphai,betai)
LogPrior=0 
d=length(betai)
#betai=rep(1,8)
#LogP=function(betai){
#d=length(betai)
#LogPrior=0 
for(i in 1:d){
LogPrior=LogPrior-(betai[i]/4)-0.9*(1-exp(-betai[i]/4))	
}
#LogPrior
#}
#LogPrior(betai)
##loglikelihood
SigmaDet=det(Sigma)
SigmaInv=solve(Sigma)
Sigmainv=(SigmaInv + t(SigmaInv))/2
LogLike = -1.0 *(T -p)*log(SigmaDet)
for(t in (p+1):T){
	vt=Vi[t]
	e=t(Erri[t,]/sqrt(vt))
LogLike = LogLike-(e%*%Sigmainv%*%t(e))
}
LogJ=sum(log(betai))
#score=LogLike/2+LogPrior(betai)+LogJ
score=LogLike+LogPrior+LogJ
score
}
########
MH2=function(alphai,betai,design,Erri,Vi,p,scorei){
T=nrow(Erri)
r=ncol(Erri)
#proposal function
etabeta = log(betai) + 0.01*rnorm(length(betai),0,1)
beta=exp(etabeta)
alpha=alphai
score=MH(alphai=alpha,betai=beta,design,Erri,Vi,p)
acceptance=score-scorei
if(log(runif(1,0,1))< acceptance){
alpha=alpha;beta=beta;score=score
}else{
alpha=alphai;beta=betai;score=scorei
}
list(beta=beta,score=score,alpha=alpha)
}

#####################main
T=nrow(y)
r=ncol(y)
PHI=array(0,c(p,T,NIter))
V=array(0,c(T,NIter))
ERR=array(0,c(T,r,NIter))
# initialization of alpha and beta

Gibbs=function(y,Pinput,p,dis,m0,C0,s0,n0,design,NIter){
T=nrow(y)
r=ncol(y)
PHI=array(0,c(p,T,NIter))
V=array(0,c(T,NIter))
ERR=array(0,c(T,r,NIter))
# initialization of alpha and beta
Beta=array(1,c(Pinput,NIter))
Alpha=array(1.999,c(Pinput,NIter))
alphai=Alpha[,1]
betai=Beta[,1]  
for(k in 1:NIter){
print(k)
##corr matrix
Sigma=corr.matrix(design,scales=betai)
#forward filtering with unknown variance
res=tvarm1(y,p,dis,m0,C0,s0,n0,Sigma)###produce[CC,mm,nn,ss]
##sample variance
Vi=myV(res$nn,res$ss,p,dis)###produce [V]
###forward filtering/ backward sampling with estimated variance
ress=tvarm2(y,p,dis,m0,C0,Vi,Sigma)##produce [Phii, Erri]
V[,k]=Vi
PHI[,,k]=Phii=ress$Phii
ERR[,,k]=Erri=ress$Erri
####Metropolis algorithm
scorei=MH(alphai,betai,design,Erri,Vi,p)##produces [scorei]
for(w in 1:200){
hh=MH2(alphai,betai,design,Erri,Vi,p,scorei)##produces [scorei,betai]
betai=hh$beta;scorei=hh$score;alphai=hh$alpha
}
Beta[,k]=hh$beta
}    
list(V,PHI,ERR,Beta) 
}  
    	
gap=Gibbs(y,Pinput,p,dis,m0,C0,s0,n0,design,NIter)
###########################PREDICTION

pred=function(gap,newinput,design,y,newy){
V=apply(gap[[1]][,-(1:Nburnin)],1,mean)	
THETA=apply(gap[[2]][,,-(1:Nburnin)],2,rowMeans)
ERR=apply(gap[[3]][,,-(1:Nburnin)],2,rowMeans)
BETA=apply(gap[[4]][,-(1:Nburnin)],1,mean)
Thetai=THETA
T=ncol(Thetai)
r=nrow(design)
CMoutput=y
newinput=X2
newoutput=matrix(0,T,nrow(newinput))
newoutput[1:Ptvar,]=1	
Sigma=corr.matrix(design,design,scales=BETA)
SigmaInv=solve(Sigma)
rho=corr.matrix(newinput,design,scales=BETA)#rep(0,r)
for(t in (Ptvar+1):T){
#err = CMoutput(t, :) - t(Thetai[:, t]) * CMoutput[(t-1):(t-Ptvar-1),]
#Err[t+1,]=arx[t+1,]-t(Phi[,(t+1)])%*%arx[t:(t-p+1),]##t in (p+1):T
err = CMoutput[t,] - t(Thetai[, t]) %*% CMoutput[(t-1):(t-p),]
mut = (t(Thetai[,t])%*%newoutput[(t-1):(t-p),]) + c(t(rho)%*%SigmaInv%*%c(err))
sigma2t = V[t]*diag(diag(1,r) - t(rho)%*% SigmaInv %*% rho)
newoutput[t,]= c(mut) + sqrt(sigma2t)*rnorm(1,0,1)#sqrt(sigma2t)*rnorm(0, 1)
}
}
#plot(quantile(Output.Phi(2, :, (Nburnin+1):NIter), 0.5,3), '-.k');
