load("CONN")##mean
CON=CONN
load("CON1")##sd
load("ALL20")
load("ALL21")##sd
load("idyno/ALL")##idyno
load("idyno/ALL00")##idyno
N=100
s=75
ninp=3
ndes=75
ntime=96
s2=c(1,2,5)##sub, o2,nh4
h=rep(NA,N)
for(i in 1:length(CON)){
h[i]=nrow(CON[[i]])
}
h[which(h==96)]
##extract complete cases
id=which(h==96)
h2=CON[id]###
h3=h2[1:s]###select 75% of sample
h3a=h2[-(1:s)]
library(abind)
dat1=abind(h3,along=3)
dat1a=abind(h3a,along=3)
out1=dat1[,6,]##biomass
out2=apply(dat1[,7:8,],3,rowSums)##total particle
inp1=dat1[,s2,]##training
##test
outtest1=dat1a[,6,]###test_mass
outtest2=apply(dat1a[,6:7,],3,rowSums)###test_nparticle
inp2=dat1a[,s2,]###test
##transforming the outputs
logout1=log(out1)
logout2=log(out2)
logtest1=log(outtest1)
logtest2=log(outtest2)
#normalize input data
min1=apply(inp1,2,min)
max1=apply(inp1,2,max)
min2=apply(inp2,2,min)
max2=apply(inp2,2,max)
input1=inp1#matrix(NA,nrow=ndes,ncol=ninp)
input2=inp2#matrix(NA,nrow=ndes,ncol=ninp)
for(i in 1:ntime){
input1[i,,]=(inp1[i,,]-min1)/(max1-min1)
input2[i,,]=(inp2[i,,]-min1)/(max1-min1)
}
X=t(input1[1,,])
y=logout1
X2=t(input2[1,,])
y2=logtest1
###############initiliaze model
###########
Ntime=ntime
Ptvar=p=2
m0 = matrix(1,nrow=ninp,ncol=1)
C0 = 10*diag(1,ninp)
n0 = 1
s0 = 10
psi=rep(1,ninp)
dis=c(.9,.9)
# iterations: each iteration consists of 200 MH draws.
NIter = 500#1100
Nburnin =100#100
Pinput=ninp
r=75
Beta=array(1,c(Pinput,NIter))
Alpha=array(1.999,c(Pinput,NIter))
alphai=Alpha[,1]
betai=Beta[,1] 
design=X#d1
#########

###########
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
Sigma=mycor(X,X,betai)
##############
tvarm1=function(input,y,p,dis,m0,C0,s0,n0,Sigma){
input=input
arx=as.matrix(y)
r=ncol(y)
##scale the output
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
for(t in p:(T)){
F=(input[t-1,,])
#F=rbind((input1[t,,]),arx[(t):(t-1),])
G=diag(psi)
R[[t]]=G%*%C[[t-1]]%*%t(G)/bw
Q[[t]]=(t(F)%*%R[[t]]%*%F)+ (s[[t-1]]*Sigma)
A[[t]]=R[[t]]%*%F%*%solve(Q[[t]]) 
at=G%*%m[[t-1]]
e[[t]]=arx[t,]-t(F)%*%at#G%*%m[[t-1]]
m[[t]]=at+A[[t]]%*%e[[t]]
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
#C[p:(t+1)]=C[(p-1):t];m[p:(t+1)]=m[(p-1):t]
#s[p:(t+1)]=s[(p-1):t];n[p:(t+1)]=n[(p-1):t]
CC=abind(C,along=3)
mm=abind(m,along=2)
nn=abind(n,along=1)
ss=abind(s,along=1)
list(CC=CC,nn=nn,ss=ss,mm=mm)
}

#res=tvarm1(input,y,p,dis,m0,C0,s0,n0,Sigma)
#####variance
myV=function(nn,ss,p,dis){
bv=dis[1]
T=length(nn)
V=rep(0,T)
nt=nn[T]
dt=nn[T]*ss[T]
VT=rgamma(1,shape=nt/2,rate=dt/2)##inverse variance
V[T]=1/VT##variance
for(t in (T-1):(p)){
nt=(1-bv)*nn[t]
dt=nn[t]*ss[t]
#VT=(bv*VT)+rgamma(1,shape=nt/2,rate=dt/2)
VT=(bv*VT)+rgamma(1,shape=nt/2,scale=2/dt)	
V[t]=1/VT	
}
V[1]=V[p]
V=V
return(V)
}
#res2=myV(res$nn,res$ss,p,dis)
#Vi=res2
################var2
#Vi=myV(res$nn,res$ss,p,dis)
tvarm2=function(input,y,p,dis,m0,C0,Vi,Sigma){
input=input
arx=as.matrix(y)
r=ncol(y)
for(i in 1:r){
arx[,i]=y[,i]-rep(colMeans(y)[i],each=T)
}
T=Ntime
bv=dis[1]; bw=dis[2]
Theta=matrix(m0,ninp,T)
Err=matrix(0,nrow=T,ncol=r)
G=diag(psi)
###forward filtering
#m=list();C=list()#;e=list()
C=list()
m=list()
length(m)=length(C)=T
mt=m0
Ct=C0
for(t in (p):(T)){
F=(input[t-1,,])
#F=rbind((input1[t,,]),arx[(t-1):(t-p),])
#F=arx[(t-1):(t-p),]
R=G%*%Ct%*%t(G)/bw
Q=(t(F)%*%R%*%F)+ (Vi[[t-1]]*Sigma)
A=R%*%F%*%solve(Q) 
at=G%*%mt
e=arx[t-1,]-t(F)%*%at
mt=(at)+A%*%e#G%*%mt+A%*%e
m[[t]]=mt
Ct=(R-A%*%Q%*%t(A))
Ct=(Ct+t(Ct))/2
C[[t]]=Ct
}
#C[1:p]=list(C[[p+1]])#
#m[1:p]=list(m[[p+1]])#
C[[1]]=(C[[p]])#
m[[1]]=(m[[p]])#
#C[p:(t+1)]=C[(p-1):t];m[p:(t+1)]=m[(p-1):t]
CC=abind(C,along=3)
mm=abind(m,along=2)
###backward sampling
h=1
Theta[,T]=mvrnorm(h,mm[,T],CC[,,T])##
for(t in (T-1):(p)){
mt=((1-bw)*mm[,t])+(bw*Theta[,(t+1)])##posterior
Err[t+1,]=arx[t+1,]-t(Theta[,(t+1)])%*%F#arx[t:(t-p+1),] #updated error	
Ct=(1-bw)*CC[,,t]
Theta[,t]=mvrnorm(h,mt,Ct)
}
Err[1,]=Err[p,]
#for(t in (1):(p)){
#Err[t,]=arx[t,]-t(Theta[,(t)])%*%F#arx[p:1,] #updated error
#}
Thetai=Theta
Erri=Err
list(Thetai=Thetai,Erri=Erri)
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
dd=length(betai)
#betai=rep(1,8)
#LogP=function(betai){
#d=length(betai)
#LogPrior=0 
for(i in 1:dd){
LogPrior=LogPrior-(betai[i]/4)-0.9*(1-exp(-betai[i]/4))	
}
#LogPrior
#}
#LogPrior(betai)
##loglikelihood
SigmaDet=det(Sigma)##eqiivalent to log det
SigmaInv=solve(Sigma)
Sigmainv=(SigmaInv + t(SigmaInv))/2
LogLike = -1.0 *(T -p)*log(2.132446e-56+SigmaDet)
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
print("acceptance")
}else{
alpha=alphai;beta=betai;score=scorei
print("rejection")
}
list(beta=beta,score=score,alpha=alpha)
}

#####################main
T=nrow(y)
r=ncol(y)
THETA=array(0,c(ninp,T,NIter))
V=array(0,c(T,NIter))
ERR=array(0,c(T,r,NIter))
# initialization of alpha and beta

Gibbs=function(y,Pinput,p,dis,m0,C0,s0,n0,design,NIter){
T=nrow(y)
r=ncol(y)
THETA=array(0,c(ninp,T,NIter))
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
res=tvarm1(input,y,p,dis,m0,C0,s0,n0,Sigma)###produce[CC,mm,nn,ss]
##sample variance
Vi=myV(res$nn,res$ss,p,dis)###produce [V]
###forward filtering/ backward sampling with estimated variance
ress=tvarm2(input,y,p,dis,m0,C0,Vi,Sigma)##produce [Thetai, Erri]
V[,k]=Vi
THETA[,,k]=Thetai=ress$Thetai
ERR[,,k]=Erri=ress$Erri
####Metropolis algorithm
scorei=MH(alphai,betai,design,Erri,Vi,p)##produces [scorei]
for(w in 1:200){
hh=MH2(alphai,betai,design,Erri,Vi,p,scorei)##produces [scorei,betai]
betai=hh$beta;scorei=hh$score;alphai=hh$alpha
}
Beta[,k]=hh$beta
}    
list(V,THETA,ERR,Beta) 
}  
    	
gap=Gibbs(y,Pinput,p,dis,m0,C0,s0,n0,design,NIter)
###########################PREDICTION
pred=function(gap,newinput,design,y,newy){
	nsim=1000
V=apply(gap[[1]][,-(1:Nburnin)],1,mean)	
THETA=apply(gap[[2]][,,-(1:Nburnin)],2,rowMeans)
ERR=apply(gap[[3]][,,-(1:Nburnin)],2,rowMeans)
BETA=apply(gap[[4]][,-(1:Nburnin)],1,mean)
Thetai=THETA
T=ncol(Thetai)
CMoutput=y
newoutput=matrix(0,T,3)
newoutput[1:Ptvar,]=y2[1:Ptvar,]
for(t in (p):T){
newinput=input2[t,,]
r=nrow(newinput)
Sigma=corr.matrix(design,design,scales=BETA)
SigmaInv=solve(Sigma)
rho=corr.matrix(newinput,design,scales=BETA)#rep(0,r)

#err = CMoutput(t, :) - t(Thetai[:, t]) * CMoutput[(t-1):(t-Ptvar-1),]
#Err[t+1,]=arx[t+1,]-t(Phi[,(t+1)])%*%arx[t:(t-p+1),]##t in (p+1):T
#err = CMoutput[t,] - t(Thetai[, t]) %*% input1[t,,]#CMoutput[(t-1):(t-ninp),]
err=ERR[t,]
mut = (t(Thetai[,t])%*%newinput[,]) + c(t(rho)%*%SigmaInv%*%c(err))
sigma2t = V[t]*diag(diag(1,r) - t(rho)%*% SigmaInv %*% rho)
#newoutput[t,]= c(mut) + sqrt(sigma2t)*rnorm(nsim,0,1)
newoutput[t,]= mut+rep(colMeans(y2))+ mvrnorm(1,mut,diag(sqrt(sigma2t)))
}
}

 cbind(y2,newoutput)
 cor(y2,newoutput)^2
          [,1]      [,2]      [,3]
[1,] 0.9193788 0.8946327 0.8298150
[2,] 0.9152627 0.8952615 0.8214542
[3,] 0.9186043 0.8946806 0.8282269

#plot(quantile(Output.Phi(2, :, (Nburnin+1):NIter), 0.5,3), '-.k');
