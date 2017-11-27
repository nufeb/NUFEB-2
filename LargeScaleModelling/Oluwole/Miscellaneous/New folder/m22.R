#library(doMPI)
#cl=startMPIcluster()
#registerDoMPI(cl)

####
set.seed(1)
load("ALL2")
ALL3=ALL2
load("alpha3")##nutrient
load("design")
ndes2=280
ALL2=ALL2[1:ndes2]
alpha3=alpha3[1:ndes2]
pow=c(1e+16,1e+06,1,1,1)###to make noise positive
para=c("KsHET","Ko2HET","Kno2HET","Kno3HET","Knh4AOB","Ko2AOB","Kno2NOB","Ko2NOB","MumHET","MumAOB","MumNOB","etaHET","bHET","bAOB","bNOB","bEPS","YHET","YAOB","YNOB","YEPS","Y1","Do2","Dnh4","Dno2","Dno3","Ds","diffT")
nvar=27
nout=5
ndes=4500
library(DiceKriging)
library(MuFiCokriging)
library(abind)
imp2=c(7:11,13)##c(nHET,nAOB,nNOB,nEPs,nDEAD,nparticle)
imp1=c(15,3,12,13,1)##c(D,eqv.dia,fractal,nparticle,mass)
alpha1=list();alpha2=list()
for(i in 1:length(ALL2)){
alpha1[[i]]=ALL2[[i]][-1,,]#### response==one step ahead
}
for(i in 1:length(ALL2)){
alpha2[[i]]=ALL2[[i]][-nrow(ALL2[[i]]),,] #####initial value
}
############nutrients
gg3=gg2=gg1=rep(NA,ndes2)
for(i in 1:ndes2){
gg3[i]=nrow(alpha3[[i]])
}
for(i in 1:ndes2){
gg2[i]=nrow(alpha2[[i]])
}
for(i in 1:ndes2){
gg1[i]=nrow(alpha1[[i]])
}
mx=pmin(gg1,gg2,gg3)
X=Y=Z=list()
for(i in 1:ndes2){
X[[i]]=alpha1[[i]][1:mx[i],,]
Y[[i]]=alpha2[[i]][1:mx[i],,]
Z[[i]]=alpha3[[i]][1:mx[i],]
}
alpha1=X;alpha2=Y;alpha3=Z
####
bb1=abind(alpha1,along=1)##response
bb2=abind(alpha2,along=1)##initial==y0
aa1=abind(alpha3,along=1)
bb3=bb1[,,1]#mean response
bb4=bb2[,,1]#mean initial
bb5=bb2[,,2]##variance==nugget/noise var
jj=which(is.na(bb5[,11]))###find index for the uncompleted simulation (time<73)
bb3=bb3[-jj,]
bb4=bb4[-jj,]
bb5=bb5[-jj,]
###remove extreme values from floc mass
#p=which(bb3[,1]>8e-10)
p=which(bb3[,1]>1e-10)
length(p)
bb3=bb3[-p,]
bb4=bb4[-p,]
bb5=bb5[-p,]
#aa1=aa1[-jj,]
##log transformed the floc mass and number of particles
bb3[,c(1,13)]=log(bb3[,c(1,13)])
bb4[,c(1,13)]=log(bb4[,c(1,13)])
#design2=design[-t,]##remove uncompleted cases design
design2=design[1:ndes2,]
des=list()
for(i in 1:length(alpha1)){
des[[i]]=matrix(rep(design2[i,],nrow(alpha1[[i]]),nrow=nrow(alpha1[[i]])),ncol=nvar,byrow=TRUE)
}
des2=abind(des,along=1)
des2=cbind(des2,aa1[,-1])####combine the nutrients with parameters
#des2=des2[-jj,]
###normalized/standardized the input "parameters"
mmin=apply(des2,2,min)
mmax=apply(des2,2,max)
center=mmin
scale=(mmax-mmin)
des3=scale(des2,center=center,scale=scale)
des3=des3[-jj,]
des3=des3[-p,]
output=bb3[,imp1]#5 outputs
s=sample(1:nrow(output),ndes)
out=matrix(NA,nrow=length(s),ncol=nout)
ini=out
bb4a=bb4[,imp1]###choose initial for the 5 outputs
for(i in 1:nout){##scale the outputs
out[,i]=cbind(output[s,i])#*pow[i])
ini[,i]=cbind(bb4a[s,i])#*pow[i])
}
input=cbind(des3[s,],ini)
noise=bb5[,imp1]######select the noise variance for just the 5 choosen outputs
nam3=c("y1","y2","y3","y4","y5")
nam2=c("s","o2","no2","no3","nh4")
inp=input
##########noise transformation using Boukouvalas (2009) et al approach
#nrep=10
#d=nrep-1
#nois=t(noise[s,])*pow^2
#h1=log(t(nois))
#h2=replace(h1,which(h1<0),0)
#r=h2+ (d+(log(2)*d)-digamma(d/2))^(-1)
#sigs=trigamma(d/2)
nug=list()
for(i in 1:nout){
nug[[i]]=noise[s,i]####noise variance/nugget values
}
##
colnames(inp)=c(para,nam2,nam3)
load("polyset")
f=list();resp=list()
for(i in 1:nout){
quadterms=polySet((length(para)+6),1,1,1)
quad1 <-formula(paste("~", makeScope(quadterms,c(para,nam2,nam3[i]))))
f[[i]]=quad1
}
for(i in 1:nout){
resp[[i]]=out[,i]
}
rownames(inp)=NULL
nest=NestedDesign(inp,nlevel=nout,n=rep(length(s),nout))
rm(bb1,bb2,bb3,bb4,bb5,noise,nois,output,X,Y,Z,alpha1,alpha2,alpha3)
library(doParallel)
cl=makeCluster(nout)
registerDoParallel(cl)

dirname=paste(getwd(),"result22",sep="/")
save(scale,file=paste(dirname,"scale",sep="/"))###spa
save(center,file=paste(dirname,"center",sep="/"))###spa
#save(pow,file=paste(dirname,"pow",sep="/"))###sp

modd=foreach(g=1:nout,.verbose=FALSE,.errorhandling="stop") %dopar% {
library(DiceKriging)
library(MuFiCokriging)
km(formula=f[[g]],nest$PX,response=resp[[g]],covtype="exp",control=list(multistart=nout),noise.var=nug[[g]])
}
save(modd,file=paste(dirname,"modd",sep="/"))###spa
#####sensitivity
#load("modd")
library(sensitivity)
n <- 1000
nvar=27+5;nout=5
nam=c(para,nam2,nam3)
X1 <- data.frame(matrix(runif((nvar+nout) * n), nrow = n))
X2 <- data.frame(matrix(runif((nvar+nout) * n), nrow = n))
colnames(X1)=colnames(X2)=nam

sens=list()
for(i in 1:nout){
sens[[i]]=sobolGP(modd[[i]],X1=X1,X2=X2,type="UK")
}
save(sens,file=paste(dirname,"sens",sep="/"))###spa
stopCluster(cl)

