ALL=ALL0=ALL1=list()
npoints2=5
#npoint=100
#for (j in 1:npoints){
#setwd(sub("j",j,paste(getwd(),"inputj",sep="/")))
path=getwd()
for (num in 1:npoints2){
setwd(sub("num",num,paste(getwd(),"inputnum",sep="/")))
load("dd")
ncomp2=40
ncomp=length(dd)
comp2=comp=matrix(NA,nrow=ncomp,ncol=4)
nparticle=rep(NA,ncomp)
mas=nparticle
D=rep(NA,ncomp)
#seg=rep(NA,ncomp)
A=rep(NA,ncomp);V=rep(NA,ncomp)
Frac=rep(NA,ncomp)

####no of particle
for(j in 1:ncomp){
nparticle[j]=nrow(dd[[j]])####natom
}

##########################Biofilm composition
for(j in 1:ncomp){
for (type in c(1:4)){
gdat=dd[[j]]
comp2[j,type]=length(which(gdat[,1]==type))/nparticle[j]
comp[j,type]=length(which(gdat[,1]==type))#/nparticle[j]
}}
##########diversity
for (j in 1:ncomp){
N=nparticle[j]
D[j]=1-sum((comp[j,])*(comp[j,]-1))/(N*(N-1))
}
##########biofilm thickness/height and surface roughness
x=c(0.000000e-04,1.3600000e-04)
y=c(0.000000e-04,1.360000e-04)
z=c(0.000000e-04,1.3600000e-04)
Lx=x[2]-x[1]
Ly=y[2]-y[1]
dmax=10##no of partitions cells
delta=0.5*(max(x)-min(x))/dmax##center of first cell along x
delta2=0.5*(max(y)-min(y))/dmax##center of first cell along y
delta3=0.5*(max(z)-min(z))/dmax##center of first cell along z

xx=seq(delta,max(x),by=2*delta)
yy=seq(delta2,max(y),by=2*delta2)
zz=seq(delta3,max(z),by=2*delta3)
grid=expand.grid(xx,rev(yy))
grid2=expand.grid(xx,rev(yy),(zz))
bl=list()
surf=list()
gh2=list()
for(k in 1:ncomp){
dat=as.matrix(dd[[k]][,3:5])
R=R2=rep(NA,nrow(grid))
P=dat
#gh=list()
for(i in 1:nrow(grid)){
##find occupied blocks
gh=which((P[,1]>(grid[i,1]-delta)&P[,1]<(grid[i,1]+delta))&(P[,2]>(grid[i,2]-delta2)&P[,2]<(grid[i,2]+delta2)))
bl[[i]]=grid2[which(grid2[,1]==grid[i,1]&grid2[,2]==grid[i,2]),]##extract blocks above baseline block
gd=bl[[i]]
for(j in 1:nrow(gd)){
gh2[[j]]=length(which((P[,1]>(gd[j,1]-delta)&P[,1]<(gd[j,1]+delta))&(P[,2]>(gd[j,2]-delta2)&P[,2]<(gd[j,2]+delta2))&(P[,3]>(gd[j,3]-delta3)&P[,3]<(gd[j,3]+delta3))))
}
ind=which(gh2!=0)

R[i]=max(dat[gh,3])#maximum z-values==max height of biofilms
R2[i]=max(gd[ind,3])
}
h=R[which(is.finite(R))]
hh=sum(h)/(length(h))###biofilms mean surface height
w=sqrt(sum((h-hh)^2)/(length(h)))##surface roughness
#surf[[k]]=c(hh,w)
##h2
h2=R2[which(is.finite(R2))]
hh2=sum(h2)/(length(h2))###biofilms mean surface height
w2=sqrt(sum((h2-hh2)^2)/(length(h2)))##surface roughness
surf[[k]]=c(hh,w,hh2,w2)
}
library(abind)
surface=abind(surf,along=0)
res=cbind(comp2,D,nparticle,surface)
#res=cbind(res,res2)
colnames(res)=c("comp1","comp2","comp3","comp4","div","nparticle","h1","w1","h2","w2")
write.table(res,file=sub("num",num,"outputnum.csv"))
#write.table(d1,file=sub("num",num,"shear_datnum.csv"))
ALL[[num]]=res
library(reshape2)
b1=do.call(rbind,lapply(ALL,melt))
b2=acast(b1,b1$Var1~b1$Var2,mean)
b3=acast(b1,b1$Var1~b1$Var2,sd)
ALL0=b2
ALL1=b3
setwd("../")
}
save(ALL,file=paste(path,"ALL",sep="/"))
save(ALL0,file=paste(path,"ALL0",sep="/"))
save(ALL1,file=paste(path,"ALL1",sep="/"))

#setwd("../")
#}
