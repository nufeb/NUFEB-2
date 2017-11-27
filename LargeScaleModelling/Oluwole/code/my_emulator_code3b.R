 #####################Emulator
#setwd("C:/Users/olu/Desktop/NUFEB_materials/Long-Run")
fam=c("id","type","Dia","X","Y","Z","Vx","Vy","Vz","Fx","Fy","Fz")#####output names
dd=list()
r=list()
#nob=100
nfac=12
ww=1:50
dir="/data/noo11/train/test2"
#mydat[[3]][[120]]
library(doParallel)
cl <- makeCluster(length(ww))
registerDoParallel(cl)
mydat<-{foreach(nam=1:length(ww),.verbose=T,.errorhandling="stop") %dopar% {
L=paste(dir,"inputnam",sep="/")
setwd(sub("nam",nam,L))
nob=length(list.files())-6
for (i in 1:nob){
dd[i]=list(read.table(sub("i",i,(paste("data","i.xlsx",sep=""))),sep=",",col.names=fam,header=FALSE))
r[i]=list(nrow(dd[[i]]))
result=list(dd,r)
}
result
}}
stopCluster(cl)
save(mydat,file=paste(dir,"mydata",sep="/"))
#mydat[[1]][[2]][[2]]##50*2*var
## (mydat[[50]][[1]][[1]][[13]])##50*2*var*12
##nrow(mydat[[35]][[1]][[93]][1:12])##important to get no of particles
(mydat[[i]][[1]][[j]][[]])
################################09-09-2015
mydat1=mydat1_var=array(NA,c(50,100,12))
for( i in 1:50){
for(j in 1:100){
mydat1[i,j,]=colMeans(mydat[[i]][[1]][[j]][1:12])####diameter mean
}}
for( i in 1:50){
for(j in 1:100){
mydat1_var[i,j,]=apply(mydat[[i]][[1]][[j]][1:12],2,var)####diameter variance
}}
dia_m=mydat1[,,3]####diameter mean
dia_v=mydat1_var[,,3]##diameter variance
dia_m=array(dia_m,c(50,1,100))
dia_v=array(dia_v,c(50,1,100))
##########compute average and variance at each time-step
dd2=matrix(NA,nrow=nob,ncol=nfac)
dd2_var=matrix(NA,nrow=nob,ncol=nfac)
for(i in 1:nob){
dd2[i,]=rbind(colMeans(dd[[i]])) ##column mean
dd2_var[i,]=rbind(apply(dd[[i]],2,var))###column variance
colnames(dd2)=fam
}

#######################INPUT1
para=c("KsHET","Ko2HET","Kno2HET","Kno3HET","Knh4AOB","Ko2AOB","Kno2NOB","Ko2NOB","MumHET","MumAOB","MumNOB","etaHET","bHET","bAOB","bNOB","bEPS","YEPS","YHET","EPSdens","EPSratio","factor","ke","time")

val=c(0.01,0.81,0.0003,0.0003,0.001,0.0005,0.0013,0.00068,0.00006944444,0.00003472222,0.00003472222,0.6,0.00000462962,0.00000127314,0.00000127314,0.00000196759, 0.18,0.61,30,1.25,1.5,5e+10)

set.seed(3)
library(lhs)
nvar=22
npoints=50
hyper=maximinLHS(npoints,nvar, 2)
design= matrix(NA,nrow=npoints, ncol=nvar)
for(j in 1:nvar){
design[,j]= qunif(hyper[,j], min=0.8*val[j], max=1.2*val[j])
}
design2=array(rep(design,100),c(50,22,100))
#######Read in the input data;sequence of time-steps
#input=seq(0,352000,2000)
input=t(matrix(rep(seq(0,198000,2000),50),ncol=50))
input=array(input,c(50,1,100))
library(abind)
inp=abind(design2,input,along=2)
dat=abind(inp,dia_m,dia_v,along=2)
dirname=getwd()
save(dat,file=paste(dirname,"dat",sep="/"))
dat=aperm(dat,c(1,3,2))
set.seed(10)
s=sample(1:nob,75,FALSE)#####sample 150 observations, leave out 25 observations for CV

train=dat[,s,]
test=dat[,-s,]
train=array(train,c(50*length(s),25))
test=array(test,c(50*25,25))
obs=train[,24]*10^6 #######emulation of average diameter
obs2=train[,25]*10^14    #######emulation of variance diameter
input=as.data.frame(train[,1:23])
names(input)=para

Y=obs
DM=input
DM_new=as.data.frame(test[,1:23])
names(DM_new)=para
#mod=lm(Y~log(DM))
mod=lm(Y~.,data=DM)
lm_pred=predict(mod,newdata=DM_new)
Obs=mod$residual

#mm=mlegp(X=(DM),Z=Obs,constantMean=1)##
#dd=predict(mm,newData=as.matrix(DM_new),se.fit=TRUE)
#jitter(Obs,amount=1e-6)######to prevent numerical problems
#interpolant(DM_new[200,],d=Obs[1:100],xold=DM,Ainv=NULL,func=RB,scales=c(1,1),give.full.list=TRUE)
#interpolant.quick(DM_new,d=Obs[1:100],xold=DM,Ainv=NULL,func=RB,scales=c(1,1)) 
kala=jitter(Obs,amount=1e-6)
################my GP method
DM=DM[1:100,22:23]####log time
DM_new=DM_new[,22:23]
source("GP_code1.R",echo=FALSE)
GP_pred=prediction(DM,Obs[1:100],DM_new)

emulex=list((lm_pred+GP_pred[,1]),GP_pred[,2])
lammp=test[,24]*10^6;emu=emulex[[1]]
prop=1-(sum((lammp-emu)^2)/sum((lammp- mean(lammp))^2))
source("myplot.R",echo=FALSE)
myplot(observed=lammp,predicted=emulex[[1]],Z=emulex[[2]],sds=2,band=TRUE)

myplot(observed=lammp[1:20],predicted=emulex[[1]][1:20],Z=emulex[[2]][1:20],sds=2,band=FALSE)

#mm=mlegp(X=(DM[1:100,23]),Z=Obs[1:100],constantMean=1)##
#dd=predict(mm,newData=as.matrix(DM_new),se.fit=TRUE)

###############plot
#id=1045:1075
id=645:675
lammps=lammp[id]
emus=emulex[[1]][id]
Z=emulex[[2]][id]
U=emus+2*Z
L=emus-2*Z
ind=order(lammps)
x1 = min(L)
x1 = min(x1,lammps)
x2 = max(U)
x2 = max(x2, lammps)

pdf("plot4")

plot(lammps,emus,main="Plot of mean diameter for Lammp Vs Emulator with 95% C.I",xlab="Observed X exp(-6)",ylab="Predicted X exp(-6)",lty=1,col="green",cex.lab=1.6,cex.axis=1.6,cex.main=1.6,xlim=c(x1,x2),ylim=c(x1,x2))
lines(lammps[ind],lammps[ind], lty = 1)
lines(lammps[ind],U[ind], col="red",lty=1)
lines(lammps[ind],L[ind], col="red",lty=1)
legend(.6,1.3,c("Predicted","Observed","95% C.I"),fill=c("green","black","red"),bty="0",border="black",cex=1.75,horiz=FALSE)
dev.off()
##########################variance

train=dat[,s,]
test=dat[,-s,]
train=array(train,c(50*length(s),25))
test=array(test,c(50*25,25))
obs=train[,24]*10^6 #######emulation of average diameter
obs2=train[,25]*10^14    #######emulation of variance diameter
input=as.data.frame(train[,1:23])
names(input)=para

Y=obs2
DM=input
DM_new=as.data.frame(test[,1:23])
names(DM_new)=para
#mod=lm(Y~log(DM))
mod=lm(Y~.,data=DM)
lm_pred=predict(mod,newdata=DM_new)
Obs=mod$residual

#mm=mlegp(X=(DM),Z=Obs,constantMean=1)##
#dd=predict(mm,newData=as.matrix(DM_new),se.fit=TRUE)
################my GP method
DM=DM[1:100,23]####log time
DM_new=DM_new[,23]
source("GP_code1.R",echo=FALSE)
GP_pred=prediction(DM,Obs[1:100],DM_new)

emulex=list((lm_pred+GP_pred[,1]),GP_pred[,2])
lammp=test[,25]*10^14;emu=emulex[[1]]
prop=1-(sum((lammp-emu)^2)/sum((lammp- mean(lammp))^2))
source("myplot.R",echo=FALSE)
myplot(observed=lammp,predicted=emulex[[1]],Z=emulex[[2]],sds=2,band=TRUE)

myplot(observed=lammp[1:20],predicted=emulex[[1]][1:20],Z=emulex[[2]][1:20],sds=2,band=FALSE)

#mm=mlegp(X=(DM[1:100,23]),Z=Obs[1:100],constantMean=1)##
#dd=predict(mm,newData=as.matrix(DM_new),se.fit=TRUE)

###############plot
###############plot
id=1045:1050
#id=645:675
lammps=lammp[id]
emus=emulex[[1]][id]
Z=emulex[[2]][id]
U=emus+2*(Z)
L=emus-2*(Z)
ind=order(lammps)
x1 = min(L)
x1 = min(x1,lammps)
x2 = max(U)
x2 = max(x2, lammps)

pdf("plot1")

plot(lammps,emus,main="Plot variance diameter of particle with 95% C.I",xlab="Observed X exp(-14)",ylab="Predicted X exp(-14)",lty=1,col="green",cex.lab=1.6,cex.axis=1.6,cex.main=1.6,xlim=c(x1,x2),ylim=c(x1,x2))
lines(lammps[ind],lammps[ind], lty = 1)
lines(lammps[ind],U[ind], col="red",lty=1)
lines(lammps[ind],L[ind], col="red",lty=1)
legend(.6,1.3,c("Predicted","Observed","95% C.I"),fill=c("green","black","red"),bty="0",border="black",cex=1.75,horiz=FALSE)
dev.off()

############16-09-2015
plot(lammps,emus,xlim=range(lammps),pch=19,xlab="Observed X exp(-14)",ylab="Predicted X exp(-14)")
polygon(c(lammps,rev(lammps),lammps[1]),c(U,rev(L),U[1]),col="gray",border="NA")
lines(lammps,emus,lwd=2)
points()


#########################22-09-2015####read in Lammp output directly
index=grep("ITEM: ATOMS id type diameter x y z vx vy vz fx fy fz",readLines("snapshot.bubblemd"),value=FALSE,fixed=TRUE)
ind=length(readLines("snapshot.bubblemd"))
#ind=countLines("snapshot.bubblemd")
dd=list()
index=c(index,ind)
for (i in 1:length(index)){
dd[i]=list(read.csv("snapshot.bubblemd",sep="",skip=index[i]-1,header=TRUE,colClasses=numeric(),nrows=index[i+1]-index[i]-9))
}
##############
