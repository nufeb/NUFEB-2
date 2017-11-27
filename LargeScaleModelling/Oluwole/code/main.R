para=c("KsHET","Ko2HET","Kno2HET","Kno3HET","Knh4AOB","Ko2AOB","Kno2NOB","Ko2NOB","MumHET","MumAOB","MumNOB","etaHET","bHET","bAOB","bNOB","bEPS","YEPS","YHET","EPSdens","EPSratio","factor","ke","time")

val=c(0.01,0.81,0.0003,0.0003,0.001,0.0005,0.0013,0.00068,0.00006944444,0.00003472222,0.00003472222,0.6,0.00000462962,0.00000127314,0.00000127314,0.00000196759, 0.18,0.61,30,1.25,1.5,5e+10)

set.seed(3)
library(lhs)
nvar=22
nob=100
npoints=50
setwd("/frontiers-shared/shared/data/seg-dat/noo11/test2b")
load("dat")

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
lm_pred=predict(mod,newdata=DM_new,se=TRUE)
lm_pred=cbind(lm_pred$fit,lm_pred$se.fit)
Obs=mod$residual

################my GP method
DM=DM[1:100,22:23]####log
DM_new=DM_new[,22:23]
source("GP_code1.R",echo=FALSE)
GP_pred=prediction(DM,Obs[1:100],DM_new)
emulex=cbind((lm_pred[,1]+GP_pred[,1]),(lm_pred[,2]+2*GP_pred[,2]))
lammp=test[,24]*10^6;emu=emulex[,1]
prop=1-(sum((lammp-emu)^2)/sum((lammp- mean(lammp))^2))



#opts(DM, scales=scales_start,Obs[1:100],RB)$par
#optimal.scales(DM, scales=scales_start,d=Obs[1:100],func=RB)
B=diag(scales,nrow=length(scales))
#corr_matrix=exp(-(as.matrix(dist(x[,1:7])))^2*B)###correct
###############################################################12-10-2015
##########kriging
s2=sample(1:nrow(DM),1500,FALSE)
 mod2=km(~.,design=DM[s2,],response=Y[s2],covtype="gauss",nugget.estim=TRUE,
estim.method="MLE",optim.method="BFGS",control=list(po.size=20,trace=TRUE))

#s3=sample(1:nrow(DM),2500,FALSE)
 mod3=km(~.,design=DM,response=Y,covtype="gauss",nugget.estim=TRUE,
estim.method="MLE",optim.method="BFGS",control=list(pop.size=50,trace=TRUE))
mod33=km(~.,design=DM[,13:16],response=Y,covtype="gauss",nugget.estim=TRUE,
estim.method="MLE",optim.method="BFGS",control=list(pop.size=50,trace=TRUE))
#pp=predict(mod3,newdata=data.frame(test),se.compute=TRUE,type="UK")
pp=predict(mod3,newdata=data.frame(DM_new),se.compute=TRUE,type="UK")
lammp=test[,24]*10^6;emu=list(pp$mean,pp$lower95,pp$upper95)
prop=1-(sum((lammp-emu)^2)/sum((lammp- mean(lammp))^2))

###############plot
id=1:1250
#id=645:695
lammps=lammp[id]
emus=emu[[1]][id]
U=emu[[3]][id]
L=emu[[2]][id]
ind=order(lammps)
x1 = min(L)
x1 = min(x1,lammps)
x2 = max(U)
x2 = max(x2, lammps)

pdf("myplot3")
par(mar=c(5.1,5.4,4.1,2.1))
plot(lammps,lammps,main="Mean floc diameter for Lammp Vs Emulator",xlab=expression("Observed"~(10^{-6})),ylab=expression("Predicted"~(10^{-6})),lty=1,col="black",cex.lab=1.6,cex.axis=1.6,cex.main=1.6,xlim=c(x1,x2),ylim=c(x1,x2))
lines(lammps[ind],emus[ind], lty = 1,col="green")
lines(lammps[ind],U[ind], col="red",lty=1)
lines(lammps[ind],L[ind], col="red",lty=1)
legend(.6,max(pp$upper95),c("Predicted","Observed","95% C.I"),fill=c("green","black","red"),bty="0",border="black",cex=1.75,horiz=FALSE)
dev.off()
#############plot as a function of time
h=20
gg=dat[h,,1:23];colnames(gg)=para####at design point5
pp2=predict(mod3,newdata=data.frame(gg[,1:23]),se.compute=TRUE,type="UK")
lammp=dat[h,,24]*10^6;time=gg[,23];emu=list(pp2$mean,pp2$lower95,pp2$upper95)
id=1:100
lammps=lammp[id]
emus=emu[[1]][id]
U=emu[[3]][id]
L=emu[[2]][id]
ind=order(emus)
x1 = min(L)
x1 = min(x1,lammps)
x2 = max(U)
x2 = max(x2,lammps)

pdf("myplot5")
par(mar=c(5.1,5.4,4.1,2.1))
plot(time,lammps,main="Mean floc diameter for Lammp Vs Emulator",xlab=expression("Time"~(seconds)),ylab=expression("Predicted"~(10^{-6})),lty=1,col="black",cex.lab=1.6,cex.axis=1.6,cex.main=1.6,xlim=range(time),ylim=c(x1,x2)+c(0,.03))

for(i in 1:length(lammps)){
lines(rep(time[i],2),c(L[i],U[i]), col="red")
}
points(time,emus,col="green")
legend(.9,max(pp2$upper95)+.04,c("Predicted","Observed","95% C.I"),fill=c("green","black","red"),bty="0",border="black",cex=1.75,horiz=FALSE)

dev.off()
#####################sensitivity analysis
set.seed(12)
library(sensitivity)
library(lhs)
mm=function(Xnew,m)
predict.km(m,Xnew,"UK",se.compute=FALSE,checknames=FALSE)$mean

nvar=23
npoints=2000
hyper=maximinLHS(npoints,nvar, dup=2)
design= matrix(NA,nrow=npoints, ncol=nvar)
for(j in 1:(nvar-1)){
design[,j]= qunif(hyper[,j], min=0.8*val[j], max=1.2*val[j])
}
tt=qunif(hyper[,23],min=min(time),max=max(time))###time
design[,23]=tt
X1=data.frame(design)
hyper=maximinLHS(npoints,nvar, dup=2)
design= matrix(NA,nrow=npoints, ncol=nvar)
for(j in 1:(nvar-1)){
design[,j]= qunif(hyper[,j], min=0.8*val[j], max=1.2*val[j])
}
tt=qunif(hyper[,23],min=min(time),max=max(time))###time
design[,23]=tt
X2=data.frame(design)

#X1=data.frame(matrix(runif(npoints*nvar),nrow=npoints))
#X2=data.frame(matrix(runif(npoints*nvar),nrow=npoints))
names(X1)=names(X2)=para
#sens=sobol2007(model=mm,X1,X2,nboot=5,type="UK",mod3=mod3)
sens2=sobolGP(model=mod3,type="UK",MCmethod="sobol2007",X1,X2,nboot=100,nsim=100,candidate)

mod4=km(~.^2,design=DM,response=Y,covtype="matern5_2",nugget.estim=TRUE,
estim.method="MLE",optim.method="BFGS",control=list(pop.size=50,trace=TRUE))
mod3b=km(~.,design=DM,response=Y,covtype="matern5_2",nugget.estim=TRUE,
estim.method="MLE",optim.method="BFGS",control=list(pop.size=50,trace=TRUE))
save(mod3b,file=paste(dirname,"mod3b",sep="/"))
save(mod4,file=paste(dirname,"mod4",sep="/"))

mm2=function(Xnew,m)
predict.km(m,Xnew,"UK",se.compute=FALSE,checknames=FALSE)$mean
sens1d=sobolGP(mod3,X1=X1,X2=X2,nboot=100)
sens2d=sobolGP(mod3b,X1=X1,X2=X2,nboot=100)
sens1e=sobolGP(mod3b,type="UK",X1=X1,X2=X2,nboot=100)

#############
sens2a=soboljansen(model=mm,X1=X1,X2=X2,nboot=1000,m=mod3)
#sens2b=soboljansen(model=mm,X1=X1,X2=X2,nboot=100,m=mod4)
sens2b=soboljansen(model=mod,X1=X1,X2=X2,nboot=1000)##linear fit
sens2c=soboljansen(model=mm,X1=X1,X2=X2,nboot=1000,m=mod3b)
save(sens2a,file=paste(dirname,"sens2a",sep="/"))
cere1=data.frame(sens2a[11])##mod3,linear Guassian cov
cere2=data.frame(sens2b[11])####mod4,quadratic; Gaussian cov
cere3=data.frame(sens2c[11])####mod3b,Linear, Matern5_2

x1=pmax(cere1[1:23,1],0)
x2=pmax(cere2[1:23,1],0)
x3=pmax(cere3[1:23,1],0)
ddd=cbind(x1,x2,x3)
par(mar=c(5.1,4.5,0,.5)+0.0)
xvals=barplot(t(ddd)
,ylab="Total sensitivity",main="Sensitivity analysis",offset=0,cex.lab=1.2,col=c("green","blue","pink"),density=c(40,60,90),cex.main=1.9,space=0.3,angle=c(45,90,135),cex.axis=1.7,cex.names=1.7,args.legend=list(x=c(1,8),y=c(0.8,1.2),cex=1.75),legend=c("Linear+Gaussian cov","Quadratic+Gaussian","Linear+matern5_2"))
text(xvals,par("usr")[3]-.03,srt=45,adj=1,labels=para,xpd=TRUE,cex=1.8)
