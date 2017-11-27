#######Dicekriging emulator
load("alldat")
npart=4###(diameter,x, y z)
natom=5###AOB,NOB,EPS, HET, inert
nvar=23
nob=40###100
npoints=jj=50##1000
runtime=86400###352000
step=2000
para=c("KsHET","Ko2HET","Kno2HET","Kno3HET","Knh4AOB","Ko2AOB","Kno2NOB","Ko2NOB","MumHET","MumAOB","MumNOB","etaHET","bHET","bAOB","bNOB","bEPS","YEPS","YHET","EPSdens","EPSratio","factor","ke","time")

val=c(0.01,0.81,0.0003,0.0003,0.001,0.0005,0.0013,0.00068,0.00006944444,0.00003472222,0.00003472222,0.6,0.00000462962,0.00000127314,0.00000127314,0.00000196759, 0.18,0.61,30,1.25,1.5,5e+10)

set.seed(10)
library(DiceKriging)
library(abind)
#train=array(NA,c(jj,length(s),natom,npart))
train=list()
hh=nrow(alldat[[1]])#####be careful
for(num in 1:jj){
s=sample(1:nrow(alldat[[num]]),nob,FALSE)#####sample 100 time steps
train[[num]]=alldat[[num]][s,,]
doll=abind(train,along=1)
#train=array(train,c(jj*length(s),natom,npart))
}
train=doll[,,1]###mean diameter only/5 particles

###############input
library(lhs)
hyper=maximinLHS(npoints,nvar, dup=2)
design= matrix(NA,nrow=npoints, ncol=nvar)
for(j in 1:(nvar-1)){
design[,j]= qunif(hyper[,j], min=0.8*val[j], max=1.2*val[j])
}

input=array(rep(design,length(s)),c(dim(design),length(s)))
input=aperm(input,c(1,3,2))
input=array(input,c(npoints*length(s),nvar))
#input[,23]=rep(seq(0,runtime,step)[s],each=npoints)
cap=seq(0,runtime,step)[s]
cap=rep(cap,each=npoints)
input[,23]=cap
#t(matrix(rep(seq(0,runtime,step),npoints),ncol=npoints))
#################CV data
test=abind(alldat,along=1)[,,1]###mean diameter only/ for 5 particles
#input2=aperm(array(rep(design,nrow(test)),c(dim(design),nrow(test))),c(1,3,2))
#input2=array(input2,c(nrow(test),nvar))
m1=array(rep(design,ncol(design)*npoints),c(dim(design),hh))
m2=aperm(m1,c(1,3,2))
#m2=aperm(m1,c(3,1,2))
input2=array(m2,c(npoints*hh,nvar))
#input2[,23]=rep(seq(0,runtime,step),each=npoints)
time=matrix(NA,nrow=npoints,ncol=hh)
for(num in 1:npoints){
#time[[num]]=seq(0,(nrow(alldat[[num]])-1)*step,step)
time[num,]=seq(0,(nrow(alldat[[num]])-1)*step,step)
}
input2[,23]=c((time))
##########################################


m1=array(rep(design,hh),c(dim(design),hh))
m2=aperm(m1,c(3,1,2))
test=abind(alldat,along=1)
test2=array(test,c(hh,jj,5,4))
test=test2[,,,1]
input2=abind(test,m2,along=3)
cap=seq(0,runtime,step)
cap=array(rep(cap,each=jj),c(jj,hh,1))
cap=aperm(cap,c(2,1,3))
input2[,,28]=cap####time
s=sample(1:nrow(alldat[[num]]),nob,FALSE)#####sample 100 time steps
train=input2[s,,]
train=array(train,c(length(s)*npoints,nvar+5))
test=array(input2,c(hh*npoints,nvar+5))
###########emulation
Y=train[,1:5] #######emulation of average diameter
input=as.data.frame(train[,6:28])
names(input)=para
DM=input
DM_new=as.data.frame(input2[,6:28])
names(DM_new)=para
#########################fit independent kriging models
library(doParallel)
cl <- makeCluster(5)
registerDoParallel(cl)
myresult<-{foreach(i=1:5,.verbose=T,.errorhandling="stop") %dopar% {
library(DiceKriging)
#s2=sample(1:nrow(DM),5000,FALSE)
s2=sample(1:nrow(DM),2000,FALSE)###subsample training data
mod1=km(~.,design=DM[s2,],response=Y[s2,i],covtype="gauss",nugget.estim=TRUE,estim.method="MLE",optim.method="BFGS",control=list(po.size=50,trace=TRUE))
#mod1=km(~.,design=DM[s2,],response=Y[s2,i],covtype="matern5_2",nugget=0,nugget.estim=TRUE,estim.method="MLE",optim.method="BFGS",control=list(po.size=50,trace=TRUE))
pp=predict(mod1,newdata=data.frame(DM_new),se.compute=TRUE,type="UK")

lammp=test[,i]*10^6;emu=cbind(pp$mean,pp$lower95,pp$upper95)
prop=1-(sum((lammp-emu[,1])^2)/sum((lammp- mean(lammp))^2))
result=list(lammp,emu,prop,mod1,pp$upper95)
}}
stopCluster(cl)
save(myresult,file=paste(dirname,"myresult",sep="/"))

###############Plots
result=myresult[[1]]####for the 1st particle
#id=1:1250
#id=645:695
id=sample(1:nrow(DM_new),120,FALSE)
lammps=result[[1]][id]
emus=result[[2]][id,1]
U=result[[2]][id,3]
L=result[[2]][id,2]
ind=order(lammps)
x1 = min(L)
x1 = min(x1,lammps)
x2 = max(U)
x2 = max(x2, lammps)
rr=list(ind[1:30],ind[31:60],ind[61:90],ind[91:120])
for(i in 1:4){
#ind=order(lammps[[rr]])
lam=lammps[rr[[i]]]
em=emus[rr[[i]]]
UU=U[rr[[i]]]
LL=L[rr[[i]]]
pdf(paste("plot","i"))
par(mar=c(5.1,5.4,4.1,2.1))
plot(lam,lam,main=paste("Lammp Vs Emulator","i"),xlab=expression("Observed"~(10^{-6})),ylab=expression("Predicted"~(10^{-6})),lty=1,col="black",cex.lab=1.6,cex.axis=1.6,cex.main=1.6,xlim=c(x1,x2),ylim=c(x1,x2))
lines(lam,em, lty = 1,col="green")
lines(lam,UU, col="red",lty=1)
lines(lam,LL, col="red",lty=1)
pp$upper95=result[[5]]
legend(.6,max(pp$upper95),c("Predicted","Observed","95% C.I"),fill=c("green","black","red"),bty="0",border="black",cex=1.75,horiz=FALSE)
dev.off()
}
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

###########Sensitivity
library(sensitivity)

