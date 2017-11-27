#######Dicekriging emulator
load("alldat")
load("design")
npart=4###(diameter,x, y z)
natom=5###AOB,NOB,EPS, HET, inert
nvar=22
nob=40###100
npoints=jj=50##1000
runtime=86400###352000
step=2000
para=c("KsHET","Ko2HET","Kno2HET","Kno3HET","Knh4AOB","Ko2AOB","Kno2NOB","Ko2NOB","MumHET","MumAOB","MumNOB","etaHET","bHET","bAOB","bNOB","bEPS","YEPS","YHET","EPSdens","EPSratio","factor","ke","time")

val=c(0.01,0.81,0.0003,0.0003,0.001,0.0005,0.0013,0.00068,0.00006944444,0.00003472222,0.00003472222,0.6,0.00000462962,0.00000127314,0.00000127314,0.00000196759, 0.18,0.61,30,1.25,1.5,5e+10)

set.seed(3)
library(DiceKriging)
library(abind)
hh=nrow(alldat[[1]])#####be careful
###############input
#library(lhs)
#hyper=maximinLHS(npoints,nvar, dup=2)
#design= matrix(NA,nrow=npoints, ncol=nvar)
f#or(j in 1:(nvar-1)){
#design[,j]= qunif(hyper[,j], min=0.8*val[j], max=1.2*val[j])
#}
##########################################
m1=array(rep(design,hh),c(dim(design),hh))
m2=aperm(m1,c(3,1,2))
test=abind(alldat,along=1)
test2=array(test,c(hh,jj,natom,npart))
test=test2[,,,1]###choose diameter/ for all 5 particles
input2=abind(test,m2,along=3)
time=cap=seq(0,runtime,step)
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
DM_new=as.data.frame(test[,6:28]);test2=test[,1:5]
names(DM_new)=para
#########################fit independent kriging models
library(doParallel)
cl <- makeCluster(natom)
registerDoParallel(cl)
myresult<-{foreach(i=1:natom,.verbose=T,.errorhandling="stop") %dopar% {
library(DiceKriging)
#s2=sample(1:nrow(DM),5000,FALSE)
s2=sample(1:nrow(DM),1000,FALSE)###subsample training data
mod1=km(~.,design=DM[s2,],response=Y[s2,i],covtype="gauss",nugget.estim=TRUE,estim.method="MLE",optim.method="BFGS",control=list(po.size=50,trace=TRUE))
#mod1=km(~.,design=DM[s2,],response=Y[s2,i],covtype="matern5_2",nugget=0,nugget.estim=TRUE,estim.method="MLE",optim.method="BFGS",control=list(po.size=50,trace=TRUE))
pp=predict(mod1,newdata=data.frame(DM_new),se.compute=TRUE,type="UK")

lammp=test2[,i];emu=cbind(pp$mean,pp$lower95,pp$upper95)
prop=1-(sum((lammp-emu[,1])^2)/sum((lammp- mean(lammp))^2))
result=list(lammp,emu,prop,mod1,pp$upper95)
}}
stopCluster(cl)
save(myresult,file=paste(dirname,"myresult",sep="/"))

###############Plot
power=10^6
particles=c("HET","AOB","NOB","EPS","INERT")
for(p in 1:natom){
#p=3###NOB
result=myresult[[p]]####for the 1st particle
res1=array(result[[1]],c(hh,npoints,1))###lammp
res2=array(result[[2]],c(hh,npoints,3))##emulator

id=sample(1:nrow(DM_new),120,FALSE)
lammps=result[[1]][id]*power
emus=result[[2]][id,1]*power
U=result[[2]][id,3]*power
L=result[[2]][id,2]*power
ind=order(lammps)
rr=list(ind[1:30],ind[31:60],ind[61:90],ind[91:120])

for(i in 1:4){
dirname=file.path(getwd(),paste("plot_",particles[p],i,".pdf",sep=""))
pdf(file=dirname)
mytitle=paste("Mean diameter for Lammp Vs Emulator",particles[p])
lam=lammps[rr[[i]]]
em=emus[rr[[i]]]
UU=U[rr[[i]]]
LL=L[rr[[i]]]
x1 = min(UU,LL,lam,em)
x2 = max(UU,LL,lam,em)
par(mar=c(5.1,5.4,4.1,2.1))
plot(lam,lam,main=mytitle,xlab=expression("Observed"~(10^{-6})),ylab=expression("Predicted"~(10^{-6})),lty=1,col="black",cex.lab=1.6,cex.axis=1.6,cex.main=1.6,xlim=c(x1,x2),ylim=c(x1,x2))
lines(lam,em, lty = 1,col="green")
lines(lam,UU, col="red",lty=1)
lines(lam,LL, col="red",lty=1)
legend(x1,x2,c("Predicted","Observed","95% C.I"),fill=c("green","black","red"),bty="0",border="black",cex=1.75,horiz=FALSE)
dev.off()
}
#############plot as a function of time
#h=1:npoints

for(h in 14:18){
dirname=file.path(getwd(),paste("tplot_",particles[p],h,".pdf",sep=""))
pdf(file=dirname)
mytitle=paste("Mean diameter for Lammp Vs Emulator",particles[p])
l=res1[,h,1]*power
e=res2[,h,1]*power
uu=res2[,h,2]*power
ll=res2[,h,3]*power

x1 = min(min(l),lammps)
x2 = max(max(uu),lammps)
par(mar=c(5.1,5.4,4.1,2.1))
plot(time,l,main=mytitle,xlab=expression("Time"~(seconds)),ylab=expression("Predicted"~(10^{-6})),lty=1,col="black",cex.lab=1.6,cex.axis=1.6,cex.main=1.6,xlim=range(time),ylim=c(x1,x2)+c(0,.03))
for(i in 1:length(time)){
lines(rep(time[i],2),c(ll[i],uu[i]), col="red")
}
points(time,e,col="green")
legend(x1,x2+.04,c("Predicted","Observed","95% C.I"),fill=c("green","black","red"),bty="0",border="black",cex=1.75,horiz=FALSE)
dev.off()
}
}
###########Sensitivity
library(sensitivity)

