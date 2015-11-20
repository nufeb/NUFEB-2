####################Main R scripts
#####sed -i -e 's/\r$//' olu4.sh  for running the main script
npoints=245###1000##npoints is the number of design points
npoints2=50####50#NOTE: number of stochasticity
npart=10###("diameter","mass","x","y","z","sub","no2","no3","o2","nh4")
natom=5###AOB,NOB,EPS, HET, inert
##
source("rlammp3.R",echo=TRUE)###read design points
source("rlammp3b.R",echo=TRUE)##read stochasticity
#######Dicekriging emulator
#load("ALL")
library(abind)
load("design")
nvar=21
nob=70###100##represents number of time points to sample
#npoints=50##1000
ntime=jj=nrow(ALL[[1]])
runtime=172800###352000
step=2000
para=c("KsHET","Ko2HET","Kno2HET","Kno3HET","Knh4AOB","Ko2AOB","Kno2NOB","Ko2NOB","MumHET","MumAOB","MumNOB","etaHET","bHET","bAOB","bNOB","bEPS","YEPS","YHET","EPSratio","factor","ke")
set.seed(3)
library(DiceKriging)
#library(abind)
###############input
des=list()
for(i in 1:ntime){
#des[[i]]=design
des[[i]]=design[1:npoints,]
}
des2=abind(des,along=0)
des2=aperm(des2,c(2,1,3))
test2=abind(ALL,along=0)
test=test2[,,,1]###choose diameter/ for all 5 particles/output data
#test=test2[,,,2]###choose mass/ for all 5 particles/output data
other_input=test2[,,,6:10]###choose c("sub","no2","no3","o2","nh4")
other1=other_input[,,1,]###type1 for all 5 different inputs
other2=other_input[,,2,]###type1
other3=other_input[,,3,]###type1
other4=other_input[,,4,]###type1
other5=other_input[,,5,]###type1
input2=abind(test,des2,other1,other2,other3,other4,other5,along=3)
time=cap=seq(0,runtime,step)
cap=array(rep(cap,ntime),c(ntime,npoints,1))
cap=aperm(cap,c(2,1,3))
input22=array(NA,c(dim(input2)[1:2],dim(input2)[3]+1))
input22[,,dim(input22)[3]]=cap####time
input22[,,1:(dim(input22)[3]-1)]=input2
input2=input22

s=sample(1:ntime,nob,FALSE)#####sample 100 time steps
train=input2[,s,]##sample "nob" time points
train=array(train,c(npoints*length(s),dim(input2)[3]))
test=array(input2,c(npoints*ntime,dim(input2)[3]))
###########emulation
namm=c("sub","no2","no3","o2","nh4")
nam=c(para,namm,"time")
Y=train[,1:5]*10^16 #######emulation of average diameter
input=as.data.frame(train[,-c(1:5)])
inp1=input[,c(1:21,22:26,47)];inp2=input[,c(1:21,27:31,47)];inp3=input[,c(1:21,32:36,47)];inp4=input[,c(1:21,37:41,47)];inp5=input[,c(1:21,42:46,47)]
names(inp1)=names(inp2)=names(inp3)=names(inp4)=names(inp5)=nam
DM=list(inp1,inp2,inp3,inp4,inp5)
input2=as.data.frame(test[,-c(1:5)]);test2=test[,1:5]*10^16
inp1=input2[,c(1:21,22:26,47)];inp2=input2[,c(1:21,27:31,47)];inp3=input2[,c(1:21,32:36,47)];inp4=input2[,c(1:21,37:41,47)];inp5=input2[,c(1:21,42:46,47)]
names(inp1)=names(inp2)=names(inp3)=names(inp4)=names(inp5)=nam
DM_new=list(inp1,inp2,inp3,inp4,inp5)
#names(DM_new)=c(para,nam,"time")
##################################stepwise
load("polyset")
quadterms=polySet(length(nam),2,2,2)
quad1 <-paste("~", makeScope(quadterms,nam))
form1=~time
scope1 <- list(lower =form1,upper =quad1)

#varsel<- rxStepControl(method = "stepwise", scope = scope1,maxSteps =1000,k=log(500000))
#mod1c<- rxLinMod(formula=form1,data=da2[1:2200000,],variableSelection=varsel,reportProgress =3,dropMain = FALSE,coefLabelStyle = "R",verbose=1,blocksPerRead=10000)


#########################fit independent kriging models
library(doParallel)
cl <- makeCluster(natom)
registerDoParallel(cl)
myresult<-{foreach(i=1:natom,.verbose=T,.errorhandling="stop") %dopar% {
library(DiceKriging)
library(abind)
mod1=lm(Y[,i]~.,data=DM[[i]])
stepmod=step(mod1,direction="both",scope=scope1,trace=1,steps=300)
form=formula(noquote(paste("~",formula(stepmod)[3])))
#s2=sample(1:nrow(DM),5000,FALSE)
s2=sample(1:nrow(DM[[i]]),5000,FALSE)###subsample training data
#dd=DM[[i]][s2,]
olu=abind(DM,along=0)
olu2=colMeans(olu);names(olu2)=nam
mod1=km(form,design=olu2[s2,],response=Y[s2,i],covtype="gauss",nugget.estim=TRUE,estim.method="MLE",optim.method="BFGS",control=list(po.size=50,trace=TRUE))
#mod1=km(~.,design=DM[s2,],response=Y[s2,i],covtype="matern5_2",nugget=0,nugget.estim=TRUE,estim.method="MLE",optim.method="BFGS",control=list(po.size=50,trace=TRUE))
pp=predict(mod1,newdata=data.frame(DM_new[[i]]),se.compute=TRUE,type="UK")

lammp=test2[,i];emu=cbind(pp$mean,pp$lower95,pp$upper95)
prop=1-(sum((lammp-emu[,1])^2)/sum((lammp- mean(lammp))^2))
result=list(lammp,emu,prop,mod1,pp$upper95)
}}
stopCluster(cl)
save(myresult,file=paste(dirname,"myresult",sep="/"))

###############Plot
dir.create("myplot")
dirn=paste(getwd(),"myplot",sep="/")
power=10^6
particles=c("HET","AOB","NOB","EPS","INERT")
for(p in 1:natom){
#p=3###NOB
result=myresult[[p]]####for the 1st particle
res1=array(result[[1]],c(npoints,ntime,1))###lammp
res2=array(result[[2]],c(npoints,ntime,3))##emulator

id=sample(1:nrow(DM_new),120,FALSE)
lammps=result[[1]][id]*power
emus=result[[2]][id,1]*power
U=result[[2]][id,3]*power
L=result[[2]][id,2]*power
ind=order(lammps)
rr=list(ind[1:30],ind[31:60],ind[61:90],ind[91:120])

for(i in 1:4){
dirname=file.path(dirn,paste("plot_",particles[p],i,".pdf",sep=""))
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
dirname=file.path(dirn,paste("tplot_",particles[p],h,".pdf",sep=""))
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
}}

###########Sensitivity
library(sensitivity)

