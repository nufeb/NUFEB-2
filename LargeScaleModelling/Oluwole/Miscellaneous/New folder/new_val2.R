###########################validation scripts
set.seed(1)
m=list()
dirname=getwd()
load("result11/mod")
load("result11/scale");load("result11/center")#;load("pow")
load("ALL2")
load("alpha3")
load("design")
imp2=c(7:11,13)##c(nHET,nAOB,nNOB,nEPs,nDEAD,nparticle)
imp1=c(15,3,12,13,1)#
library(abind)
library(MuFiCokriging)
first=list();id=list()
for(i in 1:length(alpha3)){
first[[i]]=ALL2[[i]][1,imp1,1]
id[[i]]=nrow(alpha3[[i]])
}
idd=abind(id)
ind=idd[281:300]
f=abind(first,along=0)
for(k in 1:20){
#k=9
#imp2=c(7:11,13)##c(nHET,nAOB,nNOB,nEPs,nDEAD,nparticle)
#imp1=c(15,3,12,13,1)#
ndes2=280#274##283
nvar=27+5
nout=5
ntime=(ind[k])#72
ALL3=ALL2[(ndes2+k):(ndes2+k)]
pow=c(1,1e+06,1,1,1)
para=c("KsHET","Ko2HET","Kno2HET","Kno3HET","Knh4AOB","Ko2AOB","Kno2NOB","Ko2NOB","MumHET","MumAOB","MumNOB","etaHET","bHET","bAOB","bNOB","bEPS","YHET","YAOB","YNOB","YEPS","Y1","Do2","Dnh4","Dno2","Dno3","Ds","diffT")
nam3=c("y1","y2","y3","y4","y5")
nam2=c("s","o2","no2","no3","nh4")
des=(matrix(design[(ndes2+k):(ndes2+k),],nrow=length(ALL3)))
nut=alpha3[(ndes2+k):(ndes2+k)]
nut=abind(nut,along=0)
nut=nut[,,-1]
des=matrix(rep(des,ntime),nrow=ntime,byrow=TRUE)
des=cbind(des,nut)
inp=scale(des,scale=scale,center=center)
alpha1=list();alpha2=list()
for(i in 1:length(ALL3)){
alpha1[[i]]=ALL3[[i]][1,,1]
}
y0=abind(alpha1,along=0)
y0=matrix(rep(y0[,imp1],ntime),nrow=ntime,byrow=TRUE)

j=1
ini=matrix(apply(f,2,min),ncol=nout)#y0
ini[,c(4,5)]=log(ini[,c(4,5)])##
ini[,2]=ini[,2]*pow[2]
input=t(matrix(abind(inp[j,],ini[j,])))
colnames(input)=c(para,nam2,nam3)

##iterative use of emulator (dynamically)
##j==time points
nsim=100
result=list();mstar2=list()
m1=predict(mod,newdata=input,type="UK",se.compute=TRUE)
mu=abind(m1$mux,along=2)
varr=abind(m1$varx,along=2)#variance
V1=varr#/pow
V2=0
npoints=1
library(MASS)
mstar=mvrnorm(nsim,c(mu),diag(c(V1+V2),ncol=nout,nrow=nout))
mstar2[[j]]=mstar
result[[j]]=cbind(mu,V1+V2)
#
for(j in 2:ntime){
g1=matrix(rep(t(inp[j,]),nsim),nrow=nsim,ncol=nvar,byrow=TRUE)
g2a=mstar2[[j-1]]
input=cbind(g1,g2a)
colnames(input)=c(para,nam2,nam3)
m1=predict(mod,newdata=input,type="UK",se.compute=TRUE)
p1=array(abind(m1$mux,along=2),c(npoints,nsim,nout))
p2=array(abind(m1$varx,along=2),c(npoints,nsim,nout))
mu=t(matrix(apply(p1,3,rowMeans)))
mu2=array(rep(mu,nsim),c(npoints,nout,nsim));mu2=aperm(mu2,c(1,3,2))
varr=apply(p2,c(1,3),mean)#expectation of variance
varr2=apply((p1-mu2)^2,3,rowMeans)#variance of expectation
V1=varr#/pow
V2=varr2/pow
mstar=mvrnorm(nsim,c(mu),diag(c(V1+V2),ncol=nout,nrow=nout))
mstar2[[j]]=mstar
result[[j]]=cbind(mu,V1+V2)
}
m[[k]]=result
}
save(m,file=paste(dirname,"m",sep="/"))###spa
k=1
#res=abind(result,along=1)
res=abind(m[[k]],along=1)
emu=res[,1:5]
vv=res[,6:10]
ALL3=ALL2[(ndes2+1):300]
lam=ALL3[[k]][2:(ind[k]+1),imp1,1]
lam2=ALL3[[k]][2:(ind[k]+1),imp1,2]
l=list();l2=list()
for(i in 1:nout){
l[[i]]=lam[,i]*pow[i]
}
lammp=abind(l,along=2)
lammp[,c(4,5)]=log(lammp[,c(4,5)])
P=rep(NA,nout)
for(i in 1:nout){
P[i]=1-(sum((lammp[,i]-emu[,i])^2)/sum((lammp[,i]- mean(lammp[,1]))^2))
}

#[1] -0.5777962  1.0000000  0.9956696  0.9588905  0.9844120
#[1] -0.7139528  0.3814705  0.9968068  0.9876416  0.9901352
# -0.8094391  0.4479065  0.9978268  0.9886532  0.9902947
#ALL: nondynamic#0.9828075 0.9983645 0.9986431 0.9965313 0.9999741
com=lam=ALL3[[k]][2:nrow(ALL3[[k]]),4:6,1]



#############################
pow=c(1,1,1,1,1)
ntime2=72
dev=sqrt(vv)
emu2=res[1:ntime2,1:5]
tit=c("Species diversity",expression("Floc eqv. diameter"~(1e-06~m)),"Fractal dimension","No of particle",expression("Floc total mass"~(gram)))
dir.create("p1g")
dirn=paste(getwd(),"p1g",sep="/")
for(i in 1:(nout)){
dirname=file.path(dirn,paste("p1g_",i,".jpeg",sep=""))
jpeg(file=dirname,quality=100)
mytitle=tit[i]
lammps=lammp[,i]*pow[i]
emu=emu2[,i]*pow[i]
sd=dev[,i]*pow[i]
U=emu+(2*sd)
L=emu-(2*sd)
#time=(1:ntime)*10^4
time=(1:ntime2)*10^4/(24*60*60)###days
x1 = min(U,lammps,emu,L)
x2 = max(U,lammps,emu,L)
par(mar=c(5.1,5.4,4.1,2.1))
plot(time,lammps,main=mytitle,ylab=expression("Output"),lty=1,col="black",xlab="time (day)",cex.lab=1.7,cex.axis=1.7,cex.main=1.7,xlim=range(time),ylim=c(x1,x2)+c(0,.03))
lines(time,L, col="red",lwd=3,lty=2)
lines(time,U, col="red",lwd=3,lty=2)
lines(time,emu,col="green",lwd=3)
legend(1,(x2+.03),c("Predicted","Observed","95% C.I"),fill=c("green","black","red"),bty="0",border="black",cex=1.75,horiz=FALSE)
dev.off()
}
#[1] 0.7253680 1.0000000 0.9943814 0.9864976 0.9999772
##########################################
dir.create("p1a1")
dirn=paste(getwd(),"p1a1",sep="/")
g=list()
g=list(c(-.1,.1),c(-15,15),c(-.1,.1),c(-3,3),c(-5,3))
pos=c("topright","topright","topright","topleft","topleft")
for(i in 1:(nout)){
#i=1
dirname=file.path(dirn,paste("p1a1_",i,".jpeg",sep=""))
jpeg(file=dirname,quality=100)

lammps=lammp[,i]*pow[i]
emu=emu2[,i]*pow[i]
sd=dev[,i]*pow[i]
U=emu+(2*sd)
L=emu-(2*sd)
x1 = min(min(density(U)$y),min(density(lammps)$y),min(density(emu)$y),min(density(L)$y))
x2 = max(max(density(U)$y),max(density(lammps)$y),max(density(emu)$y),max(density(L)$y))
y1 = min(U,lammps,emu,L)
y2 = max(U,lammps,emu,L)
plot(density(lammps),main=tit[i],xlab="",cex.lab=1.5,cex.axis=1.6,cex.main=1.7,ylim=c(x1,x2)+c(0,.006),xlim=c(y1,y2)+g[[i]])
lines(density(emu),col="green")
lines(density(L), col="red",lwd=3,lty=2)
lines(density(U), col="red",lwd=3,lty=2)
legend(pos[i],c("Predicted","Observed","95% C.I"),
fill=c("green","black","red"),bty="0",border="black",cex=1.75,horiz=FALSE)
dev.off()
}
###############
library(rgl)
library(plot3D)
com1=com*1e05
scatter3D(x=com1[,1],y=com1[,2],z=com1[,3],colvar=emu2[,1],
pch = 16, cex = 1.5, clab = c("Fractal","dimension"),
main = "Floc morphology", ticktype = "detailed", theta = 10, d = 2,
colkey = list(length = 0.5, width = 0.5, cex.clab = 0.75

x=com1[,1];y=com1[,2];z=com1[,3]
rgl_init()
plot3d(x, y, z, col="blue", type ="p")
############

########################sensitivity
load("C:/Users/olu/Desktop/paper_analysis/new_floc2/result22/sens")
nam=c(para,nam2)
library(sensitivity)
dat=rbind(sens[[1]]$S$mean,sens[[2]]$S$mean,sens[[3]]$S$mean,sens[[4]]$S$mean,sens[[5]]$S$mean)
dat1=dat[,1:32]
######biofilm
load("C:\\Users\\olu\Desktop\\paper_analysis\\biofilm2\\result22/sens")
dat=rbind(sens[[1]]$S$mean,sens[[2]]$S$mean,sens[[3]]$S$mean)
dat2=dat[,1:32]
dat=rbind(dat1,dat2)
par(mar=c(5.1,4.5,1.9,.5)+0.0)
nam4=c("Floc total mass","Floc equivalent diameter","Floc fractal dimension","Floc total number of particle","Floc species diversity","Biofilm average height","Biofilm surface roughness","Biofilm segregation index")
xvals=barplot((dat),ylab="Sensitivity indices",main="Sensitivity analysis",offset=0,cex.lab=1.6,col=c("green","blue","pink","red","grey","yellow","blueviolet","brown"),density=c(40,60,90,150,180,235,270,300),cex.main=1.8,space=0.3,angle=c(45,90,135,180,225),cex.axis=1.6,cex.names=1.6,args.legend=list(x=c(3,8),y=c(0.9,2),bty="n",cex=1.6),legend=paste(nam4))

text(xvals,par("usr")[3]-.01,srt=45,adj=1,labels=nam,xpd=TRUE,cex=1.6)
###
ll=colMeans(emu2-2*sqrt(vv))
uu=colMeans(emu2+2*sqrt(vv))
barplot2(colMeans(emu),col=c("green","blue","pink","red","grey"),plot.ci=TRUE,ci.u=uu,ci.l=ll)
########################mass histogram

(hist(exp(bb3[,1])),main="Histogram",xlab="Floc total mass",cex.main=1.5,cex.lab=1.5,cex.axis=1.5,col="pink")
