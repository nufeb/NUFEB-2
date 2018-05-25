################
################code for pre-processing LAMMPS data
#npoints=2
#for(xp in 1:npoints){
#setwd(sub("xp",xp,paste(getwd(),"inputxp",sep="/")))

#ncomp=70
ALL=list()
npoints2=10
path=getwd()
#library(doMPI)
#cl=startMPIcluster()
#registerDoMPI(cl)
#library(Rmpi)
cnam=c("type","diameter","mass","x","y","z","sub","o2","no2","no3","nh4")
sele=c("diameter","mass","x","y","z","sub","o2","no2","no3","nh4")
for (num in 1:npoints2){
#ALL=foreach(num=1:npoints2,.verbose=T,.errorhandling="stop") %:% 
setwd(sub("num",num,paste(getwd(),"inputnum",sep="/")))
#con<-file('snapshot.bubblemd')
#open(con)
#con2=readLines(con)
#index=grep("ITEM: ATOMS id type diameter mass x y z sub o2 no2 no3 nh4",con2,value=FALSE,fixed=TRUE)
con=readLines("snap2.txt")
con2=readLines("tot.txt")
ind2=as.numeric(unlist(strsplit(con2, "snapshot", fixed = TRUE)) [1])
ind=length(con)
h=unlist(strsplit(con, ":", fixed = TRUE)) 
index=rep(NA,ind)
for (n in 1:ind){
index[n]=as.numeric(h[(3*n)-2])
}
dd=list()
index=c(index,ind2+9)
library(doParallel)
cl=makeCluster(10)
registerDoParallel(cl)
#cl=startMPIcluster()
#registerDoMPI(cl)
library(data.table)
library(readr)
#dat=read_table("snapshot.bubblemd",col_names=TRUE,skip=index[1]-1)
#dat=as.matrix(dat)
#dd=foreach(i=1:(length(index)-1),.verbose=T,.errorhandling="stop") %dopar% {
dd=foreach(i=1:(length(index)-1),.verbose=FALSE,.errorhandling="stop") %do% {
library(data.table)
library(readr)
dat=fread("snapshot.bubblemd",sep2=" ",header=FALSE,skip=index[i]-1,nrows=index[i+1]-index[i]-9)
#dat1=dat[(index[i]-8):(index[i+1]-18),]
#dat1=matrix(as.numeric(unlist(tstrsplit(dat1,split=" ",fixed=TRUE))),ncol=12)
dat=as.matrix(dat)
dat=dat[,-1]
colnames(dat)=cnam
result=dat
}
stopCluster(cl)
rm(dat)
#mpi.quit()
dirname2=getwd()
save(dd,file=paste(dirname2,"dd",sep="/"))

####################################
ncomp2=15#(length(index)-15)
ncomp=length(dd)
comp=matrix(NA,nrow=ncomp,ncol=6)
mas=rep(NA,ncomp)
D=rep(NA,ncomp)
seg=rep(NA,ncomp)
A=rep(NA,ncomp);V=rep(NA,ncomp)
d1=matrix(NA,nrow=ncomp,ncol=4);nutrient=matrix(NA,nrow=ncomp,ncol=5)
Frac=rep(NA,ncomp);nparticle=rep(NA,ncomp)
for(j in 1:ncomp){
diam=dd[[j]][,2]##diameter
A[j]=sum(pi*diam^2)##particle surface area ###sum over j (each time step)
V[j]=sum((pi*diam^3)/6)
}
#BD=sqrt(A/pi)
BD2=(6*V/pi)^(1/3) ######Floc diameter using total volume
##############################
######Floc eqv diameter using distance technique
for(j in 1:ncomp){
dat=as.matrix(dd[[j]][,1:6])
com=colSums(dd[[j]][,3]*dd[[j]][,4:6])/sum(dd[[j]][,3])##com of floc
kala=rep(NA,nrow(dat))
for(k in 1:nrow(dat)){
kala[k]=sqrt(sum((com-dat[k,4:6])^2))###square distance
}
id=which(kala==max(kala))
diam=2*(max(kala)+(dat[id,2]/2))####this is diameter
coms=c(diam,com)
d1[j,]=coms
}
#########fractal calculation
for(j in 1:ncomp){
Ra=sqrt(sum(dd[[j]][,3]*(dd[[j]][,2])^2)/sum(dd[[j]][,3]))
ll=nrow(dd[[j]])
Rm=sum(dd[[j]][,2]/2)/ll
Frac[j]=log(Ra/Rm)/log(ll)
}
#vv=parLapply(cl,dd,myfunc2)
#########################Total Number of particle
for(j in 1:ncomp){
nparticle[j]=nrow(dd[[j]])####natom
}
##########################Floc composition
for(j in 1:ncomp){
for (type in c(1:4,6)){
gdat=dd[[j]]
comp[j,type]=length(which(gdat[,1]==type))
}}
comp=comp[,-5]
#############Simpson diversity index
for (j in 1:ncomp){
N=nparticle[j]
D[j]=1-sum((comp[j,])*(comp[j,]-1))/(N*(N-1))
}
############################nutrients
for(j in 1:ncomp){
dat2=as.matrix(dd[[j]][,7:11])
nutrient[j,]=colSums(dat2)
}
##############################Floc total mass
for(j in 1:ncomp){   

mas[j]=sum(dd[[j]][,3])

}

###########segregation index
#for(j in 1:ncomp){
for(j in 1:ncomp2){###just 40 timestep
dat=as.matrix(dd[[j]][,1:6])
S=rep(NA,nrow(dat))
kala=rep(NA,nrow(dat))
for(i in 1:nrow(dat)){
for (k in 1:nrow(dat)){
kala[k]=sqrt(sum((dat[i,4:6]-dat[(k),4:6])^2))###square distance
}
kala1=kala[-i]
#
id=which(kala1<(10*dat[i,2]))## 10 particle length
dat2=dat[id,]
rho=rep(NA,length(id))
for(p in 1:nrow(dat2)){####test p(ci,cj)==p(ci)
if(dat2[p,1]==dat[i,1]){
rho[p]=1
}else{
rho[p]=0
}
}
S[i]=mean(rho)
}
seg[j]=mean(S)
}
############surface roughness

#################

res=cbind(mas,BD2,d1,comp,Frac,nparticle,seg,D,nutrient)
ALL[[num]]=res
setwd("../")
}
save(ALL,file=paste(path,"ALL",sep="/"))
#closeCluster(cl)
#mpi.quit()
#setwd("../")
#}
#####################################################END
