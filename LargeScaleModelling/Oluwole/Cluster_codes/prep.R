##################post processing
setwd(getwd())
library(iDynoR)
library(abind)
library(reshape2)
library(XML)
library(xml2)
dirname=getwd()
ALL0=ALL=list()
ALL1=ALL00=list()
lap=1:100
npoints2=5
T=49
A=list()
B=list()
res=list()
for(num in lap){
setwd(sub("num",num,"inputnum"))
obj=list.files(pattern=sub("num",num,"multi3D6newnum"))
for(i in 1:npoints2){
dirr=obj[i]
for(t in 0:(T-1)){
#r1=readSimResultFile(paste(dirr),"agent_State",2*t)
#r1b=mean(agent_returnSpeciesResultData(r1)[[1]][,13])##mean biofilm radius
r2=readSimResultFile(paste(dirr),"agent_Sum",2*t)
r2b=agent_returnSpeciesResultData(r2)[[1]]##total pop, mass
r3=readSimResultFile(paste(dirr),"env_Sum",2*t)
r3a=env_returnGlobalProductionRates(r3)
r3b=env_returnConcentrationAndRateChange(r3)##COD concentration
r4=readSimResultFile(paste(dirr),"env_State",2*t)
r4a=env_returnMeanBiofilmThickness(r4)
#class(r3b)=class(r3a)=NULL
class(r4a)=NULL
r3a=unlist(r3a)
r3b=unlist(r3b)
r3a=c(as.numeric(r3a[[2]]),as.numeric(r3a[[4]]),as.numeric(r3a[[8]]))#02,nh4,cod
r3b=c(as.numeric(r3b[[2]]),as.numeric(r3b[[4]]),as.numeric(r3b[[8]]))##02,nh4,cod
r4a=as.numeric(r4a[[6]])
r3a=signif(r3a,3)
r3b=signif(r3b,3)
r4a=signif(r4a,3)
#dat=c(r2b,r1b,r3a,r3b,r4a)
dat=c(r2b,r4a,r3b)##
A[[t+1]]=dat
}
nam=c("pop","mass","growthrate","meanthickness","o2_con","nh4_con","cod_con")
res=abind(A,along=0)
colnames(res)=nam
B[[i]]=res

library(reshape2)
b1=do.call(rbind,lapply(B,melt))
b2=acast(b1,b1$Var1~b1$Var2,mean)
b3=acast(b1,b1$Var1~b1$Var2,sd)
ALL0=b2
ALL1=b3
path=getwd()
save(ALL1,file=paste(path,"ALL1",sep="/"))
save(ALL0,file=paste(path,"ALL0",sep="/"))
}
ALL00[[num]]=ALL1##sd
ALL[[num]]=ALL0##mean
setwd("../")
}
path=getwd()
save(ALL,file=paste(path,"ALL",sep="/"))
save(ALL00,file=paste(path,"ALL00",sep="/"))
#############################################END
