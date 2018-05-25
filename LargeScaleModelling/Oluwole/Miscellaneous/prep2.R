####################Futher processing of LAMMPS data
gdat=list()
lap=1:300
#t=c(13,90)
#t=c(13,31, 45,47,59,74,78,80,81,90,97,98,99,101,103,104,107,108,109,111,126,128,137,142,143,149,153,155,157,167,168,174,178,179,182,183,201,211,213,215,216,222,227,235,237,244,245,250,254,260,266,267,268,276,277,279,285,290)
lap2=lap
#miss=c(2,6,9,11,20,25,36,37:39,42,48,50,56,59,68,92,95)
#lap=lap[-miss]
#t=c(100,102,105,106,10,110,112,113,114,117,118,119,120,121:124,127,130,131:136,139,141,144,146,148,14,150,152,156,158,159,15,160,161,163,165,166,169,16,170,171,172,173,175,177,17,180,184,188,190,191,193,195,197,198,19,1,200,203,204,207,209,20,210,214,217,218,219,21,220,221,224,225,228,22,230,231,232,233,234,238,239,240,241,246,247,249,24,251,252,253,255,256,257,258,261,262,263,264,265,26,270,272,274,275,27,280,282,284,286,287,288,289,28,291,292,295,296,297,299,29,32,33,35,36,37,3,40,41,43,44,49,4,50,51,52,58,5,60,61,62,65,67,69,71,73,75,76,77,79,7,82,83,84,85,89,8,91,93,94,95,96)
#lap=sort(t)
for(num in lap2){
#for (num in 1:npoints){
setwd(sub("num",num,paste(getwd(),"inputnum",sep="/")))
load("ALL")
gdat[[num]]=ALL
setwd("../")
}
library(snow)
library(doParallel)
cl=makeCluster(10)
registerDoParallel(cl)
ALL2={foreach(k=lap2,.verbose=FALSE,.errorhandling="stop") %dopar% {
library(abind)
#ga=abind(gdat[[k]],along=3)
#g1=apply(ga,2,rowMeans)###average out stochasticity
ga=gdat[[k]]
nam=paste("v",1:20,sep="")
library(reshape2)
for(i in 1:length(ga)){colnames(ga[[i]])=nam}
b1=do.call(rbind,lapply(ga,melt))##matrix of unequal dimension
b2=acast(b1,b1$Var1~b1$Var2,mean)
b3=acast(b1,b1$Var1~b1$Var2,var)
#g2=apply(ga,c(1,2),var)
g3=abind(b2,b3,along=3)
}}
stopCluster(cl)
dirname=getwd()
save(ALL2,file=paste(dirname,"ALL2",sep="/"))
################r-values
gdat2=list()
npoints=300
npoints2=10
for(num in 1:npoints){
setwd(sub("num",num,paste(getwd(),"inputnum",sep="/")))
library(doParallel)
cl=makeCluster(10)
registerDoParallel(cl)
ALL={foreach(fold=1:npoints2,.verbose=TRUE,.errorhandling="stop") %do% {
vv=read.table(gsub("fold",fold,"inputfold/r-values/r-values.txt"))
}}
stopCluster(cl)
gdat2[[num]]=ALL
setwd("../")
}
save(gdat2,file=paste(dirname,"gdat2",sep="/"))
vv=list();vv2=list()
for(i in 1:300){
v=gdat2[[i]]
for(j in 1:length(v)){
v2=as.numeric(unlist(strsplit(as.matrix(v[[j]][-1,]),split=",")))
vv[[j]]=matrix(v2,ncol=6)
b1=do.call(rbind,lapply(vv,melt))
b2=acast(b1,b1$Var1~b1$Var2,mean)
}
vv2[[i]]=b2
}
alpha3=vv2
save(alpha3,file=paste(dirname,"alpha3",sep="/"))
###################################################END
