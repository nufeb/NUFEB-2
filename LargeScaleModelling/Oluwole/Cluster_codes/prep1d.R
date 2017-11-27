#####ve_concentration.csv  biomass.csv

path2=getwd()
CONC=CONC0=CONC1=list()
CON=CON0=CON1=list()
npoints2=5
npoints=100
for(j in 1:npoints){
setwd(sub("j",j,paste(getwd(),"inputj",sep="/")))
path=getwd()
for (num in 1:npoints2){
setwd(sub("num",num,paste(getwd(),"inputnum",sep="/")))
con=as.matrix(read.csv("Results/ave_concentration.csv",header=TRUE))[,1:5]
ll=nrow(con)
bio=as.matrix(read.csv("Results/biomass.csv",header=TRUE,nrows=ll))[,1:4]
height=as.matrix(read.csv("Results/ave_height.csv",header=FALSE,nrows=ll))[,2]
rough=as.matrix(read.csv("Results/roughness.csv",header=FALSE,nrows=ll))[,2]
div=as.matrix(read.csv("Results/diversity.csv",header=FALSE,nrows=ll))[,2]
types=as.matrix(read.csv("Results/ntypes.csv",header=TRUE,nrows=ll))[,1:4]
#a=rbind(con[1,],con)
#b=rbind(types[1,],types)
#biom=matrix(as.numeric(bio),dim(bio))[-nrow(bio),-ncol(bio)]
mass=rowSums(bio)##total biomass
res=cbind(con,mass,types,height,rough,div)
#res=cbind(con,mass)
CONC[[num]]=res
library(reshape2)
b1=do.call(rbind,lapply(CONC,melt))
b2=acast(b1,b1$Var1~b1$Var2,mean)
b3=acast(b1,b1$Var1~b1$Var2,sd)
CONC0=b2
CONC1=b3
setwd("../")
}
print(j)
print(num)
save(CONC,file=paste(path,"CONC",sep="/"))
save(CONC0,file=paste(path,"CONC0",sep="/"))
save(CONC1,file=paste(path,"CONC1",sep="/"))
###
CON[[j]]=CONC
CON1[[j]]=CONC1
CON[[j]]=CONC0
setwd("../")
}

save(CON,file=paste(path2,"CON",sep="/"))
save(CON0,file=paste(path2,"CON0",sep="/"))
save(CON1,file=paste(path2,"CON1",sep="/"))
#############END
