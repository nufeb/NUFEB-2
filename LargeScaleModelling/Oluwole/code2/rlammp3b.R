######main3b.R
library(abind)
ALL=list()
for(num in 1:npoints){
setwd(sub("num",num,paste(getwd(),"inputnum",sep="/")))
load("alldat")
test=abind(alldat,along=0)
#test2=array(test,c(length(alldat),ntime,natom,npart))
t#est2=array(test,c(length(alldat),dim(alldat[[1]])))
test3=apply(test,c(3,4),colMeans)######average out stochasticity
ALL[[num]]=test3
setwd("..")
}
dirname=getwd()
save(ALL,file=paste(dirname,"ALL",sep="/"))
