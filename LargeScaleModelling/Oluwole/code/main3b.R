######main3b.R
ALL=list()
#setwd("..")
#setwd(sub("num",num,paste(getwd(),"inputnum",sep="/")))
for(num in 1:nexp){
setwd(sub("num",num,paste(getwd(),"inputnum",sep="/")))
load("alldat")
test=abind(alldat,along=1)
test2=array(test,c(hh,jj,natom,npart))
test3=apply(test2,c(3,4),rowMeans)######average out stochasticity
ALL[[num]]=test3
}
