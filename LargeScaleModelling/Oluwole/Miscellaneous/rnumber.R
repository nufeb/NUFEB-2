#####################Generate copies of lammp input script and place them
####in diferent folders
#options(scipen=10)
nvar=27
npoints=300######nunber of experimental design points
npoints2=10###number of stochastic repetitions
val=c(0.01,0.81,0.0003,0.0003,0.001,0.0005,0.0013,0.00068,0.000069,0.0000087962,0.000009375,
0.6,0.0000046296296,0.000001273148,0.000001273148,0.000001967592,0.61,0.33,0.083,0.18,0.4,0.000000002,0.0000000014,0.0000000012,0.0000000012,0.0000000005,0.0001)


para=c("KsHET","Ko2HET","Kno2HET","Kno3HET","Knh4AOB","Ko2AOB","Kno2NOB","Ko2NOB","MumHET","MumAOB","MumNOB","etaHET","bHET","bAOB","bNOB","bEPS","YHET","YAOB","YNOB","YEPS","Y1","Do2","Dnh4","Dno2","Dno3","Ds","diffT")

ppara=c("variable KsHET equal","variable Ko2HET equal","variable Kno2HET equal","variable Kno3HET equal","variable Knh4AOB equal","variable Ko2AOB equal","variable Kno2NOB equal","variable Ko2NOB equal","variable MumHET equal","variable MumAOB equal","variable MumNOB equal","variable etaHET equal","variable bHET equal","variable bAOB equal","variable bNOB equal","variable bEPS equal","variable YHET equal","variable YAOB equal","variable YNOB equal","variable YEPS equal","variable Y1 equal","variable Do2 equal","variable Dnh4 equal","variable Dno2 equal","variable Dno3 equal","variable Ds equal","variable diffT equal")
dirname=getwd()
set.seed(3)
library(lhs)
hyper=maximinLHS(npoints,nvar, 2)
design= matrix(NA,nrow=npoints, ncol=nvar)
for(j in 1:nvar){
design[,j]= qunif(hyper[,j], min=0*val[j], max=2*val[j])##+-100%
}
text <- readLines("Inputscript.lammps",encoding="UTF-8")
nrep=npoints
save(design,file=paste(dirname,"design",sep="/"))

#########################################index points
id=c(rep(NA,length(ppara)))
for(i in 1:length(ppara)){
id[i]=grep(ppara[i],text,value=FALSE,fixed=TRUE)
}
for(rep in 1:nrep){
for(j in 1:nvar){
text[id[j]]=noquote(paste(ppara[j],design[rep,j]))
}
writeLines(text,sub("rep",rep,paste("Inputscript", "rep.lammps",sep="")))
}

#options(scipen=0)
for( rep in 1: nrep){
dir.create(sub("rep",rep,paste("input", "rep",sep="")))
}

for( rep in 1: nrep){
file.copy(c(sub("rep",rep,paste("Inputscript", "rep.lammps",sep="")),c("initial_cells_floc.in","AUTO_RESTART_PLEASE")), sub("rep",rep,paste("input", "rep",sep="")))

}
for( rep in 1: nrep){
file.remove(sub("rep",rep,paste("Inputscript", "rep.lammps",sep="")))
}

###stochasticity
for(fold in 1:npoints){
source("stoch.R",echo=TRUE)
}
############################################################END
rm(list=ls())
                                                 

