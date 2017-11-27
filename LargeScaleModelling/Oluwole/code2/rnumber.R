#####################Generate copies of lammp input script and place them
####in diferent folders
nvar=20
npoints=100######nunber of experimental design points
npoints2=7###number of stochastic repetitions
##we have remove ke
para=c("KsHET","Ko2HET","Kno2HET","Kno3HET","Knh4AOB","Ko2AOB","Kno2NOB","Ko2NOB","MumHET","MumAOB","MumNOB","etaHET","bHET","bAOB","bNOB","bEPS","YEPS","YHET","EPSratio","factor")
#val=c(0.01,0.81,0.0003,0.0003,0.001,0.0005,0.0013,0.00068,0.00006944444,0.00003472222,0.00003472222,0.6,0.00000462962,0.00000127314,0.00000127314,0.00000196759, 0.18,0.61,1.25,1.5,5e+10)

val=c(0.00967788609690964,0.888134026480615,0.0003,0.0003,0.00108113164926693,0.0005,0.00115419633757137,0.000674405442599356,
0.00006944444,0.00003472222,0.00003472222,0.613162181432918,0.00000462962,0.00000127314,0.00000127314,0.00000196759,0.167053869057409,0.692894242844563,1.25,1.5)

ppara=c("variable KsHET equal 0.00967788609690964","variable Ko2HET equal 0.888134026480615","variable Kno2HET equal 0.0003","variable Kno3HET equal 0.0003","variable Knh4AOB equal 0.00108113164926693","variable Ko2AOB equal 0.0005","variable Kno2NOB equal 0.00115419633757137","variable Ko2NOB equal 0.000674405442599356","variable MumHET equal 0.00006944444","variable MumAOB equal 0.00003472222","variable MumNOB equal 0.00003472222","variable etaHET equal 0.613162181432918","variable bHET equal 0.00000462962","variable bAOB equal 0.00000127314","variable bNOB equal 0.00000127314","variable bEPS equal 0.00000196759","variable YEPS equal 0.167053869057409","variable YHET equal 0.692894242844563","variable EPSratio equal 1.25","variable factor equal 1.5")


dirname=getwd()
set.seed(3)
library(lhs)
hyper=maximinLHS(npoints,nvar, 2)
design= matrix(NA,nrow=npoints, ncol=nvar)
for(j in 1:nvar){
design[,j]= qunif(hyper[,j], min=0.8*val[j], max=1.2*val[j])
}
text <- readLines("Inputscript.lammps",encoding="UTF-8")
nrep=npoints
save(design,file=paste(dirname,"design",sep="/"))
design=rbind(val,design)
#########################################index points
id=c(rep(NA,length(ppara)))
for(i in 1:length(ppara)){
id[i]=grep(ppara[i],text,value=FALSE,fixed=TRUE)
}
#id=c(49:52,54:55,57:58,62:65,68,70,72,74,76:79,81,83)
for(rep in 1:nrep){
for(j in 1:nvar){
text[id[j]]=gsub(design[rep,j] ,design[rep+1,j], text[id[j]])
#text[id[j]]=gsub(val[j] ,design[rep,j], text[id[j]])
#writeLines(text,sub("rep",rep,paste("Input", "rep.lammps",sep="")))
}
writeLines(text,sub("rep",rep,paste("Inputscript", "rep.lammps",sep="")))
#source("stoch1.R",echo=TRUE)
}


for( rep in 1: nrep){
dir.create(sub("rep",rep,paste("input", "rep",sep="")))
}

for( rep in 1: nrep){
file.copy(c(sub("rep",rep,paste("Inputscript", "rep.lammps",sep="")),"flat_surface.in"), sub("rep",rep,paste("input", "rep",sep="")))
source("stoch1.R",echo=TRUE)
}
for( rep in 1: nrep){
file.remove(sub("rep",rep,paste("Inputscript", "rep.lammps",sep="")))
}
########initial conditions
#for(fold in 1:npoints){
#source("stoch1.R",echo=TRUE)
#}

###stochasticity
for(fold in 1:npoints){
source("stoch.R",echo=TRUE)
}
############################################################END
rm(list=ls())
                                                 

