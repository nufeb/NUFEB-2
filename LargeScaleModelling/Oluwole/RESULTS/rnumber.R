#####################Generate copies of lammp input script and place them
####in diferent folders
para=c("KsHET","Ko2HET","Kno2HET","Kno3HET","Knh4AOB","Ko2AOB","Kno2NOB","Ko2NOB","MumHET","MumAOB","MumNOB","etaHET","bHET","bAOB","bNOB","bEPS","YEPS","YHET","EPSdens","EPSratio","factor","ke")

val=c(0.01,0.81,0.0003,0.0003,0.001,0.0005,0.0013,0.00068,0.00006944444,0.00003472222,0.00003472222,0.6,0.00000462962,0.00000127314,0.00000127314,0.00000196759, 0.18,0.61,30,1.25,1.5,5e+10)

set.seed(3)
library(lhs)
nvar=22
npoints=50
hyper=maximinLHS(npoints,nvar, 2)
design= matrix(NA,nrow=npoints, ncol=nvar)
for(j in 1:nvar){
design[,j]= qunif(hyper[,j], min=0.8*val[j], max=1.2*val[j])
}
text <- readLines("Inputscript.lammps",encoding="UTF-8")
nrep=npoints
#########################################
design=rbind(val,design)
id=c(49:52,54:55,57:58,62:65,68,70,72,74,76:79,81,83)
for(rep in 1:nrep){
for(j in 1:nvar){
text[id[j]]=gsub(design[rep,j] ,design[rep+1,j], text[id[j]])
#text[id[j]]=gsub(val[j] ,design[rep,j], text[id[j]])
#writeLines(text,sub("rep",rep,paste("Input", "rep.lammps",sep="")))
}
writeLines(text,sub("rep",rep,paste("Inputscript", "rep.lammps",sep="")))
}


for( rep in 1: nrep){
dir.create(sub("rep",rep,paste("input", "rep",sep="")))
}

for( rep in 1: nrep){
file.copy(c(sub("rep",rep,paste("Input", "rep.lammps",sep="")),"flat_surface.in"), sub("rep",rep,paste("input", "rep",sep="")))
}
for( rep in 1: nrep){
file.remove(sub("rep",rep,paste("Input", "rep.lammps",sep="")))
}
############################################################END
