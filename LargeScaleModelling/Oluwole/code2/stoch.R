####in diferent folders
#fold=8
#npoints2=6#NOTE: npoints2 is the number of stochastic repetitions
nvar2=4
setwd(paste(getwd(),"/",gsub("fold",fold,"inputfold"),sep=""))
#para2=c("fix dt1 HET death 1 v_bHET v_factor","fix dt2 AOB death 1 v_bAOB v_factor","fix dt3 NOB death 1 v_bNOB v_factor","fix d1 HET divide 1000 v_EPSdens 2.0","fix d2 AOB divide 1000 v_EPSdens 2.0","fix d3 NOB divide 1000 v_EPSdens 2.0","ix e1 HET eps_extract 1000 v_EPSratio v_EPSdens")
para2=c("fix d1 HET divide 1000 v_EPSdens 2.0","fix d2 AOB divide 1000 v_EPSdens 2.0","fix d3 NOB divide 1000 v_EPSdens 2.0","fix e1 HET eps_extract 1000 v_EPSratio v_EPSdens")
#val2=c(1952467,1234312,1325352,1242242,1242242,1242242,1242242)
val2=c(1201075,1139394,1183134,1271254)

set.seed(3)
nre2=npoints2
hyper2=maximinLHS(npoints2,nvar2, 2)
design2= matrix(NA,nrow=npoints2, ncol=nvar2)
for(j in 1:nvar2){
design2[,j]=round(qunif(hyper2[,j], min=0.8*val2[j], max=1.2*val2[j]))
}
text2 <- readLines( gsub("fold",fold,paste("Inputscriptfold", ".lammps",sep="")),encoding="UTF-8")

save(design2,file=paste(getwd(),"design2",sep="/"))
#########################################
design2=rbind(val2,design2)
index=c(rep(NA,length(para2)))
for(i in 1:length(para2)){
index[i]=grep(para2[i],text2,value=FALSE,fixed=TRUE)
}
#id=c(49:52,54:55,57:58,62:65,68,70,72,74,76:79,81,83)
for(re2 in 1:nre2){
for(j in 1:nvar2){
text2[index[j]]=gsub(design2[re2,j] ,design2[(re2)+1,j], text2[index[j]])
}
writeLines(text2,gsub("re2",re2,paste("Inputscriptre2", ".lammps",sep="")))
}
for(re2 in 1: nre2){
dir.create(sub("re2",re2,paste("input", "re2",sep="")))
}

for( re2 in 1: nre2){
file.copy(c(sub("re2",re2,paste("Inputscriptre2", ".lammps",sep="")),"flat_surface.in"), sub("re2",re2,paste("input", "re2",sep="")))
}
for(re2 in 1: nre2){
file.remove(sub("re2",re2,paste("Inputscriptre2", ".lammps",sep="")))
}
setwd("..")
#######################end

