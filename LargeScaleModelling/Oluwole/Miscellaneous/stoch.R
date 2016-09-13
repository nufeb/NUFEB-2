#NOTE: npoints2 is the number of stochastic repetitions
nvar2=2
setwd(paste(getwd(),"/",gsub("fold",fold,"inputfold"),sep=""))
para2=c("fix d1 all divide 500 v_EPSdens 2.0","fix e1 HET eps_extract 500 v_EPSratio v_EPSdens")
val2=c(3123,1242242)
set.seed(3)
nre2=npoints2
hyper2=maximinLHS(npoints2,nvar2, 2)
design2= matrix(NA,nrow=npoints2, ncol=nvar2)
for(j in 1:nvar2){
design2[,j]=round(qunif(hyper2[,j], min=0.90*val2[j], max=1.1*val2[j]))
}
text2 <- readLines( gsub("fold",fold,paste("Inputscriptfold", ".lammps",sep="")),encoding="UTF-8")

save(design2,file=paste(getwd(),"design2",sep="/"))
#########################################
design2=rbind(val2,design2)
index=c(rep(NA,length(para2)))
for(i in 1:length(para2)){
index[i]=grep(para2[i],text2,value=FALSE,fixed=TRUE)
}
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
file.copy(c(sub("re2",re2,paste("Inputscriptre2", ".lammps",sep="")),"initial_cells_floc.in","AUTO_RESTART_PLEASE"), sub("re2",re2,paste("input", "re2",sep="")))
}
for(re2 in 1: nre2){
file.remove(sub("re2",re2,paste("Inputscriptre2", ".lammps",sep="")))
}
setwd("..")
#######################end

