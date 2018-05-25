#######
nvar=7
npoints=100######nunber of experimental design points
npoints2=5###number of stochastic repetitions
 #val=c("0.002","1.5e-3","1e-3","2.81e-4","0.005")#c(0.000035,0.00028,0.61,0.22)
  val=c("30e-3","10e-3","6e-3","0.00006944444","0.00001158333","0.63","0.24")
 vall=as.numeric(val)
 ppara=c("sub l","o2 l","nh4 l","het 0.00006944444","aob 0.00001158333","het 0.63","aob 0.24")
 #ppara=c("o2 l 2.81e-4 2.81e-4 2.81e-4 2.81e-4 2.81e-4 2.81e-4 2.81e-4")
 dirname=getwd()
 #set.seed(3)
 library(lhs)
 hyper=maximinLHS(npoints,nvar, 2)
 design= matrix(NA,nrow=npoints, ncol=nvar)
 for(j in 1:nvar){
 design[,j]=round(qunif(hyper[,j], min=0.5*vall[j], max=1.5*vall[j]),6)##+-50%
 }
 text <- readLines("atom.in",encoding="UTF-8")
 nrep=npoints
 save(design,file=paste(dirname,"design",sep="/"))

 for(rep in 1:nrep){
 for(j in 1:length(val)){
 text=gsub(val[j],design[rep,j],text)
 }
 writeLines(text,sub("rep",rep,paste("atom", "rep.in",sep="")))
 text <- readLines("atom.in",encoding="UTF-8")
 }
 #options(scipen=0)
 for( rep in 1: nrep){
 dir.create(sub("rep",rep,paste("input", "rep",sep="")))
 }

 for( rep in 1: nrep){
 file.copy(c(sub("rep",rep,paste("atom","rep.in",sep="")),"Inputscript.lammps"),to=sub("rep",rep,paste("input", "rep",sep="")))
 }

 for( rep in 1: nrep){
 setwd(paste(sub("rep",rep,"inputrep")))
 file.rename(from=sub("rep",rep,paste("atom", "rep.in",sep="")),to=sub("rep",rep,paste("atom", ".in",sep="")))
 setwd("..")
 }

 ###stochasticity
for(fold in 1:npoints){
 source("stoch.R",echo=TRUE)
 }
                     
