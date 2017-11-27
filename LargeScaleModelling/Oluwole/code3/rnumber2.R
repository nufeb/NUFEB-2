############
#<param name="Sbulk" unit="g.L-1">10e-3</param>
#<param name="Sin" unit="g.L-1">10e-3</param>
#<param name="Sbulk" unit="g.L-1">6e-3</param>
#<param name="Sin" unit="g.L-1">6e-3</param>
#<param name="Sbulk" unit="g.L-1">10e-3</param>
#<param name="Sin" unit="g.L-1">10e-3</param>

ppara=c("\t\t\t\t<param name=\"Sbulk\" unit=\"g.L-1\">10e-3</param>","\t\t\t\t<param name=\"Sin\" unit=\"g.L-1\">10e-3</param>","\t\t\t\t<param name=\"Sbulk\" unit=\"g.L-1\">6e-3</param>", 
"\t\t\t\t<param name=\"Sin\" unit=\"g.L-1\">6e-3</param>","\t\t\t\t<param name=\"Sbulk\" unit=\"g.L-1\">10e-3</param>",
"\t\t\t\t<param name=\"Sin\" unit=\"g.L-1\">10e-3</param>") 

nvar=3#*2
npoints=100######nunber of experimental design points
npoints2=5###number of stochastic repetitions
val=c("10e-3","6e-3","10e-3")
val=rep(val,each=2)
vall=as.numeric(val)
dirname=getwd()
library(lhs)
hyper=maximinLHS(npoints,nvar, 2)
design= matrix(NA,nrow=npoints, ncol=2*nvar)
for(j in 1:nvar){
design[,j]=round(qunif(hyper[,j], min=0.5*vall[j], max=1.5*vall[j]),6)##+-50%
}
#text <- readLines("atom.in",encoding="UTF-8")
text=readLines("multi_species_multi3D.xml",encoding="UTF-8")
nrep=npoints
des=matrix(rep(design,each=2),nrow=npoints,ncol=2*nvar,byrow=TRUE)
design=des
save(design,file=paste(dirname,"design",sep="/"))

ind=list()
for(j in 1:length(val)){
ind[[j]]=(grep(ppara[j],text))#[[1]]
}
ind=unique(sort(unlist(ind)))
b=list()
for(rep in 1:nrep){
for(j in 1:length(val)){
text[ind][j]=gsub(val[j],design[rep,j],text[ind][j])
}
writeLines(text,sub("rep",rep,paste("multiEPS3D","rep.xml",sep="")))
text=readLines("multi_species_multi3D.xml",encoding="UTF-8")
}

for(rep in 1:nrep){
dir.create(sub("rep",rep,paste("input", "rep",sep="")))
}

for(rep in 1:nrep){
file.copy(c(sub("rep",rep,paste("multiEPS3D","rep.xml",sep=""))),to=sub("rep",rep,paste("input", "rep",sep="")))
}



