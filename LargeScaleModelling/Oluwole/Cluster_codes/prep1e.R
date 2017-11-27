#####combine

path2=getwd()
ALL2=ALL20=ALL21=list()
npoints=100
for(j in 1:npoints){
setwd(sub("j",j,paste(getwd(),"inputj",sep="/")))
load("ALL");load("ALL0");load("ALL1")

ALL2[[j]]=ALL
ALL21[[j]]=ALL1
ALL20[[j]]=ALL0
setwd("../")
}

save(ALL2,file=paste(path2,"ALL2",sep="/"))
save(ALL20,file=paste(path2,"ALL20",sep="/"))
save(ALL21,file=paste(path2,"ALL21",sep="/"))
#############END
