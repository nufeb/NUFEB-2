######################
v=readLines("new.txt")
v2=c(unlist(strsplit(v, "input", fixed =FALSE)))
v3=strsplit(v2,"/res")
index=list()
for (n in 1:(length(v3)/2)){
index[[n]]=as.numeric(v3[[(2*n)]])
}
path=getwd()
lap=unlist(index)
ALL=list()
for(num in lap){
#for (num in 1:npoints){
setwd(sub("num",num,paste(getwd(),"inputnum",sep="/")))
load("res")
k=which(lap==num)
ALL[[k]]=res
setwd("../")
}
save(ALL,file=paste(path,"ALL",sep="/"))
