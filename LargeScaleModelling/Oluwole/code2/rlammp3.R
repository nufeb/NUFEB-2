###########R script for reading and processing of Lammp snapshot.bubblemd 
#library(doParallel)
#cl <- makeCluster(npoints)
#registerDoParallel(cl)
#mydat<-{foreach(xp=1:npoints,.verbose=T,.errorhandling="stop") %do% {
#mydat=list()
for(xp in 91:npoints){
setwd(sub("xp",xp,paste(getwd(),"inputxp",sep="/")))
#################################start of main code for each folder
path=getwd()
#cnam=c("type","diameter","mass")
cnam=c("type","diameter","mass","x","y","z","sub","no2","no3","o2","nh4")
sele=c("diameter","mass","x","y","z","sub","no2","no3","o2","nh4")
alldat=list()
alldat2=list()
alldat3=list()
for (num in 1:npoints2){
setwd(sub("num",num,paste(getwd(),"inputnum",sep="/")))
#ITEM: ATOMS id type diameter mass x y z sub no2 no3 o2 nh4
#index=grep("ITEM: ATOMS id type diameter mass",readLines("snapshot.bubblemd"),value=FALSE,fixed=TRUE)
#index=grep("ITEM: ATOMS id type diameter x y z vx vy vz fx fy fz",readLines("snapshot.bubblemd"),value=FALSE,fixed=TRUE)
index=grep("ITEM: ATOMS id type diameter mass x y z sub no2 no3 o2 nh4",readLines("snapshot.bubblemd"),value=FALSE,fixed=TRUE)
ind=length(readLines("snapshot.bubblemd"))
#ind=countLines("snapshot.bubblemd")
dd=list()
index=c(index,ind+9)
for (i in 1:(length(index)-1)){
#for (i in 1:6){
dd[i]=list(read.csv("snapshot.bubblemd",sep="",skip=index[i]-1,header=TRUE,colClasses=numeric(),nrows=index[i+1]-index[i]-9)[,2:(length(cnam)+1)]);names(dd[[i]])=cnam
}
setwd(path)
###
ww=array(NA,c(length(dd),natom,npart))###5 particle types
ww2=array(NA,c(length(dd),natom,npart))###5 particle types
ww3=rep(NA,length(dd))
for(j in 1:length(dd)){
for(k in 1:natom){
#ww[j,k,]=rbind(colMeans(subset(dd[[j]],type==k,select=sele)))##type k
ww[j,k,]=rbind(colMeans(subset(dd[[j]],type==k,select=sele)))
ww2[j,k,]=rbind(apply(subset(dd[[j]],type==k,select=sele),2,var))##type k
ww3[j]=nrow(dd[[j]])####natom
alldat[[num]]=ww;alldat2[[num]]=ww2;alldat3[[num]]=ww3
}}
setwd(path)
}
dirname=path
save(alldat,file=paste(dirname,"alldat",sep="/"));save(alldat2,file=paste(dirname,"alldat2",sep="/"))
save(alldat3,file=paste(dirname,"alldat3",sep="/"))
setwd("..")
#############################end of main code for each folder
}
#}
#stopCluster(cl)

#######################################END

