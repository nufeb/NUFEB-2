#####################Emulator
#setwd("C:/Users/olu/Desktop/NUFEB_materials/Long-Run")
fam=c("id","type","Dia","X","Y","Z","Vx","Vy","Vz","Fx","Fy","Fz")#####output names
dd=list()
r=list()
#nob=100
nfac=12
ww=1:50
dir="/data/noo11/train/test2"
#mydat[[3]][[120]]
library(doParallel)
cl <- makeCluster(length(ww))
registerDoParallel(cl)
mydat<-{foreach(nam=1:length(ww),.verbose=T,.errorhandling="stop") %dopar% {
L=paste(dir,"inputnam",sep="/")
setwd(sub("nam",nam,L))
nob=length(list.files())-6
for (i in 1:nob){
dd[i]=list(read.table(sub("i",i,(paste("data","i.xlsx",sep=""))),sep=",",col.names=fam,header=FALSE))
r[i]=list(nrow(dd[[i]]))
result=list(dd,r)
}
result
}}
stopCluster(cl)
save(mydat,file=paste(dir,"mydata",sep="/"))
#mydat[[1]][[2]][[2]]##50*2*var
## (mydat[[50]][[1]][[1]][[13]])##50*2*var*12
##nrow(mydat[[35]][[1]][[93]][1:12])##important to get no of particles
(mydat[[i]][[1]][[j]][[]])
mydat1=array(NA,c(50,100,12))
for( i in 1:50){
for(j in 1:100){
#mydat1[,i,j]=colMeans(mydat[[i]][[1]][[j]][[3]])####diameter
mydat1[i,j,]=colMeans(mydat[[i]][[1]][[j]][1:12])
}}
##########compute average and variance at each time-step
dd2=matrix(NA,nrow=nob,ncol=nfac)
dd2_var=matrix(NA,nrow=nob,ncol=nfac)
for(i in 1:nob){
dd2[i,]=rbind(colMeans(dd[[i]])) ##column mean
dd2_var[i,]=rbind(apply(dd[[i]],2,var))###column variance
colnames(dd2)=fam
}
###########
#######Read in the input data;sequence of time-steps
input=seq(0,352000,2000)

set.seed(10)
s=sample(1:nob,150,FALSE)#####sample 150 observations, leave out 27 observations for CV
obs=dd2[,3]*10^6     #######emulation of average diameter
input=as.data.frame(input)


Y=obs[s]
DM=input[s,1]
DM_new=input[-s,1][-1]
mod=lm(Y~log(DM))
lm_pred=predict(mod,newdata=data.frame(DM=DM_new))
Obs=mod$residual

#mm=mlegp(X=(DM),Z=Obs,constantMean=1)##
#dd=predict(mm,newData=as.matrix(DM_new),se.fit=TRUE)
################my method
source("C:\\Users\\olu\\Desktop\\NUFEB_materials\\GP_code1.R",echo=FALSE)
GP_pred=prediction(DM,Obs,DM_new)

emulex=list((lm_pred),GP_pred[,2])
lammp=observed[-s][-1];emu=emulex[[1]]
prop=1-(sum((lammp-emu)^2)/sum((lammp- mean(lammp))^2))

