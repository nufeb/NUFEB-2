##THIS R CODE READS IN THE CLIMATE INPUT DATA
##########climgen relative to baseline (2005-2014) average
climgen_path=path
##############################################NEW

if(rcp=="other") {
 data1<- list.files(path=climgen_path,pattern=".nc")
} else {
data1<- list.files(path=climgen_path,pattern=sub("rcp",rcp,"rcp"))}

output="output"
cld=paste(climgen_path,data1[[1]],sep="/")
v3=paste(output,data1[[1]],sep="_")
v3=sub("output_cld_pat_CMIP3",output,paste(v3))
pre=paste(climgen_path,data1[[2]],sep="/")
tmp=paste(climgen_path,data1[[3]],sep="/")
wet=paste(climgen_path,data1[[4]],sep="/")
#############################
new_co2 <- list.files(path=climgen_path,pattern = ".txt")
####################### 
n=j
npixel=59199
load("rainfed/k3")
library(ncdf)
load("rainfed/N")
load("rainfed/grid_out")
###################################################
###COMPUTING SEASONAL VECTORS
#########################################
library(ncdf)
grid_len=npixel
cld=open.ncdf(cld)
kenny=cld[8]$dim$MONTH$vals
pre=open.ncdf(pre)
tmp=open.ncdf(tmp)
wet=open.ncdf(wet)
npixel=59199
w=sort(c(seq(1,120,12),seq(2,120,12),seq(12,120,12)))
s=sort(c(seq(6,120,12),seq(7,120,12),seq(8,120,12)))
sp=sort(c(seq(3,120,12),seq(4,120,12),seq(5,120,12)))
au=sort(c(seq(9,120,12),seq(10,120,12),seq(11,120,12)))

library(doParallel)
cl <- makeCluster(length(n))
registerDoParallel(cl)
tt<-{foreach(z=1:length(n),.combine =cbind) %do%{
T=matrix(0,nrow=59199,ncol=length(n))
T=k3+((z-1)*201600)}}

ind=c(tt)
my_list=list(cld,pre,tmp,wet)
load("rainfed/k3")

############COMPUTE SEASONAL "CLD_MONTHLY"
m=1
h=c("CLD_MONTHLY","PRECIP_MONTHLY","TMP_MONTHLY","WET_MONTHLY")
skip=12*(2005-kenny[1])+1
olu=c(get.var.ncdf(my_list[[m]],varid=h[m],start=c(skip,1,1),count=c(120*length(1),-1,-1)))
e=array(olu,dim=c(120,length(1),280,720,4))
X1=cbind(c(colMeans(e[s,,,,m])),c(colMeans(e[w,,,,m])),c(colMeans(e[sp,,,,m])),c(colMeans(e[au,,,,m])))
X2=array(X1,c(length(1),201600,4))
X2a=aperm(X2,c(2,1,3))
X2b=array(X2a,c(201600*length(1),4))
X2c=X2b[ind[1:59199],]
if(j==0) {
  X2d=rbind(X2c,X2c,X2c,X2c,X2c,X2c,X2c,X2c)
} else {
  X2d=X2c
}
if(j==0){
skip=skip+120
}else{
skip=skip+(120*n)
}
olu=c(get.var.ncdf(my_list[[m]],varid=h[m],start=c(skip,1,1),count=c(120*length(n),-1,-1)))
e=array(olu,dim=c(120,length(n),280,720,4))
X1=cbind(c(colMeans(e[s,,,,m])),c(colMeans(e[w,,,,m])),c(colMeans(e[sp,,,,m])),c(colMeans(e[au,,,,m])))
X2=array(X1,c(length(n),201600,4))
X2a=aperm(X2,c(2,1,3))
X3a=array(X2a,c(201600*length(n),4))
X3b=X3a[ind,]
cha=X3b-X2d
inp1=cbind(cha,X2d)
###########COMPUTE SEASONAL "PRECIP_MONTHLY"
m=2
skip=12*(2005-kenny[1])+1
olu=c(get.var.ncdf(my_list[[m]],varid=h[m],start=c(skip,1,1),count=c(120*length(1),-1,-1)))
e=array(olu,dim=c(120,length(1),280,720,4))
X1=cbind(c(colMeans(e[s,,,,m])),c(colMeans(e[w,,,,m])),c(colMeans(e[sp,,,,m])),c(colMeans(e[au,,,,m])))
X2=array(X1,c(length(1),201600,4))
X2a=aperm(X2,c(2,1,3))
X2b=array(X2a,c(201600*length(1),4))
X2c=X2b[ind[1:59199],]
if(j==0) {
  X2d=rbind(X2c,X2c,X2c,X2c,X2c,X2c,X2c,X2c)
} else {
  X2d=X2c
}
if(j==0){
skip=skip+120
}else{
skip=skip+(120*n)
}
olu=c(get.var.ncdf(my_list[[m]],varid=h[m],start=c(skip,1,1),count=c(120*length(n),-1,-1)))
e=array(olu,dim=c(120,length(n),280,720,4))
X1=cbind(c(colMeans(e[s,,,,m])),c(colMeans(e[w,,,,m])),c(colMeans(e[sp,,,,m])),c(colMeans(e[au,,,,m])))
X2=array(X1,c(length(n),201600,4))
X2a=aperm(X2,c(2,1,3))
X3a=array(X2a,c(201600*length(n),4))
X3b=X3a[ind,]
cha=X3b-X2d
inp2=cbind(cha,X2d)
###########COMPUTE SEASONAL "TMP_MONTHLY"
m=3
h=c("CLD_MONTHLY","PRECIP_MONTHLY","TMP_MONTHLY","WET_MONTHLY")
skip=12*(2005-kenny[1])+1
olu=c(get.var.ncdf(my_list[[m]],varid=h[m],start=c(skip,1,1),count=c(120*length(1),-1,-1)))
e=array(olu,dim=c(120,length(1),280,720,4))
X1=cbind(c(colMeans(e[s,,,,m])),c(colMeans(e[w,,,,m])),c(colMeans(e[sp,,,,m])),c(colMeans(e[au,,,,m])))
X2=array(X1,c(length(1),201600,4))
X2a=aperm(X2,c(2,1,3))
X2b=array(X2a,c(201600*length(1),4))
X2c=X2b[ind[1:59199],]
if(j==0) {
  X2d=rbind(X2c,X2c,X2c,X2c,X2c,X2c,X2c,X2c)
} else {
  X2d=X2c
}
if(j==0){
skip=skip+120
}else{
skip=skip+(120*n)
}
olu=c(get.var.ncdf(my_list[[m]],varid=h[m],start=c(skip,1,1),count=c(120*length(n),-1,-1)))
e=array(olu,dim=c(120,length(n),280,720,4))
X1=cbind(c(colMeans(e[s,,,,m])),c(colMeans(e[w,,,,m])),c(colMeans(e[sp,,,,m])),c(colMeans(e[au,,,,m])))
X2=array(X1,c(length(n),201600,4))
X2a=aperm(X2,c(2,1,3))
X3a=array(X2a,c(201600*length(n),4))
X3b=X3a[ind,]
cha=X3b-X2d
inp3=cbind(cha,X2d)
#########COMPUTE SEASONAL WET_MONTHLY
m=4
skip=12*(2005-kenny[1])+1
olu=c(get.var.ncdf(my_list[[m]],varid=h[m],start=c(skip,1,1),count=c(120*length(1),-1,-1)))
e=array(olu,dim=c(120,length(1),280,720,4))
X1=cbind(c(colMeans(e[s,,,,m])),c(colMeans(e[w,,,,m])),c(colMeans(e[sp,,,,m])),c(colMeans(e[au,,,,m])))
X2=array(X1,c(length(1),201600,4))
X2a=aperm(X2,c(2,1,3))
X2b=array(X2a,c(201600*length(1),4))
X2c=X2b[ind[1:59199],]
if(j==0) {
  X2d=rbind(X2c,X2c,X2c,X2c,X2c,X2c,X2c,X2c)
} else {
  X2d=X2c
}
if(j==0){
skip=skip+120
}else{
skip=skip+(120*n)
}
olu=c(get.var.ncdf(my_list[[m]],varid=h[m],start=c(skip,1,1),count=c(120*length(n),-1,-1)))
e=array(olu,dim=c(120,length(n),280,720,4))
X1=cbind(c(colMeans(e[s,,,,m])),c(colMeans(e[w,,,,m])),c(colMeans(e[sp,,,,m])),c(colMeans(e[au,,,,m])))
X2=array(X1,c(length(n),201600,4))
X2a=aperm(X2,c(2,1,3))
X3a=array(X2a,c(201600*length(n),4))
X3b=X3a[ind,]
cha=X3b-X2d
inp4=cbind(cha,X2d)
dat=cbind(inp1,inp2,inp3,inp4)

dat=dat[,c(1:4,9:12,17:20,25:28,5:8,13:16,21:24,29:32)]

load("rainfed/other_inp2")

if(rcp=="RCP3PD") {
 clim=other_inp2[[1]][,n,]
} else {
if(rcp=="RCP45") {
 clim=other_inp2[[2]][,n,]
} else {
if(rcp=="RCP6") {
 clim=other_inp2[[3]][,n,]
} else {
if(rcp=="RCP85") {
 clim=other_inp2[[4]][,n,]
}else{
source("rainfed/co2_file.R",echo=TRUE)
clim=cbind(co2_dat,other_inp2[[1]][,n,3:4])
}}}}
w=j
s2=cbind(dat,clim)
#############################################
#########################end
