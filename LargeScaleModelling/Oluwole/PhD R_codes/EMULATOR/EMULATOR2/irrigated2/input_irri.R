###cimgen relative to baseline (2005-2014) average
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
load("irrigated2/k3")
library(ncdf)
load("irrigated2/N")
load("irrigated2/grid_out")
###################################################
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
load("irrigated2/k3")

#####################################
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
###########
m=2
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
inp2=cbind(cha,X2d)
###########
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
#########
m=4
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
inp4=cbind(cha,X2d)
dat=cbind(inp1,inp2,inp3,inp4)
dat=dat[,c(1:4,9:12,17:20,25:28,5:8,13:16,21:24,29:32)]
load("irrigated2/other_inp")

if(rcp=="RCP3PD") {
 clim=other_inp[[1]]
} else {
if(rcp=="RCP45") {
 clim=other_inp[[2]]
} else {
if(rcp=="RCP6") {
 clim=other_inp[[3]]
} else {
if(rcp=="RCP85") {
 clim=other_inp[[4]]
}else{
source("irrigated2/co2_file.R",echo=TRUE)
clim=cbind(other_inp[[1]][,3:4],co2_dat)
}}}}
w=j
cm2=array(clim,c(npixel,8,4))
#######################co2 reverse relative to 1st decade
####dim=c(59199,8,4)
co2_9=cm2[,8,1]+cm2[,8,2]
co2_ini=rep(cm2[,1,2],8)#####co2
co2_all=c(cm2[,1:8,2],co2_9)
co2_b=co2_all[59200:532791]-co2_ini####cco2
my_co2=c(co2_b,co2_ini)
cm2[,,1:2]=my_co2
##########################
cm3=cbind(c(cm2[,n,1]),c(cm2[,n,2]),c(cm2[,n,3]),c(cm2[,n,4]))
clim=cm3
load("irrigated2/inni")
inn=array(inni,c(npixel,8,4))
wala=cbind(c(inn[,n,1]),c(inn[,n,2]),c(inn[,n,3]),c(inn[,n,4]))
s2=cbind(wala,dat,clim)
s2=array(s2,c(npixel,length(n),40))
close.ncdf(pre)
close.ncdf(tmp)
close.ncdf(wet)
close.ncdf(cld)


rm(au,cha,cld,clim,cm2,cm3,co2_9,co2_all,co2_b,co2_ini,wala,dat,ind,inn,inp1,inp2,inp3,inp4,my_co2,new_co2,olu,other_inp,pre,s,tmp,sp,T,X1,X2,X2a,X2b,X2c,X2d,X3a,X3b,z,wet)
#########################end

