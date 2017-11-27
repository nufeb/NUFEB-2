##############################wls
if(out=="current"){
load("irrigated2/masked/crop_err4c")
load("irrigated2/masked/N2")
T=N2
} else {
load("irrigated2/crop_N")
load("irrigated2/crop_err4c")
T=crop_N
}
##########
crop_pred4b=list(rep(0,7))
for(man in 1:7){
v4=sub("man",man,(sub("j",j,"p_j_man")))
load(paste("irrigated2/pred",v4,sep="/"))
crop_pred4b[[man]]=eval(parse(text = as.name(v4)))
}
#rm(man)
library(mlegp)
m=i
w=j##############correction
###NEW CROP ANALYSIS REL TO 1ST DEC" 2nd STAGE EMUL
#######################################
m=1
crop_err4c[[m]]=aperm(crop_err4c[[m]],c(1,2,3,5,4))
D_1=crop_pred4b[[1]][[m]]
D_2=crop_pred4b[[2]][[m]]
D_3=crop_pred4b[[3]][[m]]
D_4=crop_pred4b[[4]][[m]]
D_5=crop_pred4b[[5]][[m]]
D_6=crop_pred4b[[6]][[m]]
D_7=crop_pred4b[[7]][[m]]
library(abind)
s1=abind(D_1,D_2,D_3,D_4,D_5,D_6,D_7,along=4)

#s1=array(crop_pred4b[[m]],c(59199,12))
l=186
#e1=array(crop_err4c[[m]],c(l,8,12,7))
e1=crop_err4c[[m]]
######
E=(e1[,w,,,])
E=array(E,c(186,6*2*7))
Y2=array(s1,c(59199,6*2*7))
#m1=prcomp(Y2[N1[[m]],-c(5,11)],scale=F,center=F)
m1=prcomp(Y2,scale=F,center=F)
l_cereal=m1$rotation ###for 84 by 84
y_star=prediction[T[[m]],m]
Gamma=Y2[T[[m]],]%*%l_cereal ###for
x_i=t(Gamma[,1:4])%*%Y2[T[[m]],]
x_star=t(Gamma[,1:4])%*%y_star
X=t(x_i)
X_star=t(x_star)
###############fit GP emulator
residual=t(E)
ind=which(residual[m,]!=0) 
u=which(duplicated(X))
if(length(u)==0){
fitPC = mlegp(X=X,Z=residual[,ind],constantMean=1,simplex.ntries=3,simplex.maxiter =100)
}else{
fitPC = mlegp(X=X[-u,],Z=residual[-u,ind],constantMean=1,simplex.ntries=3,simplex.maxiter =100)}
id=length(ind)
pred=rep(0,id)
for(r in 1:id){
pred[r]=predict(fitPC[[r]],X_star)
}
predY=rep(0,186)
predY[ind]=pred
E_crop1=predY
rm(data1,D_1,D_2,D_3,D_4,D_5,D_6,D_7,s1)
#######################################
m=2
crop_err4c[[m]]=aperm(crop_err4c[[m]],c(1,2,3,5,4))
D_1=crop_pred4b[[1]][[m]]
D_2=crop_pred4b[[2]][[m]]
D_3=crop_pred4b[[3]][[m]]
D_4=crop_pred4b[[4]][[m]]
D_5=crop_pred4b[[5]][[m]]
D_6=crop_pred4b[[6]][[m]]
D_7=crop_pred4b[[7]][[m]]
library(abind)
s1=abind(D_1,D_2,D_3,D_4,D_5,D_6,D_7,along=4)

#s1=array(crop_pred4b[[m]],c(59199,12))
l=186
#e1=array(crop_err4c[[m]],c(l,8,12,7))
e1=crop_err4c[[m]]
######
E=(e1[,w,,,])
E=array(E,c(186,6*2*7))
Y2=array(s1,c(59199,6*2*7))
#m1=prcomp(Y2[N1[[m]],-c(5,11)],scale=F,center=F)
m1=prcomp(Y2,scale=F,center=F)
l_cereal=m1$rotation ###for 84 by 84
y_star=prediction[T[[m]],m]
Gamma=Y2[T[[m]],]%*%l_cereal ###for
x_i=t(Gamma[,1:4])%*%Y2[T[[m]],]
x_star=t(Gamma[,1:4])%*%y_star
X=t(x_i)
X_star=t(x_star)
###############fit GP emulator
residual=t(E)
ind=which(residual[m,]!=0) 
u=which(duplicated(X))
if(length(u)==0){
fitPC = mlegp(X=X,Z=residual[,ind],constantMean=1,simplex.ntries=3,simplex.maxiter =100)
}else{
fitPC = mlegp(X=X[-u,],Z=residual[-u,ind],constantMean=1,simplex.ntries=3,simplex.maxiter =100)}
id=length(ind)
pred=rep(0,id)
for(r in 1:id){
pred[r]=predict(fitPC[[r]],X_star)
}
predY=rep(0,186)
predY[ind]=pred
E_crop2=predY
rm(D_1,D_2,D_3,D_4,D_5,D_6,D_7,s1)
#######################################
m=3
crop_err4c[[m]]=aperm(crop_err4c[[m]],c(1,2,3,5,4))
D_1=crop_pred4b[[1]][[m]]
D_2=crop_pred4b[[2]][[m]]
D_3=crop_pred4b[[3]][[m]]
D_4=crop_pred4b[[4]][[m]]
D_5=crop_pred4b[[5]][[m]]
D_6=crop_pred4b[[6]][[m]]
D_7=crop_pred4b[[7]][[m]]
library(abind)
s1=abind(D_1,D_2,D_3,D_4,D_5,D_6,D_7,along=4)

#s1=array(crop_pred4b[[m]],c(59199,12))
l=186
#e1=array(crop_err4c[[m]],c(l,8,12,7))
e1=crop_err4c[[m]]
######
E=(e1[,w,,,])
E=array(E,c(186,6*2*7))
Y2=array(s1,c(59199,6*2*7))
#m1=prcomp(Y2[N1[[m]],-c(5,11)],scale=F,center=F)
m1=prcomp(Y2,scale=F,center=F)
l_cereal=m1$rotation ###for 84 by 84
y_star=prediction[T[[m]],m]
Gamma=Y2[T[[m]],]%*%l_cereal ###for
x_i=t(Gamma[,1:4])%*%Y2[T[[m]],]
x_star=t(Gamma[,1:4])%*%y_star
X=t(x_i)
X_star=t(x_star)
###############fit GP emulator
residual=t(E)
ind=which(residual[m,]!=0) 
u=which(duplicated(X))
if(length(u)==0){
fitPC = mlegp(X=X,Z=residual[,ind],constantMean=1,simplex.ntries=3,simplex.maxiter =150)
}else{
fitPC = mlegp(X=X[-u,],Z=residual[-u,ind],constantMean=1,simplex.ntries=3,simplex.maxiter =150)}
id=length(ind)
pred=rep(0,id)
for(r in 1:id){
pred[r]=predict(fitPC[[r]],X_star)
}
predY=rep(0,186)
predY[ind]=pred
E_crop3=predY
rm(D_1,D_2,D_3,D_4,D_5,D_6,D_7,s1)
#######################################
m=4
crop_err4c[[m]]=aperm(crop_err4c[[m]],c(1,2,3,5,4))
D_1=crop_pred4b[[1]][[m]]
D_2=crop_pred4b[[2]][[m]]
D_3=crop_pred4b[[3]][[m]]
D_4=crop_pred4b[[4]][[m]]
D_5=crop_pred4b[[5]][[m]]
D_6=crop_pred4b[[6]][[m]]
D_7=crop_pred4b[[7]][[m]]
library(abind)
s1=abind(D_1,D_2,D_3,D_4,D_5,D_6,D_7,along=4)

#s1=array(crop_pred4b[[m]],c(59199,12))
l=186
#e1=array(crop_err4c[[m]],c(l,8,12,7))
e1=crop_err4c[[m]]
######
E=(e1[,w,,,])
E=array(E,c(186,6*2*7))
Y2=array(s1,c(59199,6*2*7))
#m1=prcomp(Y2[N1[[m]],-c(5,11)],scale=F,center=F)
m1=prcomp(Y2,scale=F,center=F)
l_cereal=m1$rotation ###for 84 by 84
y_star=prediction[T[[m]],m]
Gamma=Y2[T[[m]],]%*%l_cereal ###for
x_i=t(Gamma[,1:4])%*%Y2[T[[m]],]
x_star=t(Gamma[,1:4])%*%y_star
X=t(x_i)
X_star=t(x_star)
###############fit GP emulator
residual=t(E)
ind=which(residual[m,]!=0) 
u=which(duplicated(X))
if(length(u)==0){
fitPC = mlegp(X=X,Z=residual[,ind],constantMean=1,simplex.ntries=3,simplex.maxiter =150)
}else{
fitPC = mlegp(X=X[-u,],Z=residual[-u,ind],constantMean=1,simplex.ntries=3,simplex.maxiter =100)}
id=length(ind)
pred=rep(0,id)
for(r in 1:id){
pred[r]=predict(fitPC[[r]],X_star)
}
predY=rep(0,186)
predY[ind]=pred
E_crop4=predY
#######################################

##############
E_crop=list(E_crop1,E_crop2,E_crop3,E_crop4)
source("irrigated2/aggregate3.R",echo=FALSE)
rm(D_1,D_2,D_3,D_4,D_5,D_6,D_7,s1)
rm(E_crop1,E_crop2,E_crop3,E_crop4,crop_err4c,X,x_star,Gamma,x_i,Y2,y_star,predY,fitPC,residual,m1,E,e1,l_cereal,input)
