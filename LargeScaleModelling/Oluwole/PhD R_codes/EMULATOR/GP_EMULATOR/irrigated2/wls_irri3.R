##############################wls using WLS
library(MASS)
if(out=="current"){
load("irrigated2/masked/crop_err4c")
load("irrigated2/masked/N1")
T=N1
} else {
load("irrigated2/crop_N")
load("irrigated2/unmasked/crop_err4c")
T=crop_N
}
##########
v4=sub("man",man,(sub("j",j,"p_j_man")))
load(paste("irrigated2/pred",v4,sep="/"))
crop_pred4b=eval(parse(text = as.name(v4)))
m=i
w=j##############correction
###NEW CROP ANALYSIS REL TO 1ST DEC" 2nd STAGE EMUL
##########

m=1
s1=array(crop_pred4b[[m]],c(59199,4*5))
l=186
e1=crop_err4c[[m]]
E=(e1[,w,,man])
E=array(E,c(186,4*5))
Y2=array(s1,c(59199,4*5))
m1=prcomp(t(Y2[T[[m]],]),scale=TRUE,center=TRUE)
l_cereal=m1$rotation
##ada=(Y2[T[[m]],]-m1$center)/matrix(rep(m1$scale,20),ncol=20)
X=t(m1$x[,1:4])
#####=t(t(l_cereal[,1:4])%*%ada)
y_star=prediction[T[[m]],m]
X_star=t(t(predict(m1,t(as.matrix(y_star)))[,1:4]))

lambda=(m1$sdev)^2/sum((m1$sdev)^2)
#lambda=(m1$sdev)^2
d=(sqrt(colSums(lambda*(matrix(X_star,nrow=4,ncol=20)-X)^2)))
#d=sort(d)
w_i=(1/d^2)/sum((1/d^2))
test=any(is.na(w_i))
if(test==TRUE){
ind=which(is.na(w_i))
if(length(ind)==1){E_crop1=E[,ind] 
}else {E_crop1=rowMeans(E[,ind])}} else {
X=t(X)
W=diag(w_i)
beta=ginv(t(X)%*%W%*%X)%*%t(X)%*%W%*%t(E)
E_star=t(beta)%*%X_star
r1=apply(E,1,min)
r2=apply(E,1,max)
r3=rowMeans(E)
EE=pmin(c(E_star),(r2))
E_crop1=pmax(EE,(r1))}
#############
m=2
s1=array(crop_pred4b[[m]],c(59199,4*5))
l=186
e1=crop_err4c[[m]]
E=(e1[,w,,man])
E=array(E,c(186,4*5))
Y2=array(s1,c(59199,4*5))
m1=prcomp(t(Y2[T[[m]],]),scale=TRUE,center=TRUE)
l_cereal=m1$rotation
##ada=(Y2[T[[m]],]-m1$center)/matrix(rep(m1$scale,20),ncol=20)
X=t(m1$x[,1:4])
#####=t(t(l_cereal[,1:4])%*%ada)
y_star=prediction[T[[m]],m]
X_star=t(t(predict(m1,t(as.matrix(y_star)))[,1:4]))

lambda=(m1$sdev)^2/sum((m1$sdev)^2)
#lambda=(m1$sdev)^2
d=(sqrt(colSums(lambda*(matrix(X_star,nrow=4,ncol=20)-X)^2)))
#d=sort(d)
w_i=(1/d^2)/sum((1/d^2))
test=any(is.na(w_i))
if(test==TRUE){
ind=which(is.na(w_i))
if(length(ind)==1){E_crop1=E[,ind] 
}else {E_crop1=rowMeans(E[,ind])}} else {
X=t(X)
W=diag(w_i)
beta=ginv(t(X)%*%W%*%X)%*%t(X)%*%W%*%t(E)
E_star=t(beta)%*%X_star
r1=apply(E,1,min)
r2=apply(E,1,max)
r3=rowMeans(E)
EE=pmin(c(E_star),(r2))
E_crop2=pmax(EE,(r1))}
####################
m=3
s1=array(crop_pred4b[[m]],c(59199,4*5))
l=186
e1=crop_err4c[[m]]
E=(e1[,w,,man])
E=array(E,c(186,4*5))
Y2=array(s1,c(59199,4*5))
m1=prcomp(t(Y2[T[[m]],]),scale=TRUE,center=TRUE)
l_cereal=m1$rotation
##ada=(Y2[T[[m]],]-m1$center)/matrix(rep(m1$scale,20),ncol=20)
X=t(m1$x[,1:4])
#####=t(t(l_cereal[,1:4])%*%ada)
y_star=prediction[T[[m]],m]
X_star=t(t(predict(m1,t(as.matrix(y_star)))[,1:4]))

lambda=(m1$sdev)^2/sum((m1$sdev)^2)
#lambda=(m1$sdev)^2
d=(sqrt(colSums(lambda*(matrix(X_star,nrow=4,ncol=20)-X)^2)))
#d=sort(d)
w_i=(1/d^2)/sum((1/d^2))
test=any(is.na(w_i))
if(test==TRUE){
ind=which(is.na(w_i))
if(length(ind)==1){E_crop1=E[,ind] 
}else {E_crop1=rowMeans(E[,ind])}} else {
X=t(X)
W=diag(w_i)
beta=ginv(t(X)%*%W%*%X)%*%t(X)%*%W%*%t(E)
E_star=t(beta)%*%X_star
r1=apply(E,1,min)
r2=apply(E,1,max)
r3=rowMeans(E)
EE=pmin(c(E_star),(r2))
E_crop3=pmax(EE,(r1))}
############
m=4
s1=array(crop_pred4b[[m]],c(59199,4*5))
l=186
e1=crop_err4c[[m]]
E=(e1[,w,,man])
E=array(E,c(186,4*5))
Y2=array(s1,c(59199,4*5))
m1=prcomp(t(Y2[T[[m]],]),scale=TRUE,center=TRUE)
l_cereal=m1$rotation
##ada=(Y2[T[[m]],]-m1$center)/matrix(rep(m1$scale,20),ncol=20)
X=t(m1$x[,1:4])
#####=t(t(l_cereal[,1:4])%*%ada)
y_star=prediction[T[[m]],m]
X_star=t(t(predict(m1,t(as.matrix(y_star)))[,1:4]))

lambda=(m1$sdev)^2/sum((m1$sdev)^2)
#lambda=(m1$sdev)^2
d=(sqrt(colSums(lambda*(matrix(X_star,nrow=4,ncol=20)-X)^2)))
#d=sort(d)
w_i=(1/d^2)/sum((1/d^2))
test=any(is.na(w_i))
if(test==TRUE){
ind=which(is.na(w_i))
if(length(ind)==1){E_crop1=E[,ind] 
}else {E_crop1=rowMeans(E[,ind])}} else {
X=t(X)
W=diag(w_i)
beta=ginv(t(X)%*%W%*%X)%*%t(X)%*%W%*%t(E)
E_star=t(beta)%*%X_star
r1=apply(E,1,min)
r2=apply(E,1,max)
r3=rowMeans(E)
EE=pmin(c(E_star),(r2))
E_crop4=pmax(EE,(r1))}
E_crop=list(E_crop1,E_crop2,E_crop3,E_crop4)
detach("package:MASS", unload=TRUE)
source("irrigated2/aggregate3.R",echo=FALSE)

#rm(E_crop1,E_crop2,E_crop3,E_crop4,crop_err4,crop_pred4,beta,E_star,X,x_star,W,Gamma,x_i,Y2,y_star)


