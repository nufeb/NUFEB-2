##############################wls using WLS

if(out=="current"){
load("rainfed2/masked/crop_err4c")
load("rainfed2/masked/N1")
T=N1
} else {
load("rainfed2/crop_N")
load("rainfed2/unmasked/crop_err4c")
T=crop_N
}
##########
v4=sub("man",man,(sub("j",j,"p_j_man")))
load(paste("rainfed2/pred",v4,sep="/"))
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
X=m1$x[,1:4]
#####=t(t(l_cereal[,1:4])%*%ada)
y_star=prediction[T[[m]],m]
X_star=predict(m1,t(as.matrix(y_star)))[,1:4]


lambda=(m1$sdev)^2/sum((m1$sdev)^2)
#lambda=(m1$sdev)^2
d=(sqrt(colSums(lambda*(matrix(x_star,nrow=4,ncol=16)-x_i)^2)))
#d=sort(d)
w_i=(1/d^2)/sum((1/d^2))
test=any(is.na(w_i))
if(test==TRUE){
ind=which(is.na(w_i))
if(length(ind)==1){E_crop1=E[,ind] 
}else {E_crop1=rowMeans(E[,ind])}} else {
X=t(x_i)
W=diag(w_i)
beta=ginv(t(X)%*%W%*%X)%*%t(X)%*%W%*%t(E)
E_star=t(beta)%*%x_star
r1=apply(E,1,min)
r2=apply(E,1,max)
r3=rowMeans(E)
EE=pmin(c(E_star),(r2))
E_crop1=pmax(EE,(r1))}


