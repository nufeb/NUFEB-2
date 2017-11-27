######################R main emulator codes
library(MASS)
source("abind.R",echo=FALSE)
#library(abind)
##predict nutrient concentration overtime
#coef=as.matrix(read.table('coef.txt'))
new=as.matrix(read.table('input/input.txt',header=FALSE))
#new=as.matrix(read.table('input/Xstar.csv',header=TRUE))
namm0=c("nh3","no2","no3","o2","co2","Biomass","nparticle","Height","Roughness")
 X=as.matrix(read.table('X.csv',header=TRUE))##training data
Y=as.matrix(read.table('Y.csv',header=TRUE))##tr
namm=colnames(X)
s2=c(2.166465e-15,2.079999e+01,4.336633e-10,4.304684e-10) 
s1=c(2.467027e-15,2.339083e+01,8.968750e-10,6.864521e-11) 
Y=scale(Y,scale=s2,center=s1)
X=log(X)
input=new[,drop=FALSE]
input=log(input)
###compute correlation
mycor=function(X1,X2,betai){
pdf<- diag(betai, nrow = length(betai))
tx1<- t(X1)
tx2 <- t(X2)
R1 <- X1%*%pdf%*%tx1 
R2 <- X2%*%pdf%*%tx2
S1 <- t(as.matrix(diag(R1)))%x%rep(1,nrow(X2))
S2 <- as.matrix(diag(R2))%x%t(rep(1,nrow(X1)))
a1 <- t(tx1)%*%pdf%*%tx2
a2 <- t(tx2)%*%pdf%*%tx1
return(exp(Re(t(a1)+a2-S1-S2)))
}
############
nout=4
tau=0.5
#scale2=c( 0.201935506, 0.006313825, 0.000010000, 0.177466408, 0.000010000, 6.577416780, 7.721354529, 6.594165598,4.081750274, 0.004609465)
scale2=c(0.003847799,0.001000000,0.001000000,0.077725734,0.001000000,6.862597279,4.714029223,2.095145620,4.586664423,0.001000000)
id2=length(scale2)
theta=scale2[-id2]	
#load("Ainv")
load("A/A0")
load("A/A1")
load("A/A2")
load("A/A3")
load("A/A4")
load("A/A5")
load("A/A6")
load("A/A7")
Ainv=cbind(A0,A1,A2,A3,A4,A5,A6,A7)
#
pred=function(newdata){
X0=newdata
#colnames(X0)=namm
form=~nh3+no2+no3+o2+co2+Biomass+nparticle+Height+Roughness
H<- X# model.matrix(form,data=as.data.frame(X))
m<-dim(H)[2]
n<-dim(Y)[1]
n2<- dim(newdata)[1]
######model matrix 
H0<-X0#model.matrix(form,as.data.frame(X0))
A01=mycor(X0,X,betai=theta)###cross correlation
A00=mycor(X0,X0,betai=theta)##test point correlation
A00=A00+(scale2[id2]*diag(tau,dim(A00)))
iOmega=solve(t(H)%*%Ainv%*%H)
betahat=solve(t(H)%*%Ainv%*%H)%*%(t(H)%*%Ainv%*%Y)
mu_star=H0%*%betahat+t(A01)%*%Ainv%*%(Y-H%*%betahat)
r1=H0-(t(A01)%*%Ainv%*%H)
c_star=A00-(t(A01)%*%Ainv%*%A01)+(((r1)%*%(iOmega)%*%t(r1)))
Sigma=(t(Y-H%*%betahat)%*%Ainv%*%(Y-H%*%betahat))/(n-m)
svar=list()
for(i in 1:n2){
svar[[i]]=(diag(c_star)[i]*diag(Sigma))
}
K0=abind(svar,along=2)
out=list(mu=mu_star,K=t(K0))
return(out)
}
##
gap=pred(input)
out1=gap[[1]]
out2=(gap[[2]])
upper=out1+(2*sqrt(out2))
lower=out1-(2*sqrt(out2))

low=upp=drdt=out1
for(i in 1:4){
	drdt[,i]=(out1[,i]*s2[i])+s1[i]
	upp[,i]=(upper[,i]*s2[i])+s1[i]
	low[,i]=(lower[,i]*s2[i])+s1[i]
	#olu[,i]=(olu[,i]*s2[i])+s1[i]
}
VV=((upp-low)/4)
colnames(drdt)=colnames(Y)
colnames(VV)=colnames(out1)=colnames(out2)=colnames(Y)
#write.table(VV,file="output/deviation.txt",row.names=FALSE)
write.table(drdt,file="output/Rate.txt",row.names=FALSE,col.names = TRUE)
