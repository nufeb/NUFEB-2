#########################Run the emulator
source("access/CovM1M2.R",echo=TRUE)
load("access/mod")
source("access/input.R",echo=TRUE)
#load("access/form")
source("access/form.R")
form=list(c1,c2,c3,c4,c5)
newdata=DM_new
#load("access/dat")
output=list()
for(i in 1:5){
mm=mod[[i]]
beta=mm[[1]];X=mm[[2]];y=mm[[3]];z=mm[[4]];M=mm[[5]];n=mm[[6]];p=mm[[7]];mc=mm[[8]]
#class(mc)=NULL
#detach("package:DiceKriging",unload=TRUE)
aux <- covMatrix(mc,X=X, noise.var=model@noise.var)
C <- aux[[1]]
T <- chol(C)
F.newdata <- model.matrix(form[[i]], data = data.frame(newdata))
y.predict.linear <- F.newdata%*%beta
c.newdata <- covMat1Mat2(mc, X1 =as.matrix(X), X2 =as.matrix(newdata),nugget.flag = mc@nugget.flag)
Tinv.c.newdata <- backsolve(t(T), c.newdata, upper.tri=FALSE)
y.predict.c <- t(Tinv.c.newdata)%*%z
y.predict <- as.numeric(y.predict.linear + y.predict.c)
################Standard error of predictions
total.sd2 <- mc@sd2
total.sd2 <- total.sd2 + mc@nugget
s2.predict.1 <- apply(Tinv.c.newdata, 2,crossprod)#compute c(x)'*C^(-1)*c(x)for x = newdata
#if (type == "SK") {
#s2.predict <- pmax(total.sd2 - s2.predict.1, 0)
#s2.predict <- as.numeric(s2.predict)
#q95 <- qnorm(0.975)
#}
#else if (type == "UK") {
T.M <- chol(t(M)%*%M)   # equivalently : qrR <- qr.R(qr(M))
s2.predict.mat <- backsolve(t(T.M), t(F.newdata - t(Tinv.c.newdata)%*%M) , upper.tri = FALSE)
s2.predict.2 <- apply(s2.predict.mat, 2, crossprod)
s2.predict <- pmax(total.sd2 - s2.predict.1 + s2.predict.2, 0)
s2.predict <- as.numeric(s2.predict)
#if (bias.correct) s2.predict <- s2.predict * object@n/(object@n - object@p)
q95 <- qt(0.975,n-p)
#}
lower95 <- y.predict - q95*sqrt(s2.predict)
upper95 <- y.predict + q95*sqrt(s2.predict)
varr=cbind(lower95,upper95)
sd <- sqrt(s2.predict)
outs=list()
outs[[1]]=y.predict;outs[[2]]=sd;outs[[3]]=varr
output[[i]]=outs
}
#########Output
#res1=output[[1]][[1]]### EPS MASS
res2=output[[2]][[1]]###TOTAL FLOC MASS
res3=output[[3]][[1]]##FLOC DIAM (VOL)
res4=round(output[[4]][[1]])###NO of ATOMS
#res5=output[[5]][[1]]###FLOC DIAM (distance)
time=cap=seq(0,(runtime-step),step)
text <- readLines("Output/esnap",encoding="UTF-8")
#ind=c(2,4,10)
des=cbind(time,res4,res3,res2)
options(scipen=2)
tt=list()
for(i in 1:length(time)){
#for(i in 1:4){
text[2]=gsub(text[2],des[i,1], text[2])
text[4]=gsub(text[4],des[i,2], text[4])
text[10]=gsub(text[10],paste("Floc",des[i,3],des[i,4],sep=" "), text[10])
tt[[i]]=text
}
library(abind)
text=abind(tt)
dirname=paste(getwd(),"Output",sep="/")
writeLines(text,paste("snapshot", ".bubblemd",sep=""))
file.copy("snapshot.bubblemd",dirname)
file.remove("snapshot.bubblemd")

########################################################END

