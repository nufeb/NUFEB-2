#########################Run the emulator
source("access/CovM1M2",echo=TRUE)
load("access/mod")
source("access/input.R",echo=TRUE)
newdata=DM_new
#load("access/dat")
output=list()
for(i in 1:5){
mm=mod[[i]]
beta=mm[[1]];form=mm[[2]];X=mm[[3]];y=mm[[4]];z=mm[[5]];M=mm[[6]];n=mm[[7]];p=mm[[8]];mc=mm[[9]]
aux <- covMatrix(mc,X=X, noise.var=model@noise.var)
C <- aux[[1]]
T <- chol(C)
F.newdata <- model.matrix(form, data = data.frame(newdata))
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

########################################################END

