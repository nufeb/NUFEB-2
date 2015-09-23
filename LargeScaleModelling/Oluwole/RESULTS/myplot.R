##############GP plot codes
myplot=function(observed,predicted,Z,sds=1,band=FALSE){
x1 = min(predicted - sds * sqrt(Z))
x1 = min(x1,observed)
x2 = max(predicted + sds * sqrt(Z))
x2 = max(x2, observed)
plot(observed, predicted,xlim = c(x1, x2), ylim = c(x1, x2), xlab = "observed",ylab = "predicted",main="Cross validation: Emulator vs observation",cex=1.5,cex.axis==1.5)
#index = order(observed)
index=1:21
lines(observed[index], observed[index], lty = 1)
if (band) {
lines(observed[index], predicted[index] - sds * sqrt(Z[index]),col = "red")
lines(observed[index], predicted[index] + sds * sqrt(Z[index]), col = "red")
}
else {
for (i in 1:length(observed)) {
lines(rep(observed[i],2), c(predicted[i] - sds*sqrt(Z[i]),predicted[i]+sds*sqrt(Z[i])),col = "red")
}}}
##########END
## open circle== GP predictions
## black lines== observations
## red lines== C.I
#myplot(y,fit$cv[,1],fit$cv[,2]) 

        


    