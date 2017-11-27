oja1=ini[[1]][,1:2,,,]
oja2=ini[[2]][,1:2,,,]
oja3=ini[[3]][,1:2,,,]
oja4=ini[[4]][,1:2,,,]
oja1=array(oja1,c(59199,4,7,2))
oja2=array(oja2,c(59199,4,7,2))
oja3=array(oja3,c(59199,4,7,2))
oja4=array(oja4,c(59199,4,7,2))
ini=list(oja1,oja2,oja3,oja4)

load("E:/WOLE/validation/S_irri_lpj")#
load("E:/WOLE/validation/S_irri_lpj")#

cor(crop[,1],S_irri_lpj[[1]][,8,4,5])^2##.98
cor(crop[,2],S_irri_lpj[[2]][,8,4,5])^2#.98
cor(crop[,3],S_irri_lpj[[3]][,8,4,5])^2#.24
cor(crop[,4],S_irri_lpj[[4]][,8,4,5])^2#.69

cor(crop[,1],S_irri_lpj[[1]][,5,4,6])^2##.98
cor(crop[,2],S_irri_lpj[[2]][,5,4,6])^2#.998
cor(crop[,3],S_irri_lpj[[3]][,5,4,6])^2#.033
cor(crop[,4],S_irri_lpj[[4]][,5,4,6])^2#.60e ’library(help="mlegp")’.

cor(crop[,1],S_irri_lpj[[1]][,5,3,2])^2##.98
cor(crop[,2],S_irri_lpj[[2]][,5,3,2])^2#.998
cor(crop[,3],S_irri_lpj[[3]][,5,3,2])^2#.025
cor(crop[,4],S_irri_lpj[[4]][,5,3,2])^2#.63.


 cor(crop[,3],S_irri_lpj[[3]][,6,3,1])
cor(crop[,3],S_irri_lpj[[3]][,1,3,1])


source("LPJML_crop-yield_code_v3f.R",echo=TRUE)

j=1

cor(crop[,1],S_irri_lpj[[1]][,j,3,1])^2
cor(crop[,2],S_irri_lpj[[2]][,j,3,1])^2
cor(crop[,3],S_irri_lpj[[3]][,j,3,1])^2
cor(crop[,4],S_irri_lpj[[4]][,j,3,1])^2

j=2
cor(crop[,5],S_irri_lpj[[1]][,j,3,1])^2
cor(crop[,6],S_irri_lpj[[2]][,j,3,1])^2
cor(crop[,7],S_irri_lpj[[3]][,j,3,1])^2
cor(crop[,8],S_irri_lpj[[4]][,j,3,1])^2
j=3
cor(crop[,9],S_irri_lpj[[1]][,j,3,1])^2
cor(crop[,10],S_irri_lpj[[2]][,j,3,1])^2
cor(crop[,11],S_irri_lpj[[3]][,j,3,1])^2
cor(crop[,12],S_irri_lpj[[4]][,j,3,1])^2

j=4
cor(crop[,13],S_irri_lpj[[1]][,j,3,1])^2
cor(crop[,14],S_irri_lpj[[2]][,j,3,1])^2
cor(crop[,15],S_irri_lpj[[3]][,j,3,1])^2
cor(crop[,16],S_irri_lpj[[4]][,j,3,1])^2


j=8
cor(crop[,29],S_irri_lpj[[1]][,j,3,1])^2
cor(crop[,30],S_irri_lpj[[2]][,j,3,1])^2
cor(crop[,31],S_irri_lpj[[3]][,j,3,1])^2
cor(crop[,32],S_irri_lpj[[4]][,j,3,1])^2

#############################

#########################################
j=1

cor(crop[,1],S_rain_lpj[[1]][,j,3,1])^2
cor(crop[,2],S_rain_lpj[[2]][,j,3,1])^2
cor(crop[,3],S_rain_lpj[[3]][,j,3,1])^2
cor(crop[,4],S_rain_lpj[[4]][,j,3,1])^2

j=2
cor(crop[,5],S_rain_lpj[[1]][,j,3,1])^2
cor(crop[,6],S_rain_lpj[[2]][,j,3,1])^2
cor(crop[,7],S_rain_lpj[[3]][,j,3,1])^2
cor(crop[,8],S_rain_lpj[[4]][,j,3,1])^2
j=3
cor(crop[,9],S_rain_lpj[[1]][,j,3,1])^2
cor(crop[,10],S_rain_lpj[[2]][,j,3,1])^2
cor(crop[,11],S_rain_lpj[[3]][,j,3,1])^2
cor(crop[,12],S_rain_lpj[[4]][,j,3,1])^2

j=4
cor(crop[,13],S_rain_lpj[[1]][,j,3,1])^2
cor(crop[,14],S_rain_lpj[[2]][,j,3,1])^2
cor(crop[,15],S_rain_lpj[[3]][,j,3,1])^2
cor(crop[,16],S_rain_lpj[[4]][,j,3,1])^2


j=8
cor(crop[,29],S_rain_lpj[[1]][,j,3,1])^2
cor(crop[,30],S_rain_lpj[[2]][,j,3,1])^2
cor(crop[,31],S_rain_lpj[[3]][,j,3,1])^2
cor(crop[,32],S_rain_lpj[[4]][,j,3,1])^2
###############22-05-2014
load("E:/WOLE/validation/S_irri_lpj")#
load("E:\\validation\\irri_WLS/my_result3")
gaga3=array(my_result3,c(186,8,7,4))
gaga2=array(my_result2,c(186,8,7,4))
gaga1=array(my_result1,c(186,8,7,4))
gaga=abind(gaga1,gaga2,gaga3,along=5)
c1=gaga[,,,1,]
c1=aperm(c1,c(1,2,4,3))
c2=gaga[,,,2,]
c2=aperm(c2,c(1,2,4,3))
c3=gaga[,,,3,]
c3=aperm(c3,c(1,2,4,3))
c4=gaga[,,,4,]
c4=aperm(c4,c(1,2,4,3))
irri=list(c1,c2,c3,c4)
##cor(S_irri_lpj[[1]],irri[[1]])^2
cor(S_irri_lpj[[1]][,,3,],irri[[1]])^2
cor(S_irri_lpj[[2]][,,3,],irri[[2]])^2
 cor(S_irri_lpj[[3]][,,3,],irri[[3]])^2
cor(S_irri_lpj[[4]][,,3,],irri[[4]])^2

load("E:\\validation\\irri/my_result1")
gaga=array(my_result1,c(186,8,7,4))
c1=gaga[,,,1]
c2=gaga[,,,2]
c3=gaga[,,,3]
c4=gaga[,,,4]

irri=list(c1,c2,c3,c4)
cor(S_irri_lpj[[1]][,,1,],irri[[1]])^2
cor(S_irri_lpj[[2]][,,1,],irri[[2]])^2
 cor(S_irri_lpj[[3]][,,1,],irri[[3]])^2
cor(S_irri_lpj[[4]][,,1,],irri[[4]])^2


cor(S_irri_lpj[[1]][,,1:3,],irri[[1]][,,,])^2

lpj_irri=S_irri_lpj[[4]][,,1:3,]
emu_irri=irri[[4]][,,,]
cor(emu_irri,lpj_irri)^2
1-(sum((lpj_irri-emu_irri)^2)/sum((lpj_irri-mean(lpj_irri))^2))
#####

###############################################
olu=c(36,176,178,143,29)
olu2=c(which(ind==olu[1]),which(ind==olu[2]),which(ind==olu[3]),which(ind==olu[4]))

par(mfrow=c(2,2))
plot(fitPC2[[22]],main="rice residual from China",type=1)
plot(fitPC2[[90]],main="rice residual from USA",type=1)
plot(fitPC2[[73]],main="rice residual from Russia",type=1)
plot(fitPC2[[13]],main="rice residual from Canada",type=1)



plot(fitPC1[[19]],main="cereal residual from China",type=1)
plot(fitPC1[[79]],main="cereal residual from UK",type=1)
plot(fitPC1[[80]],main="cereal residual from USA",type=1)
plot(fitPC1[[62]],main="cereal residual from Russia",type=1)

plot(fitPC3[[21]],main="maize residual from China",type=1)
plot(fitPC3[[94]],main="maize residual from USA",type=1)
plot(fitPC3[[72]],main="maize residual from Russia",type=1)
plot(fitPC3[[18]],main="maize residual from Canada",type=1)


plot(fitPC4[[15]],main="oil residual from China",type=1)
plot(fitPC4[[62]],main="oil residual from USA",type=1)
plot(fitPC4[[47]],main="oil residual from Russia",type=1)
plot(fitPC4[[13]],main="oil residual from Canada",type=1)

#######################C.I for cereal
par(mfrow=c(2,1))
i=4
cr=crop[,i]
tt=which(cr!=0)
lpj=S_irri_lpj[[i]][,8,3,5]
set.seed(1234)
 df <- data.frame(x =1:length(cr),F =cr,L =cr-SE[[i]],U =cr+SE[[i]],lpj)
df[,2:5]=df[,2:4]*10^(-15)
df=df[tt,]
plot(df$lpj, df$F,main="Emulator prediction with C.I",xlab="UN country levels",ylab="Oil yield (PgC)",ylim=range(df[,2:5]),type="o",lty=1)
#polygon(c(df$x,rev(df$x)),c(df$L,rev(df$U)),col = "blue", border = FALSE)
 lines(df$x, df$U, col="red",lty=1)
 lines(df$x, df$L, col="red",lty=1)
plot(seq(tt), df$lpj,main="LPJmL",xlab="UN country levels",ylab="Oil yield (PgC)",ylim=range(df[,2:5]),type="o",lty=1,cex.axis=1)
#############################

i=4
cr=crop[,i]
tt=which(cr!=0)
lpj=S_irri_lpj[[i]][,8,3,5]
set.seed(1234)
 df <- data.frame(x =1:length(cr),F =cr,L =cr-SE[[i]],U =cr+SE[[i]],lpj)
df[,2:5]=df[,2:4]*10^(-15)
df=df[tt,]
plot(df$lpj, df$F,main="Emulator prediction with C.I",xlab="UN country levels",ylab="Oil yield (PgC)",ylim=range(df[,2:5]),type="o",lty=1)
par(new=TRUE)
plot(df$lpj, df$U, col="red",lty=1,type="l")
par(new=TRUE)
plot(df$lpj, df$L, col="red",lty=1)

#######################OR
library(plotrix)
library(ggplot2)
par(mfrow=c(2,1))
i=4
cr=crop[,i]
tt=which(cr!=0)
lpj=S_irri_lpj[[i]][,8,3,5]
set.seed(1234)
 df <- data.frame(x =1:length(cr),F =cr,L =cr-SE[[i]],U =cr+SE[[i]],lpj)
df[,2:5]=df[,2:4]*10^(-15)
df=df[tt,]
plotCI(df$lpj, df$F,ui=df$U,li=df$L,main="Emulator prediction with C.I",xlab="UN country levels",ylab="Oil yield (PgC)",col="red",err="y")

plotCI(x=seq(tt),yy=df$lpj,main="LPJmL",xlab="UN country levels",ylab="Oil yield (PgC)",ylim=range(df[,2:5]),type="o",lty=1,cex.axis=1)

ggplot(df, aes(x =df$lpj, y = df$F)) +geom_point(size = 4) +geom_errorbar(aes(ymax = df$U, ymin =df$L))
###########################
plotCI(barplot(df$F,col="gray"),df$F,df$L,add=TRUE)
df1=colSums(df1)
df2=colSums(df2)
df3=colSums(df3)
df4=colSums(df4)
par(mfrow=c(2,2))
plotCI(barplot(df1$F,col="gray"),df1$F,df1$L,add=TRUE)
wara=cbind(df1,df2,df3,df4)
############

x <- 1:185
y <- cr*10^(-13)
lo <- loess(y[tt]~x[tt])
plot(x,y)
xl <- seq(min(x),max(x), (max(x) - min(x))/1000)
lines(xl, predict(lo,xl), col='red', lwd=2)
#####################

i=1
cr=crop[,i]
tt=which(cr!=0)
lpj=S_irri_lpj[[i]][,8,3,5]
set.seed(1234)
 df <- data.frame(x =1:length(cr),F =cr,L =cr-SE[[i]],U =cr+SE[[i]],lpj)
df[,2:5]=df[,2:4]*10^(-15)
df1=df[tt,]

i=2
cr=crop[,i]
tt=which(cr!=0)
lpj=S_irri_lpj[[i]][,8,3,5]
set.seed(1234)
 df <- data.frame(x =1:length(cr),F =cr,L =cr-SE[[i]],U =cr+SE[[i]],lpj)
df[,2:5]=df[,2:4]*10^(-15)
df2=df[tt,]

i=3
cr=crop[,i]
tt=which(cr!=0)
lpj=S_irri_lpj[[i]][,8,3,5]
set.seed(1234)
 df <- data.frame(x =1:length(cr),F =cr,L =cr-SE[[i]],U =cr+SE[[i]],lpj)
df[,2:5]=df[,2:4]*10^(-15)
df3=df[tt,]

i=4
cr=crop[,i]
tt=which(cr!=0)
lpj=S_irri_lpj[[i]][,8,3,5]
set.seed(1234)
 df <- data.frame(x =1:length(cr),F =cr,L =cr-SE[[i]],U =cr+SE[[i]],lpj)
df[,2:5]=df[,2:4]*10^(-15)
df4=df[tt,]

df1=cbind(df1,rep("cereal",nrow(df1)))
df2=cbind(df2,rep("rice",nrow(df2)))
df3=cbind(df3,rep("maize",nrow(df3)))
df4=cbind(df4,rep("oil",nrow(df4)))
nam=c("x","F","L","U","lpj","crop")
colnames(df1)=colnames(df2)=colnames(df3)=colnames(df4)=nam
wara=rbind(df1,df2,df3,df4)

#wara2=rbind(colSums(df1),colSums(df2),colSums(df3),colSums(df4))
wara2=rbind(colMeans(df1),colMeans(df2),colMeans(df3),colMeans(df4))
barplot2(wara2[,c(2,5)],ci.l=cbind(wara2[,3],rep(0,4)),ci.u=cbind(wara2[,4],rep(0,4)),plot.ci=TRUE,beside=TRUE,ci.color="red",names.arg=c("Emulator","LPJmL"),col=c("pink","green","blue","grey"),ylab="Yield (PgC)",legend =c("Cereal","Rice","Maize","Oil"))

barplot2(wara2[,c(2,5)],ci.l=cbind(wara2[,3],rep(0,4)),ci.u=cbind(wara2[,4],rep(0,4)),plot.ci=TRUE,beside=TRUE,ci.color="red",names.arg=c("Emulator","LPJmL"),col=c("pink","green","blue","grey"),ylab="Yield (PgC)",legend =c("Cereal","Rice","Maize","Oil"),main="Global average yield change for LPJmL and emulator with C.I")

t1=df1[c(19,79,80,62),]
t2=df2[c(22,90,73,13),]
t3=df3[c(21,94,72,18),]
t4=df4[c(15,62,47,13),]
barplot2(t1[,c(2,5)],ci.l=cbind(t1[,3],rep(0,4)),ci.u=cbind(t1[,4],rep(0,4)),plot.ci=TRUE,beside=TRUE,ci.color="red",names.arg=c("Emulator","LPJmL"),col=c("pink","green","blue","grey"),ylab="Yield (PgC)",legend =c("China","UK","USA","Russia"),main="Change in yield for some selected countries for LPJmL and emulator with C.I")

par(mfrow=c(2,2))
barplot2(as.matrix(t1[,c(2,5)]),ci.l=cbind(t1[,3],rep(0,4)),ci.u=cbind(t1[,4],rep(0,4)),plot.ci=TRUE,beside=TRUE,ci.color="red",names.arg=c("Emulator","LPJmL"),col=c("pink","green","blue","grey"),ylab="Yield (PgC)",legend =c("China","UK","USA","Russia"),main="Change in yield for cereal LPJmL and emulator with C.I")

barplot2(as.matrix(t2[,c(2,5)]),ci.l=cbind(t2[,3],rep(0,4)),ci.u=cbind(t2[,4],rep(0,4)),plot.ci=TRUE,beside=TRUE,ci.color="red",names.arg=c("Emulator","LPJmL"),col=c("pink","green","blue","grey"),ylab="Yield (PgC)",legend =c("China","USA","Russia","Canada"),main="Change in yield for  rice LPJmL and emulator with C.I")

barplot2(as.matrix(t3[,c(2,5)]),ci.l=cbind(t3[,3],rep(0,4)),ci.u=cbind(t3[,4],rep(0,4)),plot.ci=TRUE,beside=TRUE,ci.color="red",names.arg=c("Emulator","LPJmL"),col=c("pink","green","blue","grey"),ylab="Yield (PgC)",legend =c("China","USA","Russia","Canada"),main="Change in yield for  maize LPJmL and emulator with C.I")

barplot2(as.matrix(t4[,c(2,5)]),ci.l=cbind(t4[,3],rep(0,4)),ci.u=cbind(t4[,4],rep(0,4)),plot.ci=TRUE,beside=TRUE,ci.color="red",names.arg=c("Emulator","LPJmL"),col=c("pink","green","blue","grey"),ylab="Yield (PgC)",legend =c("China","USA","Russia","Canada"),main="Change in yield for  maize LPJmL and emulator with C.I")



library(gplots)
library(psych)
library(sciplot)
######################

lpj_rain=S_irri_lpj[[4]][,,,1]
emu_rain=had_on_irri[[4]][,,,1]
cor(emu_rain,lpj_rain)^2
1-(sum((lpj_rain-emu_rain)^2)/sum((lpj_rain-mean(lpj_rain))^2))
#####

lpj_rain=S_irri_lpj[[4]][,,,2]
emu_rain=had_on_irri[[4]][,,,2]
cor(emu_rain,lpj_rain)^2
1-(sum((lpj_rain-emu_rain)^2)/sum((lpj_rain-mean(lpj_rain))^2))
#####
lpj_rain=S_irri_lpj[[4]][,,,3]
emu_rain=had_on_irri[[4]][,,,3]
cor(emu_rain,lpj_rain)^2
1-(sum((lpj_rain-emu_rain)^2)/sum((lpj_rain-mean(lpj_rain))^2))
#####
lpj_rain=S_irri_lpj[[4]][,,,4]
emu_rain=had_on_irri[[4]][,,,4]
cor(emu_rain,lpj_rain)^2
1-(sum((lpj_rain-emu_rain)^2)/sum((lpj_rain-mean(lpj_rain))^2))
#####
lpj_rain=S_irri_lpj[[4]][,,,5]
emu_rain=had_on_irri[[4]][,,,5]
cor(emu_rain,lpj_rain)^2
1-(sum((lpj_rain-emu_rain)^2)/sum((lpj_rain-mean(lpj_rain))^2))
#####
lpj_rain=S_irri_lpj[[4]][,,,6]
emu_rain=had_on_irri[[4]][,,,6]
cor(emu_rain,lpj_rain)^2
1-(sum((lpj_rain-emu_rain)^2)/sum((lpj_rain-mean(lpj_rain))^2))
#####
lpj_rain=S_irri_lpj[[4]][,,,7]
emu_rain=had_on_irri[[4]][,,,7]
cor(emu_rain,lpj_rain)^2
1-(sum((lpj_rain-emu_rain)^2)/sum((lpj_rain-mean(lpj_rain))^2))
#####