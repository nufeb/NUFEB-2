#load("E:\\WOLE\\GP/my2")
#load("E:\\WOLE\\GP/myp1")
#load("E:\\WOLE\\GP/myp3")
#load("E:\\WOLE\\GP/myp4")
#mm1=array(myp1,c(59199,8,7,4,1))###grid X decade X man X crop X rcp1###
#mm2=array(myp2,c(59199,8,7,4,1))
#mm3=array(myp3,c(59199,8,7,4,1))
#mm4=array(myp4,c(59199,8,7,4,1))
library(abind)
#mm5=abind(mm1,mm2,mm3,mm4,along=5)
load("mm5")
#prediction=mm5[,j,man,,RC]
####################
###load prediction
#source(".R")
j=1
prediction=mm5[,j,man,,RC]
source("rainfed2/wls_rain2.R",echo=F)

library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction2[,m]+E_crop[[m]]
}}
stopCluster(cl)
m=i
crop1=crpp
####
j=2
prediction=mm5[,j,man,,RC]
source("rainfed2/wls_rain2.R",echo=F)

library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction2[,m]+E_crop[[m]]
}}
stopCluster(cl)
m=i
crop2=crpp
###
j=3
prediction=mm5[,j,man,,RC]
source("rainfed2/wls_rain2.R",echo=F)

library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction2[,m]+E_crop[[m]]
}}
stopCluster(cl)
m=i
crop3=crpp
###
j=4
prediction=mm5[,j,man,,RC]
source("rainfed2/wls_rain2.R",echo=F)

library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction2[,m]+E_crop[[m]]
}}
stopCluster(cl)
m=i
crop4=crpp
###
j=5
prediction=mm5[,j,man,,RC]
source("rainfed2/wls_rain2.R",echo=F)

library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction2[,m]+E_crop[[m]]
}}
stopCluster(cl)
m=i
crop5=crpp
###
j=6
prediction=mm5[,j,man,,RC]
source("rainfed2/wls_rain2.R",echo=F)

library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction2[,m]+E_crop[[m]]
}}
stopCluster(cl)
m=i
crop6=crpp
##
j=7
prediction=mm5[,j,man,,RC]
source("rainfed2/wls_rain2.R",echo=F)

library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction2[,m]+E_crop[[m]]
}}
stopCluster(cl)
m=i
crop7=crpp
##
j=8
prediction=mm5[,j,man,,RC]
source("rainfed2/wls_rain2.R",echo=F)

library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction2[,m]+E_crop[[m]]
}}
stopCluster(cl)
m=i
crop8=crpp
crop=rbind(crop1,crop2,crop3,crop4,crop5,crop6,crop7,crop8)
my_crop=crop
#######################
