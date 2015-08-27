###############correction to residuals irrigated 24-05-214
load("E:\\WOLE\\validation\\irri\\residual/S3")
load("irrigated2/masked/crop_err4c")

V=aperm(S3,c(1,2,3,5,4))
V2=array(V,c(186,8,16,7))
W=crop_err4c[[3]]
olu=W[,,1:4,]
library(abind)
V=abind(olu,V2,along=3)
crop_err4c[[3]]=V
dirname="E:\\WOLE\\GP/irrigated2/masked"
save(crop_err4c,file=paste(dirname,"crop_err4c",sep="/"))
rm(list=ls())
########################
load("E:\\WOLE\\validation\\irri\\residual/S3_u")
load("irrigated2/unmasked/crop_err4c")

V=aperm(S3_u,c(1,2,3,5,4))
V2=array(V,c(186,8,16,7))
W=crop_err4c[[3]]
olu=W[,,1:4,]
library(abind)
V=abind(olu,V2,along=3)
crop_err4c[[3]]=V
dirname="E:\\WOLE\\GP/irrigated2/unmasked"
save(crop_err4c,file=paste(dirname,"crop_err4c",sep="/"))
rm(list=ls())
###########################masked+rainfed
setwd("E:\\WOLE\\validation")
load("rain/residual/residual_rain_gfdl")####also have problem
load("rain/residual/residual_miroc")
load("E:\\GPemulator\\masked\\rainfed2\\residual/S3")

gaga1=aperm(residual_rain_gfdl[[1]],c(1,2,3,5,4))
gaga2=aperm(residual_rain_gfdl[[2]],c(1,2,3,5,4))
#gaga3=aperm(residual_rain_gfdl[[3]],c(1,2,3,5,4))
gaga3=aperm(S3,c(1,2,3,5,4))
gaga4=aperm(residual_rain_gfdl[[4]],c(1,2,3,5,4))

ga1=array(gaga1,c(186,8,4,7,1))
ga2=array(gaga2,c(186,8,4,7,1))
ga3=array(gaga3,c(186,8,4,7,1))
ga4=array(gaga4,c(186,8,4,7,1))
library(abind)###186   8   4   7   4
err1=abind(ga1,residual_miroc[[1]],along=5)
err2=abind(ga2,residual_miroc[[2]],along=5)
err3=abind(ga3,residual_miroc[[3]],along=5)
err4=abind(ga4,residual_miroc[[4]],along=5)

err1=aperm(err1,c(1,2,3,5,4))
err2=aperm(err2,c(1,2,3,5,4))
err3=aperm(err3,c(1,2,3,5,4))
err4=aperm(err4,c(1,2,3,5,4))

r1=array(err1,c(186,8,20,7))
r2=array(err2,c(186,8,20,7))
r3=array(err3,c(186,8,20,7))
r4=array(err4,c(186,8,20,7))
crop_err4c=list(r1,r2,r3,r4)
dirname="E:\\WOLE\\GP\\rainfed2\\masked"
save(crop_err4c,file=paste(dirname,"crop_err4c",sep="/"))
###############################unmasked+rainfed
load("E:\\GPemulator\\unmasked\\rainfed2\\residual/S3_residual")
load("rain/residual/residual_rain_gfdl_un")
load("rain/residual/residual_miroc_u")
gaga1=aperm(residual_rain_gfdl[[1]],c(1,2,3,5,4))
gaga2=aperm(residual_rain_gfdl[[2]],c(1,2,3,5,4))
#gaga3=aperm(residual_rain_gfdl[[3]],c(1,2,3,5,4))
gaga3=aperm(S3_residual,c(1,2,3,5,4))

gaga4=aperm(residual_rain_gfdl[[4]],c(1,2,3,5,4))

ga1=array(gaga1,c(186,8,4,7,1))
ga2=array(gaga2,c(186,8,4,7,1))
ga3=array(gaga3,c(186,8,4,7,1))
ga4=array(gaga4,c(186,8,4,7,1))
library(abind)###186   8   4   7   4
err1=abind(ga1,residual_miroc_u[[1]],along=5)
err2=abind(ga2,residual_miroc_u[[2]],along=5)
err3=abind(ga3,residual_miroc_u[[3]],along=5)
err4=abind(ga4,residual_miroc_u[[4]],along=5)

err1=aperm(err1,c(1,2,3,5,4))
err2=aperm(err2,c(1,2,3,5,4))
err3=aperm(err3,c(1,2,3,5,4))
err4=aperm(err4,c(1,2,3,5,4))

r1=array(err1,c(186,8,20,7))
r2=array(err2,c(186,8,20,7))
r3=array(err3,c(186,8,20,7))
r4=array(err4,c(186,8,20,7))
crop_err4c=list(r1,r2,r3,r4)
dirname="E:\\WOLE\\GP\\rainfed2\\unmasked"
save(crop_err4c,file=paste(dirname,"crop_err4c",sep="/"))
#####################corrected hadgem lpj rain
load("E:/WOLE/validation/S_rain_lpj")#
load("E:\\GPemulator\\masked\\rainfed2\\residual/S3_lpj_had")
S_rain_lpj[[3]]=S3_lpj_had
dirname="E:/WOLE/validation"
save(S_rain_lpj,file=paste(dirname,"S_rain_lpj",sep="/"))

