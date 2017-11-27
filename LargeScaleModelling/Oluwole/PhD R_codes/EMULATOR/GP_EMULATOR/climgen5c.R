dirname <- getwd()
################################
if(j==0) {

j=1
source("irrigated2/input_irri.R",echo=TRUE)
if(fertilization==0){
fert=3
} else {
if(fertilization==1) {
fert=1
} else {
fert=2
}}
npixel=59199
load("irrigated2/k3")
library(ncdf)
load("irrigated2/ini")
load("irrigated2/crop_N")
load("irrigated2/grid_out")
load("irrigated2/crop_form4b")
load("irrigated2/coef4b")
m=i
###################################################
nam=c("cereal","rice","maize","oil","ini_cereal","ini_rice","ini_maize","ini_oil","scld" ,   "wcld"  ,  "spcld" ,  "acld" ,"spre",   
 "wpre" ,   "sppre" ,  "apre" ,   "stmp" ,   "wtmp" ,   "sptmp"  , "atmp",   
"swet" ,   "wwet" ,   "spwet" ,  "awet" ,   "iscld", "iwcld",   "ispcld", 
"iacld","ispre" ,"iwpre" ,  "isppre" , "iapre" ,  "istmp" ,  "iwtmp" , 
"isptmp", "iatmp" ,  "iswet"  , "iwwet" ,  "ispwet" , "iawet" ,  "cco2",   
"co2" , "soil",  "lat.lpj","lon.lpj","LA") 
################################
w=1
s2a=s2[,1,]
input=s2a
LA=rep(man,npixel)
input=cbind(input,LA)
colnames(input)=nam
input=as.data.frame(input)

if(RC==5){
ini_cereal=input[,5]=rowMeans(ini[[1]][,1:4,man,fert])
ini_rice=input[,6]=rowMeans(ini[[2]][,1:4,man,fert])
ini_maize=input[,7]=rowMeans(ini[[3]][,1:4,man,fert])
ini_oil=input[,8]=rowMeans(ini[[4]][,1:4,man,fert])
} else {
ini_cereal=input[,5]=ini[[1]][,RC,man,fert]
ini_rice=input[,6]=ini[[2]][,RC,man,fert]
ini_maize=input[,7]=ini[[3]][,RC,man,fert]
ini_oil=input[,8]=ini[[4]][,RC,man,fert]
}
ini_crop=list(ini_cereal,ini_rice,ini_maize,ini_oil)


if(fertilization==1){
 co2=input[,42]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,42]+rep(370.47,npixel)/2))
}}

if(fertilization==1){
 cco2=input[,41]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,41]+rep(0,npixel)/2))
}}

input[,41]=cco2
input[,42]=co2
#########cereal with other crops
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
prediction<-{foreach(m=1:4,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4b[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4b[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
################
#####include PCA analysis
source("irrigated2/wls_irri2.R",echo=F)

library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction2[,m]+E_crop[[m]]
crop=matrix(0,nrow=186,ncol=1)
crop[m]=tj
}}
stopCluster(cl)
m=i
#ini1=input[,(4+m)]
crop1=crpp
p1=prediction

################################
j=2
source("irrigated2/input_irri.R",echo=TRUE)
w=2
s2a=s2[,1,]
input=s2a
LA=rep(man,npixel)
input=cbind(input,LA)
colnames(input)=nam
input=as.data.frame(input)
if(RC==5){
ini_cereal=input[,5]=rowMeans(ini[[1]][,1:4,man,fert])
ini_rice=input[,6]=rowMeans(ini[[2]][,1:4,man,fert])
ini_maize=input[,7]=rowMeans(ini[[3]][,1:4,man,fert])
ini_oil=input[,8]=rowMeans(ini[[4]][,1:4,man,fert])
} else {
ini_cereal=input[,5]=ini[[1]][,RC,man,fert]
ini_rice=input[,6]=ini[[2]][,RC,man,fert]
ini_maize=input[,7]=ini[[3]][,RC,man,fert]
ini_oil=input[,8]=ini[[4]][,RC,man,fert]
}
ini_crop=list(ini_cereal,ini_rice,ini_maize,ini_oil)


if(fertilization==1){
 co2=input[,42]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,42]+rep(370.47,npixel)/2))
}}

if(fertilization==1){
 cco2=input[,41]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,41]+rep(0,npixel)/2))}
}

input[,41]=cco2
input[,42]=co2
#########cereal with other crops
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
prediction<-{foreach(m=1:4,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4b[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4b[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
################
#####include PCA analysis
source("irrigated2/wls_irri2.R",echo=F)
library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction2[,m]+E_crop[[m]]
crop=matrix(0,nrow=186,ncol=1)
crop[m]=tj
}}
stopCluster(cl)
m=i
#ini2=input[,(4+m)]
crop2=crpp
p2=prediction

################################
j=3
source("irrigated2/input_irri.R",echo=TRUE)

w=3
s2a=s2[,1,]
input=s2a
LA=rep(man,npixel)
input=cbind(input,LA)
colnames(input)=nam
input=as.data.frame(input)
if(RC==5){
ini_cereal=input[,5]=rowMeans(ini[[1]][,1:4,man,fert])
ini_rice=input[,6]=rowMeans(ini[[2]][,1:4,man,fert])
ini_maize=input[,7]=rowMeans(ini[[3]][,1:4,man,fert])
ini_oil=input[,8]=rowMeans(ini[[4]][,1:4,man,fert])
} else {
ini_cereal=input[,5]=ini[[1]][,RC,man,fert]
ini_rice=input[,6]=ini[[2]][,RC,man,fert]
ini_maize=input[,7]=ini[[3]][,RC,man,fert]
ini_oil=input[,8]=ini[[4]][,RC,man,fert]
}
ini_crop=list(ini_cereal,ini_rice,ini_maize,ini_oil)


if(fertilization==1){
 co2=input[,42]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,42]+rep(370.47,npixel)/2))
}}

if(fertilization==1){
 cco2=input[,41]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,41]+rep(0,npixel)/2))}
}

input[,41]=cco2
input[,42]=co2
#########cereal with other crops
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
prediction<-{foreach(m=1:4,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4b[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4b[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
################
#####include PCA analysis
source("irrigated2/wls_irri2.R",echo=F)
library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction2[,m]+E_crop[[m]]
crop=matrix(0,nrow=186,ncol=1)
crop[m]=tj
}}
stopCluster(cl)
m=i
#ini3=input[,(4+m)]
crop3=crpp
p3=prediction

################################
j=4
source("irrigated2/input_irri.R",echo=TRUE)

w=4
s2a=s2[,1,]
input=s2a
LA=rep(man,npixel)
input=cbind(input,LA)
colnames(input)=nam
input=as.data.frame(input)
if(RC==5){
ini_cereal=input[,5]=rowMeans(ini[[1]][,1:4,man,fert])
ini_rice=input[,6]=rowMeans(ini[[2]][,1:4,man,fert])
ini_maize=input[,7]=rowMeans(ini[[3]][,1:4,man,fert])
ini_oil=input[,8]=rowMeans(ini[[4]][,1:4,man,fert])
} else {
ini_cereal=input[,5]=ini[[1]][,RC,man,fert]
ini_rice=input[,6]=ini[[2]][,RC,man,fert]
ini_maize=input[,7]=ini[[3]][,RC,man,fert]
ini_oil=input[,8]=ini[[4]][,RC,man,fert]
}
ini_crop=list(ini_cereal,ini_rice,ini_maize,ini_oil)


if(fertilization==1){
 co2=input[,42]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,42]+rep(370.47,npixel)/2))}
}

if(fertilization==1){
 cco2=input[,41]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,41]+rep(0,npixel)/2))}
}

input[,41]=cco2
input[,42]=co2
#########cereal with other crops
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
prediction<-{foreach(m=1:4,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4b[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4b[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
################
#####include PCA analysis
source("irrigated2/wls_irri2.R",echo=F)
library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction2[,m]+E_crop[[m]]
crop=matrix(0,nrow=186,ncol=1)
crop[m]=tj
}}
stopCluster(cl)
m=i
#ini4=input[,(4+m)]
crop4=crpp
p4=prediction

################################
j=5
source("irrigated2/input_irri.R",echo=TRUE)

w=5
s2a=s2[,1,]
input=s2a
LA=rep(man,npixel)
input=cbind(input,LA)
colnames(input)=nam
input=as.data.frame(input)
if(RC==5){
ini_cereal=input[,5]=rowMeans(ini[[1]][,1:4,man,fert])
ini_rice=input[,6]=rowMeans(ini[[2]][,1:4,man,fert])
ini_maize=input[,7]=rowMeans(ini[[3]][,1:4,man,fert])
ini_oil=input[,8]=rowMeans(ini[[4]][,1:4,man,fert])
} else {
ini_cereal=input[,5]=ini[[1]][,RC,man,fert]
ini_rice=input[,6]=ini[[2]][,RC,man,fert]
ini_maize=input[,7]=ini[[3]][,RC,man,fert]
ini_oil=input[,8]=ini[[4]][,RC,man,fert]
}
ini_crop=list(ini_cereal,ini_rice,ini_maize,ini_oil)

if(fertilization==1){
 co2=input[,42]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,42]+rep(370.47,npixel)/2))}
}

if(fertilization==1){
 cco2=input[,41]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,41]+rep(0,npixel)/2))}
}

input[,41]=cco2
input[,42]=co2
#########cereal with other crops
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
prediction<-{foreach(m=1:4,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4b[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4b[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
################
#####include PCA analysis
source("irrigated2/wls_irri2.R",echo=F)
library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction2[,m]+E_crop[[m]]
crop=matrix(0,nrow=186,ncol=1)
crop[m]=tj
}}
stopCluster(cl)
m=i
#ini5=input[,(4+m)]
crop5=crpp
p5=prediction

################################
j=6
source("irrigated2/input_irri.R",echo=TRUE)

w=6
s2a=s2[,1,]
input=s2a
LA=rep(man,npixel)
input=cbind(input,LA)
colnames(input)=nam
input=as.data.frame(input)
if(RC==5){
ini_cereal=input[,5]=rowMeans(ini[[1]][,1:4,man,fert])
ini_rice=input[,6]=rowMeans(ini[[2]][,1:4,man,fert])
ini_maize=input[,7]=rowMeans(ini[[3]][,1:4,man,fert])
ini_oil=input[,8]=rowMeans(ini[[4]][,1:4,man,fert])
} else {
ini_cereal=input[,5]=ini[[1]][,RC,man,fert]
ini_rice=input[,6]=ini[[2]][,RC,man,fert]
ini_maize=input[,7]=ini[[3]][,RC,man,fert]
ini_oil=input[,8]=ini[[4]][,RC,man,fert]
}
ini_crop=list(ini_cereal,ini_rice,ini_maize,ini_oil)


if(fertilization==1){
 co2=input[,42]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,42]+rep(370.47,npixel)/2))}
}

if(fertilization==1){
 cco2=input[,41]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,41]+rep(0,npixel)/2))}
}

input[,41]=cco2
input[,42]=co2
#########cereal with other crops
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
prediction<-{foreach(m=1:4,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4b[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4b[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
################
#####include PCA analysis
source("irrigated2/wls_irri2.R",echo=F)
library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction2[,m]+E_crop[[m]]
crop=matrix(0,nrow=186,ncol=1)
crop[m]=tj
}}
stopCluster(cl)
m=i
#ini6=input[,(4+m)]
crop6=crpp
p6=prediction

################################
j=7
source("irrigated2/input_irri.R",echo=TRUE)

w=7
s2a=s2[,1,]
input=s2a
LA=rep(man,npixel)
input=cbind(input,LA)
colnames(input)=nam
input=as.data.frame(input)
if(RC==5){
ini_cereal=input[,5]=rowMeans(ini[[1]][,1:4,man,fert])
ini_rice=input[,6]=rowMeans(ini[[2]][,1:4,man,fert])
ini_maize=input[,7]=rowMeans(ini[[3]][,1:4,man,fert])
ini_oil=input[,8]=rowMeans(ini[[4]][,1:4,man,fert])
} else {
ini_cereal=input[,5]=ini[[1]][,RC,man,fert]
ini_rice=input[,6]=ini[[2]][,RC,man,fert]
ini_maize=input[,7]=ini[[3]][,RC,man,fert]
ini_oil=input[,8]=ini[[4]][,RC,man,fert]
}
ini_crop=list(ini_cereal,ini_rice,ini_maize,ini_oil)


if(fertilization==1){
 co2=input[,42]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,42]+rep(370.47,npixel)/2))}
}

if(fertilization==1){
 cco2=input[,41]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,41]+rep(0,npixel)/2))}
}

input[,41]=cco2
input[,42]=co2
#########cereal with other crops
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
prediction<-{foreach(m=1:4,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4b[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4b[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
################
#####include PCA analysis
source("irrigated2/wls_irri2.R",echo=F)
library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction2[,m]+E_crop[[m]]
crop=matrix(0,nrow=186,ncol=1)
crop[m]=tj
}}
stopCluster(cl)
m=i
#ini7=input[,(4+m)]
crop7=crpp
p7=prediction

################################
j=8
source("irrigated2/input_irri.R",echo=TRUE)

w=8
s2a=s2[,1,]
input=s2a
LA=rep(man,npixel)
input=cbind(input,LA)
colnames(input)=nam
input=as.data.frame(input)
if(RC==5){
ini_cereal=input[,5]=rowMeans(ini[[1]][,1:4,man,fert])
ini_rice=input[,6]=rowMeans(ini[[2]][,1:4,man,fert])
ini_maize=input[,7]=rowMeans(ini[[3]][,1:4,man,fert])
ini_oil=input[,8]=rowMeans(ini[[4]][,1:4,man,fert])
} else {
ini_cereal=input[,5]=ini[[1]][,RC,man,fert]
ini_rice=input[,6]=ini[[2]][,RC,man,fert]
ini_maize=input[,7]=ini[[3]][,RC,man,fert]
ini_oil=input[,8]=ini[[4]][,RC,man,fert]
}
ini_crop=list(ini_cereal,ini_rice,ini_maize,ini_oil)


if(fertilization==1){
 co2=input[,42]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,42]+rep(370.47,npixel)/2))}
}

if(fertilization==1){
 cco2=input[,41]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,41]+rep(0,npixel)/2))}
}

input[,41]=cco2
input[,42]=co2
ini_crop=list(ini_cereal,ini_rice,ini_maize,ini_oil)

##########################
#########cereal with other crops
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
prediction<-{foreach(m=1:4,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4b[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4b[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
################
#####include PCA analysis
source("irrigated2/wls_irri2.R",echo=F)
library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction2[,m]+E_crop[[m]]
crop=matrix(0,nrow=186,ncol=1)
crop[m]=tj
}}
stopCluster(cl)
m=i
ini8=input[,(4+m)]
crop8=crpp
p8=prediction

#################
predi=rbind(p1,p2,p3,p4,p5,p6,p7,p8)
my_pred=predi

crop=rbind(crop1,crop2,crop3,crop4,crop5,crop6,crop7,crop8)
my_crop=crop
#ini=ini8
#colnames(crop1)=colnames(crop2)=colnames(crop3)=colnames(crop4)=colnames(crop5)=colnames(crop6)=colnames(crop7)=colnames(crop8)=c("cereal","rice","maize","oil_max")
#crop=cbind(crop1,ini1,crop2,ini2,crop3,ini3,crop4,ini4,crop5,ini5,crop6,ini6,crop7,ini7,crop8,ini8)
#crop=as.matrix(crop)
#j=0
#source("irrigated2/aggregate3b.R",echo=F)
#colnames(ini)=c("ini_cereal","ini_rice","ini_maize","ini_oil")
#result_model=cbind(crop,ini)
#write.csv(result_model,file=paste(dirname,sub("v3",v3,(sub("j",1,(sub("j",j,"crop_ALL_8_v3_result.csv"))))),sep="/"),row.names=FALSE)


###############################END
} else {
if(fertilization==0){
fert=3
} else {
if(fertilization==1) {
fert=1
} else {
fert=2
}}
j=j
source("irrigated2/input_irri.R",echo=TRUE)

npixel=59199
#load("irrigated2/k3")
library(ncdf)
load("irrigated2/ini")
load("irrigated2/crop_N")
load("irrigated2/grid_out")
load("irrigated2/crop_form4b")
load("irrigated2/coef4b")
m=i
###################################################
nam=c("cereal","rice","maize","oil","ini_cereal","ini_rice","ini_maize","ini_oil","scld" ,   "wcld"  ,  "spcld" ,  "acld" ,"spre",   
 "wpre" ,   "sppre" ,  "apre" ,   "stmp" ,   "wtmp" ,   "sptmp"  , "atmp",   
"swet" ,   "wwet" ,   "spwet" ,  "awet" ,   "iscld", "iwcld",   "ispcld", 
"iacld","ispre" ,"iwpre" ,  "isppre" , "iapre" ,  "istmp" ,  "iwtmp" , 
"isptmp", "iatmp" ,  "iswet"  , "iwwet" ,  "ispwet" , "iawet" ,  "cco2",   
"co2" , "soil",  "lat.lpj","lon.lpj","LA") 
w=1
s2a=s2[,1,]
input=s2a
LA=rep(man,npixel)
input=cbind(input,LA)
colnames(input)=nam
input=as.data.frame(input)
if(RC==5){
ini_cereal=input[,5]=rowMeans(ini[[1]][,1:4,man,fert])
ini_rice=input[,6]=rowMeans(ini[[2]][,1:4,man,fert])
ini_maize=input[,7]=rowMeans(ini[[3]][,1:4,man,fert])
ini_oil=input[,8]=rowMeans(ini[[4]][,1:4,man,fert])
} else {
ini_cereal=input[,5]=ini[[1]][,RC,man,fert]
ini_rice=input[,6]=ini[[2]][,RC,man,fert]
ini_maize=input[,7]=ini[[3]][,RC,man,fert]
ini_oil=input[,8]=ini[[4]][,RC,man,fert]
}
ini_crop=list(ini_cereal,ini_rice,ini_maize,ini_oil)


if(fertilization==1){
 co2=input[,42]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,42]+rep(370.47,npixel)/2))
}}

if(fertilization==1){
 cco2=input[,41]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,41]+rep(0,npixel)/2))
}}

input[,41]=cco2
input[,42]=co2

################cereal with other crops
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
prediction<-{foreach(m=1:4,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4b[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4b[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)


#################################################
#####include PCA analysis
source("irrigated2/wls_irri2.R",echo=F)

######
library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction2[,m]+E_crop[[m]]
crop=matrix(0,nrow=186,ncol=1)
crop[m]=tj
}}
stopCluster(cl)
m=i
ini=as.matrix(input[,(4+m)])

source("irrigated2/aggregate3b.R",echo=F)
crop=cbind(crpp,ini)
crop=as.matrix(crop)
colnames(crop)=c("cereal","rice","maize","oil_max","ini_cereal","ini_rice","ini_maize","ini_oil")

write.csv(crop,file=paste(dirname,sub("v3",v3,(sub("j",1,(sub("j",j,"crop_j_j_v3_result.csv"))))),sep="/"),row.names=FALSE)
}
