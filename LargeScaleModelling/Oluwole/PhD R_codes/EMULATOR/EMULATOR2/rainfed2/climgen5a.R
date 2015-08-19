###############################
if(j==0) {
j=1
source("rainfed2/input_rain.R",echo=TRUE)
if(fertilization==0){
fert=3
} else {
if(fertilization==1) {
fert=1
} else {
fert=2
}}
npixel=59199
load("rainfed2/k3")
library(ncdf)
load("rainfed2/ini")
load("rainfed2/crop_N")
load("rainfed2/grid_out")
load("rainfed2/crop_form4b")
load("rainfed2/coef4b")
m=i
###################################################
nam=c("ini_tropical_cereal","ini_pulses","ini_temperate_root","ini_soybeans","scld" ,   "wcld"  ,  "spcld" ,  "acld" ,"spre",   
 "wpre" ,   "sppre" ,  "apre" ,   "stmp" ,   "wtmp" ,   "sptmp"  , "atmp",   
"swet" ,   "wwet" ,   "spwet" ,  "awet" ,   "iscld", "iwcld",   "ispcld", 
"iacld","ispre" ,"iwpre" ,  "isppre" , "iapre" ,  "istmp" ,  "iwtmp" , 
"isptmp", "iatmp" ,  "iswet"  , "iwwet" ,  "ispwet" , "iawet" ,  "cco2", "co2" , "soil",  "lat.lpj","LA") 
################################
w=1
s2a=s2[,1,]
input=s2a
LA=rep(man,npixel)
input=cbind(input,LA)
colnames(input)=nam
input=as.data.frame(input)


if(RC==5){
ini_tropical_cereal=input[,1]=rowMeans(ini[[1]][,1:4,man,fert])
ini_pulses=input[,2]=rowMeans(ini[[2]][,1:4,man,fert])
ini_temperate_root=input[,3]=rowMeans(ini[[3]][,1:4,man,fert])
ini_soybeans=input[,4]=rowMeans(ini[[4]][,1:4,man,fert])
} else {
ini_tropical_cereal=input[,1]=ini[[1]][,RC,man,fert]
ini_pulses=input[,2]=ini[[2]][,RC,man,fert]
ini_temperate_root=input[,3]=ini[[3]][,RC,man,fert]
ini_soybeans=input[,4]=ini[[4]][,RC,man,fert]
}
ini_crop=list(ini_tropical_cereal,ini_pulses,ini_temperate_root,ini_soybeans)


if(fertilization==1){
 co2=input[,38]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,38]+rep(370.47,npixel)/2))
}}

if(fertilization==1){
 cco2=input[,37]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,37]+rep(0,npixel)/2))
}}

input[,37]=cco2
input[,38]=co2
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
source("rainfed2/wls_rain.R",echo=F)

library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction[crop_N[[m]],m]+E_crop[[m]]
crop=matrix(0,nrow=npixel,ncol=1)
crop[crop_N[[m]]]=tj
crop
}}
stopCluster(cl)
m=i
ini1=input[,(m)]
crop1=crpp
################################
j=2
source("rainfed2/input_rain.R",echo=TRUE)
w=2
s2a=s2[,1,]
input=s2a
LA=rep(man,npixel)
input=cbind(input,LA)
colnames(input)=nam
input=as.data.frame(input)
if(RC==5){
ini_tropical_cereal=input[,1]=rowMeans(ini[[1]][,1:4,man,fert])
ini_pulses=input[,2]=rowMeans(ini[[2]][,1:4,man,fert])
ini_temperate_root=input[,3]=rowMeans(ini[[3]][,1:4,man,fert])
ini_soybeans=input[,4]=rowMeans(ini[[4]][,1:4,man,fert])
} else {
ini_tropical_cereal=input[,1]=ini[[1]][,RC,man,fert]
ini_pulses=input[,2]=ini[[2]][,RC,man,fert]
ini_temperate_root=input[,3]=ini[[3]][,RC,man,fert]
ini_soybeans=input[,4]=ini[[4]][,RC,man,fert]
}
ini_crop=list(ini_tropical_cereal,ini_pulses,ini_temperate_root,ini_soybeans)


if(fertilization==1){
 co2=input[,38]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,38]+rep(370.47,npixel)/2))
}}

if(fertilization==1){
 cco2=input[,37]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,37]+rep(0,npixel)/2))}
}

input[,37]=cco2
input[,38]=co2
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
source("rainfed2/wls_rain.R",echo=F)
library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction[crop_N[[m]],m]+E_crop[[m]]
crop=matrix(0,nrow=npixel,ncol=1)
crop[crop_N[[m]]]=tj
crop
}}
stopCluster(cl)
m=i
ini2=input[,(m)]
crop2=crpp
################################
j=3
source("rainfed2/input_rain.R",echo=TRUE)

w=3
s2a=s2[,1,]
input=s2a
LA=rep(man,npixel)
input=cbind(input,LA)
colnames(input)=nam
input=as.data.frame(input)
if(RC==5){
ini_tropical_cereal=input[,1]=rowMeans(ini[[1]][,1:4,man,fert])
ini_pulses=input[,2]=rowMeans(ini[[2]][,1:4,man,fert])
ini_temperate_root=input[,3]=rowMeans(ini[[3]][,1:4,man,fert])
ini_soybeans=input[,4]=rowMeans(ini[[4]][,1:4,man,fert])
} else {
ini_tropical_cereal=input[,1]=ini[[1]][,RC,man,fert]
ini_pulses=input[,2]=ini[[2]][,RC,man,fert]
ini_temperate_root=input[,3]=ini[[3]][,RC,man,fert]
ini_soybeans=input[,4]=ini[[4]][,RC,man,fert]
}
ini_crop=list(ini_tropical_cereal,ini_pulses,ini_temperate_root,ini_soybeans)


if(fertilization==1){
 co2=input[,38]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,38]+rep(370.47,npixel)/2))
}}

if(fertilization==1){
 cco2=input[,37]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,37]+rep(0,npixel)/2))}
}

input[,37]=cco2
input[,38]=co2
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
source("rainfed2/wls_rain.R",echo=F)
library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction[crop_N[[m]],m]+E_crop[[m]]
crop=matrix(0,nrow=npixel,ncol=1)
crop[crop_N[[m]]]=tj
crop
}}
stopCluster(cl)
m=i
ini3=input[,(m)]
crop3=crpp
################################
j=4
source("rainfed2/input_rain.R",echo=TRUE)

w=4
s2a=s2[,1,]
input=s2a
LA=rep(man,npixel)
input=cbind(input,LA)
colnames(input)=nam
input=as.data.frame(input)
if(RC==5){
ini_tropical_cereal=input[,1]=rowMeans(ini[[1]][,1:4,man,fert])
ini_pulses=input[,2]=rowMeans(ini[[2]][,1:4,man,fert])
ini_temperate_root=input[,3]=rowMeans(ini[[3]][,1:4,man,fert])
ini_soybeans=input[,4]=rowMeans(ini[[4]][,1:4,man,fert])
} else {
ini_tropical_cereal=input[,1]=ini[[1]][,RC,man,fert]
ini_pulses=input[,2]=ini[[2]][,RC,man,fert]
ini_temperate_root=input[,3]=ini[[3]][,RC,man,fert]
ini_soybeans=input[,4]=ini[[4]][,RC,man,fert]
}
ini_crop=list(ini_tropical_cereal,ini_pulses,ini_temperate_root,ini_soybeans)


if(fertilization==1){
 co2=input[,38]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,38]+rep(370.47,npixel)/2))}
}

if(fertilization==1){
 cco2=input[,37]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,37]+rep(0,npixel)/2))}
}

input[,37]=cco2
input[,38]=co2
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
source("rainfed2/wls_rain.R",echo=F)
library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction[crop_N[[m]],m]+E_crop[[m]]
crop=matrix(0,nrow=npixel,ncol=1)
crop[crop_N[[m]]]=tj
crop
}}
stopCluster(cl)
m=i
ini4=input[,(m)]
crop4=crpp
################################
j=5
source("rainfed2/input_rain.R",echo=TRUE)

w=5
s2a=s2[,1,]
input=s2a
LA=rep(man,npixel)
input=cbind(input,LA)
colnames(input)=nam
input=as.data.frame(input)
if(RC==5){
ini_tropical_cereal=input[,1]=rowMeans(ini[[1]][,1:4,man,fert])
ini_pulses=input[,2]=rowMeans(ini[[2]][,1:4,man,fert])
ini_temperate_root=input[,3]=rowMeans(ini[[3]][,1:4,man,fert])
ini_soybeans=input[,4]=rowMeans(ini[[4]][,1:4,man,fert])
} else {
ini_tropical_cereal=input[,1]=ini[[1]][,RC,man,fert]
ini_pulses=input[,2]=ini[[2]][,RC,man,fert]
ini_temperate_root=input[,3]=ini[[3]][,RC,man,fert]
ini_soybeans=input[,4]=ini[[4]][,RC,man,fert]
}
ini_crop=list(ini_tropical_cereal,ini_pulses,ini_temperate_root,ini_soybeans)


if(fertilization==1){
 co2=input[,38]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,38]+rep(370.47,npixel)/2))}
}

if(fertilization==1){
 cco2=input[,37]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,37]+rep(0,npixel)/2))}
}

input[,37]=cco2
input[,38]=co2
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
source("rainfed2/wls_rain.R",echo=F)
library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction[crop_N[[m]],m]+E_crop[[m]]
crop=matrix(0,nrow=npixel,ncol=1)
crop[crop_N[[m]]]=tj
crop
}}
stopCluster(cl)
m=i
ini5=input[,(m)]
crop5=crpp
################################
j=6
source("rainfed2/input_rain.R",echo=TRUE)

w=6
s2a=s2[,1,]
input=s2a
LA=rep(man,npixel)
input=cbind(input,LA)
colnames(input)=nam
input=as.data.frame(input)
if(RC==5){
ini_tropical_cereal=input[,1]=rowMeans(ini[[1]][,1:4,man,fert])
ini_pulses=input[,2]=rowMeans(ini[[2]][,1:4,man,fert])
ini_temperate_root=input[,3]=rowMeans(ini[[3]][,1:4,man,fert])
ini_soybeans=input[,4]=rowMeans(ini[[4]][,1:4,man,fert])
} else {
ini_tropical_cereal=input[,1]=ini[[1]][,RC,man,fert]
ini_pulses=input[,2]=ini[[2]][,RC,man,fert]
ini_temperate_root=input[,3]=ini[[3]][,RC,man,fert]
ini_soybeans=input[,4]=ini[[4]][,RC,man,fert]
}
ini_crop=list(ini_tropical_cereal,ini_pulses,ini_temperate_root,ini_soybeans)


if(fertilization==1){
 co2=input[,38]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,38]+rep(370.47,npixel)/2))}
}

if(fertilization==1){
 cco2=input[,37]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,37]+rep(0,npixel)/2))}
}

input[,37]=cco2
input[,38]=co2
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
source("rainfed2/wls_rain.R",echo=F)
library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction[crop_N[[m]],m]+E_crop[[m]]
crop=matrix(0,nrow=npixel,ncol=1)
crop[crop_N[[m]]]=tj
crop
}}
stopCluster(cl)
m=i
ini6=input[,(m)]
crop6=crpp
################################
j=7
source("rainfed2/input_rain.R",echo=TRUE)

w=7
s2a=s2[,1,]
input=s2a
LA=rep(man,npixel)
input=cbind(input,LA)
colnames(input)=nam
input=as.data.frame(input)
if(RC==5){
ini_tropical_cereal=input[,1]=rowMeans(ini[[1]][,1:4,man,fert])
ini_pulses=input[,2]=rowMeans(ini[[2]][,1:4,man,fert])
ini_temperate_root=input[,3]=rowMeans(ini[[3]][,1:4,man,fert])
ini_soybeans=input[,4]=rowMeans(ini[[4]][,1:4,man,fert])
} else {
ini_tropical_cereal=input[,1]=ini[[1]][,RC,man,fert]
ini_pulses=input[,2]=ini[[2]][,RC,man,fert]
ini_temperate_root=input[,3]=ini[[3]][,RC,man,fert]
ini_soybeans=input[,4]=ini[[4]][,RC,man,fert]
}
ini_crop=list(ini_tropical_cereal,ini_pulses,ini_temperate_root,ini_soybeans)


if(fertilization==1){
 co2=input[,38]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,38]+rep(370.47,npixel)/2))}
}

if(fertilization==1){
 cco2=input[,37]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,37]+rep(0,npixel)/2))}
}

input[,37]=cco2
input[,38]=co2
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
source("rainfed2/wls_rain.R",echo=F)
library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction[crop_N[[m]],m]+E_crop[[m]]
crop=matrix(0,nrow=npixel,ncol=1)
crop[crop_N[[m]]]=tj
crop
}}
stopCluster(cl)
m=i
ini7=input[,(m)]
crop7=crpp
################################
j=8
source("rainfed2/input_rain.R",echo=TRUE)

w=8
s2a=s2[,1,]
input=s2a
LA=rep(man,npixel)
input=cbind(input,LA)
colnames(input)=nam
input=as.data.frame(input)
if(RC==5){
ini_tropical_cereal=input[,1]=rowMeans(ini[[1]][,1:4,man,fert])
ini_pulses=input[,2]=rowMeans(ini[[2]][,1:4,man,fert])
ini_temperate_root=input[,3]=rowMeans(ini[[3]][,1:4,man,fert])
ini_soybeans=input[,4]=rowMeans(ini[[4]][,1:4,man,fert])
} else {
ini_tropical_cereal=input[,1]=ini[[1]][,RC,man,fert]
ini_pulses=input[,2]=ini[[2]][,RC,man,fert]
ini_temperate_root=input[,3]=ini[[3]][,RC,man,fert]
ini_soybeans=input[,4]=ini[[4]][,RC,man,fert]
}
ini_crop=list(ini_tropical_cereal,ini_pulses,ini_temperate_root,ini_soybeans)


if(fertilization==1){
 co2=input[,38]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,38]+rep(370.47,npixel)/2))}
}

if(fertilization==1){
 cco2=input[,37]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,37]+rep(0,npixel)/2))}
}

input[,37]=cco2
input[,38]=co2
ini_crop=list(ini_tropical_cereal,ini_pulses,ini_temperate_root,ini_soybeans)

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
source("rainfed2/wls_rain.R",echo=F)
library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction[crop_N[[m]],m]+E_crop[[m]]
crop=matrix(0,nrow=npixel,ncol=1)
crop[crop_N[[m]]]=tj
crop
}}
stopCluster(cl)
m=i
ini8=input[,(m)]
crop8=crpp
#################
ini=rbind(ini1,ini2,ini3,ini4,ini5,ini6,ini7,ini8)
colnames(ini)=c("ini_tropical_cereal","ini_pulses","ini_temperate_root","ini_soybeans")
colnames(ini1)=colnames(ini2)=colnames(ini3)=colnames(ini4)=colnames(ini5)=colnames(ini6)=colnames(ini7)=colnames(ini8)=colnames(ini)

if(coupling=="GEMINI"){
crop=cbind(crop1,ini1,crop2,ini2,crop3,ini3,crop4,ini4,crop5,ini5,crop6,ini6,crop7,ini7,crop8,ini8)
crop=as.matrix(crop)
crop=array(crop,dim=c(npixel,2*ncol(crop1),8))
j=0
source("rainfed2/aggregate.R",echo=F)
} else {
#crop=cbind(crop1,crop2,crop3,crop4,crop5,crop6,crop7,crop8)
#crop=array(crop,dim=c(npixel,ncol(crop1),8))
crop=cbind(crop1,ini1,crop2,ini2,crop3,ini3,crop4,ini4,crop5,ini5,crop6,ini6,crop7,ini7,crop8,ini8)
crop=as.matrix(crop)
crop=array(crop,dim=c(npixel,2*ncol(crop1),8))

j=0
source("rainfed2/net.R",echo=F)}
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
source("rainfed2/input_rain.R",echo=TRUE)

npixel=59199
#load("rainfed2/k3")
library(ncdf)
load("rainfed2/ini")
load("rainfed2/crop_N")
load("rainfed2/grid_out")
load("rainfed2/crop_form4b")
load("rainfed2/coef4b")
m=i
###################################################
nam=c("ini_tropical_cereal","ini_pulses","ini_temperate_root","ini_soybeans","scld" ,   "wcld"  ,  "spcld" ,  "acld" ,"spre",   
 "wpre" ,   "sppre" ,  "apre" ,   "stmp" ,   "wtmp" ,   "sptmp"  , "atmp",   
"swet" ,   "wwet" ,   "spwet" ,  "awet" ,   "iscld", "iwcld",   "ispcld", 
"iacld","ispre" ,"iwpre" ,  "isppre" , "iapre" ,  "istmp" ,  "iwtmp" , 
"isptmp", "iatmp" ,  "iswet"  , "iwwet" ,  "ispwet" , "iawet" ,  "cco2", "co2" , "soil",  "lat.lpj","LA") 
w=1
s2a=s2[,1,]
input=s2a
LA=rep(man,npixel)
input=cbind(input,LA)
colnames(input)=nam
input=as.data.frame(input)
if(RC==5){
ini_tropical_cereal=input[,1]=rowMeans(ini[[1]][,1:4,man,fert])
ini_pulses=input[,2]=rowMeans(ini[[2]][,1:4,man,fert])
ini_temperate_root=input[,3]=rowMeans(ini[[3]][,1:4,man,fert])
ini_soybeans=input[,4]=rowMeans(ini[[4]][,1:4,man,fert])
} else {
ini_tropical_cereal=input[,1]=ini[[1]][,RC,man,fert]
ini_pulses=input[,2]=ini[[2]][,RC,man,fert]
ini_temperate_root=input[,3]=ini[[3]][,RC,man,fert]
ini_soybeans=input[,4]=ini[[4]][,RC,man,fert]
}
ini_crop=list(ini_tropical_cereal,ini_pulses,ini_temperate_root,ini_soybeans)


if(fertilization==1){
 co2=input[,38]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,38]+rep(370.47,npixel)/2))
}}

if(fertilization==1){
 cco2=input[,37]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,37]+rep(0,npixel)/2))
}}

input[,37]=cco2
input[,38]=co2

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
source("rainfed2/wls_rain.R",echo=F)

######
library(doParallel)
cl <- makeCluster(length(m))
registerDoParallel(cl)
crpp<-{foreach(m=i,.combine=cbind,.verbose=TRUE) %dopar% {
tj=prediction[crop_N[[m]],m]+E_crop[[m]]
crop=matrix(0,nrow=npixel,ncol=1)
crop[crop_N[[m]]]=tj
crop
}}
stopCluster(cl)
m=i
ini=input[,(m)]

if(coupling=="GEMINI"){
crop=cbind(crpp,ini)
crop=as.matrix(crop)
crop=array(crop,dim=c(59199,ncol(crop),1))
source("rainfed2/aggregate.R",echo=F)
} else {
#crop=crpp
#crop=array(crop,dim=c(59199,ncol(crop),1))
crop=cbind(crpp,ini)
crop=as.matrix(crop)
crop=array(crop,dim=c(59199,ncol(crop),1))
source("rainfed2/net.R",echo=F)
}}
