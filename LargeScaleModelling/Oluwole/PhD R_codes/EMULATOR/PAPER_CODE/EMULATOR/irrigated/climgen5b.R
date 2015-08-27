#########R CODE FOR THE FIRST STAGE OF EMULATOR
################################
if(j==0) {

j=1
source("irrigated/input_irri.R",echo=TRUE)
if(fertilization==0){
fert=3:4
} else {
if(fertilization==1) {
fert=1:2
} else {
fert=1:4
}}
load("OIL/ini")
load("OIL/crop_form4")
load("OIL/coef4")
inip=ini
crop_form4p=crop_form4
coef4p=coef4

npixel=59199
load("irrigated/k3")
library(ncdf)
load("irrigated/ini")
load("irrigated/crop_N")
load("irrigated/grid_out")
load("irrigated/crop_form4")
load("irrigated/coef4")
m=i
###################################################
nam=c("scld" ,   "wcld"  ,  "spcld" ,  "acld" ,"spre",   
 "wpre" ,   "sppre" ,  "apre" ,   "stmp" ,   "wtmp" ,   "sptmp"  , "atmp",   
"swet" ,   "wwet" ,   "spwet" ,  "awet" ,   "iscld", "iwcld",   "ispcld", 
"iacld","ispre" ,"iwpre" ,  "isppre" , "iapre" ,  "istmp" ,  "iwtmp" , 
"isptmp", "iatmp" ,  "iswet"  , "iwwet" ,  "ispwet" , "iawet" ,  "cco2",   
"co2" , "soil",  "lat.lpj","LA") 
################################
w=1
LA=rep(man,npixel)
input=cbind(s2,LA)
colnames(input)=nam
input=as.data.frame(input)

if(RC==5){
ini_cereal=rowMeans(ini[,1,1:4,man,fert])
ini_rice=rowMeans(ini[,2,1:4,man,fert])
ini_maize=rowMeans(ini[,3,1:4,man,fert])
ini_oil=rowMeans(inip[,1,1:4,man,fert])
ini_grd=rowMeans(inip[,2,1:4,man,fert])
} else {
ini_cereal=rowMeans(ini[,1,RC,man,fert])
ini_rice=rowMeans(ini[,2,RC,man,fert])
ini_maize=rowMeans(ini[,3,RC,man,fert])
ini_oil=rowMeans(inip[,1,RC,man,fert])
ini_grd=rowMeans(inip[,2,RC,man,fert])
}
ini_crop=cbind(ini_cereal,ini_rice,ini_maize,ini_oil,ini_grd)


if(fertilization==1){
 co2=input[,34]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,34]+rep(370.47,npixel)/2))
}}

if(fertilization==1){
 cco2=input[,33]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,33]+rep(0,npixel)/2))
}}

input[,33]=cco2
input[,34]=co2
#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
predictionp<-{foreach(m=1:2,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4p[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4p[[m]]))
crop22[-crop_N[[3+m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction=cbind(prediction,predictionp)
################
#####include PCA analysis
source("irrigated/wls_irri.R",echo=F)

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
#ini1=input[,(4+m)]
crop1=crpp
################################
j=2
source("irrigated/input_irri.R",echo=TRUE)
w=2
LA=rep(man,npixel)
input=cbind(s2,LA)
colnames(input)=nam
input=as.data.frame(input)

if(fertilization==1){
 co2=input[,34]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,34]+rep(370.47,npixel)/2))
}}

if(fertilization==1){
 cco2=input[,33]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,33]+rep(0,npixel)/2))
}}

input[,33]=cco2
input[,34]=co2
#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
predictionp<-{foreach(m=1:2,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4p[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4p[[m]]))
crop22[-crop_N[[3+m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction=cbind(prediction,predictionp)
################
#####include PCA analysis
source("irrigated/wls_irri.R",echo=F)
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
#ini2=input[,(4+m)]
crop2=crpp
################################
j=3
source("irrigated/input_irri.R",echo=TRUE)

w=3
LA=rep(man,npixel)
input=cbind(s2,LA)
colnames(input)=nam
input=as.data.frame(input)

if(fertilization==1){
 co2=input[,34]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,34]+rep(370.47,npixel)/2))
}}

if(fertilization==1){
 cco2=input[,33]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,33]+rep(0,npixel)/2))
}}

input[,33]=cco2
input[,34]=co2
#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
predictionp<-{foreach(m=1:2,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4p[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4p[[m]]))
crop22[-crop_N[[3+m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction=cbind(prediction,predictionp)
################
#####include PCA analysis
source("irrigated/wls_irri.R",echo=F)
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
#ini3=input[,(4+m)]
crop3=crpp
################################
j=4
source("irrigated/input_irri.R",echo=TRUE)

w=4
LA=rep(man,npixel)
input=cbind(s2,LA)
colnames(input)=nam
input=as.data.frame(input)

if(fertilization==1){
 co2=input[,34]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,34]+rep(370.47,npixel)/2))
}}

if(fertilization==1){
 cco2=input[,33]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,33]+rep(0,npixel)/2))
}}

input[,33]=cco2
input[,34]=co2
#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
predictionp<-{foreach(m=1:2,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4p[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4p[[m]]))
crop22[-crop_N[[3+m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction=cbind(prediction,predictionp)
################
#####include PCA analysis
source("irrigated/wls_irri.R",echo=F)
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
#ini4=input[,(4+m)]
crop4=crpp
################################
j=5
source("irrigated/input_irri.R",echo=TRUE)

w=5
LA=rep(man,npixel)
input=cbind(s2,LA)
colnames(input)=nam
input=as.data.frame(input)

if(fertilization==1){
 co2=input[,34]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,34]+rep(370.47,npixel)/2))
}}

if(fertilization==1){
 cco2=input[,33]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,33]+rep(0,npixel)/2))
}}

input[,33]=cco2
input[,34]=co2
#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
predictionp<-{foreach(m=1:2,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4p[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4p[[m]]))
crop22[-crop_N[[3+m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction=cbind(prediction,predictionp)
################
#####include PCA analysis
source("irrigated/wls_irri.R",echo=F)
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
#ini5=input[,(4+m)]
crop5=crpp
################################
j=6
source("irrigated/input_irri.R",echo=TRUE)

w=6
LA=rep(man,npixel)
input=cbind(s2,LA)
colnames(input)=nam
input=as.data.frame(input)

if(fertilization==1){
 co2=input[,34]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,34]+rep(370.47,npixel)/2))
}}

if(fertilization==1){
 cco2=input[,33]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,33]+rep(0,npixel)/2))
}}

input[,33]=cco2
input[,34]=co2
#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
predictionp<-{foreach(m=1:2,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4p[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4p[[m]]))
crop22[-crop_N[[3+m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction=cbind(prediction,predictionp)
################
#####include PCA analysis
source("irrigated/wls_irri.R",echo=F)
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
#ini6=input[,(4+m)]
crop6=crpp
################################
j=7
source("irrigated/input_irri.R",echo=TRUE)

w=7
LA=rep(man,npixel)
input=cbind(s2,LA)
colnames(input)=nam
input=as.data.frame(input)

if(fertilization==1){
 co2=input[,34]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,34]+rep(370.47,npixel)/2))
}}

if(fertilization==1){
 cco2=input[,33]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,33]+rep(0,npixel)/2))
}}

input[,33]=cco2
input[,34]=co2
#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
predictionp<-{foreach(m=1:2,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4p[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4p[[m]]))
crop22[-crop_N[[3+m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction=cbind(prediction,predictionp)
################
#####include PCA analysis
source("irrigated/wls_irri.R",echo=F)
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
#ini7=input[,(4+m)]
crop7=crpp
################################
j=8
source("irrigated/input_irri.R",echo=TRUE)

w=8
LA=rep(man,npixel)
input=cbind(s2,LA)
colnames(input)=nam
input=as.data.frame(input)

if(fertilization==1){
 co2=input[,34]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,34]+rep(370.47,npixel)/2))
}}

if(fertilization==1){
 cco2=input[,33]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,33]+rep(0,npixel)/2))
}}

input[,33]=cco2
input[,34]=co2
##########################
#########cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
predictionp<-{foreach(m=1:2,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4p[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4p[[m]]))
crop22[-crop_N[[3+m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction=cbind(prediction,predictionp)
################
#####include PCA analysis
source("irrigated/wls_irri.R",echo=F)
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
#ini8=input[,(4+m)]
crop8=crpp
#################
ini1=ini2=ini3=ini4=ini5=ini6=ini7=ini8=ini_crop
#ini=rbind(ini1,ini2,ini3,ini4,ini5,ini6,ini7,ini8)
#colnames(ini)=c("ini_cereal","ini_rice","ini_maize","ini_oil","ini_grd")
#colnames(ini1)=colnames(ini2)=colnames(ini3)=colnames(ini4)=colnames(ini5)=colnames(ini6)=colnames(ini7)=colnames(ini8)=colnames(ini)

if(resolution=="LOW"){
crop=cbind(crop1,ini1,crop2,ini2,crop3,ini3,crop4,ini4,crop5,ini5,crop6,ini6,crop7,ini7,crop8,ini8)
crop=as.matrix(crop)
crop=array(crop,dim=c(npixel,2*ncol(crop1),8))
j=0
source("irrigated/aggregate.R",echo=F)
} else {
#crop=cbind(crop1,crop2,crop3,crop4,crop5,crop6,crop7,crop8)
#crop=array(crop,dim=c(npixel,ncol(crop1),8))
crop=cbind(crop1,ini1,crop2,ini2,crop3,ini3,crop4,ini4,crop5,ini5,crop6,ini6,crop7,ini7,crop8,ini8)
crop=as.matrix(crop)
crop=array(crop,dim=c(npixel,2*ncol(crop1),8))

j=0
source("irrigated/net.R",echo=F)}
###############################END
} else {
if(fertilization==0){
fert=3:4
} else {
if(fertilization==1) {
fert=1:2
} else {
fert=1:4
}}
j=j
source("irrigated/input_irri.R",echo=TRUE)
load("OIL/ini")
load("OIL/crop_N")
load("OIL/crop_form4")
load("OIL/coef4")
inip=ini
crop_Np=crop_N
crop_form4p=crop_form4
coef4p=coef4

npixel=59199
#load("irrigated/k3")
library(ncdf)
load("irrigated/ini")
load("irrigated/crop_N")
load("irrigated/grid_out")
load("irrigated/crop_form4")
load("irrigated/coef4")
m=i
###################################################
nam=c("scld" ,   "wcld"  ,  "spcld" ,  "acld" ,"spre",   
 "wpre" ,   "sppre" ,  "apre" ,   "stmp" ,   "wtmp" ,   "sptmp"  , "atmp",   
"swet" ,   "wwet" ,   "spwet" ,  "awet" ,   "iscld", "iwcld",   "ispcld", 
"iacld","ispre" ,"iwpre" ,  "isppre" , "iapre" ,  "istmp" ,  "iwtmp" , 
"isptmp", "iatmp" ,  "iswet"  , "iwwet" ,  "ispwet" , "iawet" ,  "cco2",   
"co2" , "soil",  "lat.lpj","LA") 
w=1

LA=rep(man,npixel)
input=cbind(s2,LA)
colnames(input)=nam
input=as.data.frame(input)

if(RC==5){
ini_cereal=rowMeans(ini[,1,1:4,man,fert])
ini_rice=rowMeans(ini[,2,1:4,man,fert])
ini_maize=rowMeans(ini[,3,1:4,man,fert])
ini_oil=rowMeans(inip[,1,1:4,man,fert])
ini_grd=rowMeans(inip[,2,1:4,man,fert])
} else {
ini_cereal=rowMeans(ini[,1,RC,man,fert])
ini_rice=rowMeans(ini[,2,RC,man,fert])
ini_maize=rowMeans(ini[,3,RC,man,fert])
ini_oil=rowMeans(inip[,1,RC,man,fert])
ini_grd=rowMeans(inip[,2,RC,man,fert])
}
ini_crop=cbind(ini_cereal,ini_rice,ini_maize,ini_oil,ini_grd)

if(fertilization==1){
 co2=input[,34]
} else {
if(fertilization==0) {
co2=rep(370.47,npixel)
} else {
co2=fertilization*((input[,34]+rep(370.47,npixel)/2))
}}

if(fertilization==1){
 cco2=input[,33]
} else {
if(fertilization==0) {
cco2=rep(0,npixel)
} else {
cco2=fertilization*((input[,33]+rep(0,npixel)/2))
}}

input[,33]=cco2
input[,34]=co2
################cereal with other crops
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
prediction<-{foreach(m=1:3,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4[[m]]))
crop22[-crop_N[[m]]]=0
crop1=crop22
}}
stopCluster(cl)
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
predictionp<-{foreach(m=1:2,.combine = cbind,.verbose=TRUE) %dopar% {
crop=matrix(0,ncol=1,nrow=npixel)
ras=terms(formula(paste("~",paste(crop_form4p[[m]]))[2]),input,keep.order = TRUE,specials=NULL)
x=model.matrix(ras,input)
crop22 <- as.vector(x%*%as.numeric(coef4p[[m]]))
crop22[-crop_N[[3+m]]]=0
crop1=crop22
}}
stopCluster(cl)
prediction=cbind(prediction,predictionp)

#################################################
#####include PCA analysis
source("irrigated/wls_irri.R",echo=F)

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
#ini=input[,(4+m)]
ini=ini_crop
if(resolution=="LOW"){
crop=cbind(crpp,ini)
crop=as.matrix(crop)
crop=array(crop,dim=c(59199,ncol(crop),1))
source("irrigated/aggregate.R",echo=F)
} else {
#crop=crpp
#crop=array(crop,dim=c(59199,ncol(crop),1))
crop=cbind(crpp,ini)
crop=as.matrix(crop)
crop=array(crop,dim=c(59199,ncol(crop),1))

source("irrigated/net.R",echo=F)
}}