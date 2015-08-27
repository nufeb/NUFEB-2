path="E:\\clim\\management_irrigated_rainfed_co2_on_off\\had2"

dirname="E:\\validation/rain"
out="current"
type="rainfed"
fertilization=1

i=1:4
#########################
new_dir=getwd()
source("conf1.R",echo=F)

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE)
res1=my_crop
#predd1=my_pred

#################
new_dir=getwd()
source("conf2.R",echo=F)

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res2=my_crop
#predd2=my_pred

#############
i=1:4
#########################
new_dir=getwd()
source("conf3.R",echo=F)

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res3=my_crop
#predd3=my_pred

################
i=1:4
#########################
new_dir=getwd()
source("conf4.R",echo=F)

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res4=my_crop
#predd4=my_pred

#############
i=1:4
#########################
new_dir=getwd()
source("conf5.R",echo=F)

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res5=my_crop
#predd5=my_pred

#######################
i=1:4
#########################
new_dir=getwd()
source("conf6.R",echo=F)

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res6=my_crop
#predd6=my_pred

##############
i=1:4
#########################
new_dir=getwd()
source("conf7.R",echo=F)

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res7=my_crop
#predd7=my_pred
my_result3=rbind(res1,res2,res3,res4,res5,res6,res7)
#=rbind(#predd1,#predd2,#predd3,#predd4,#predd5,#predd6,#predd7)
save(my_result3,file=paste(dirname,"my_result3",sep="/"))
#save(myp3,file=paste(dirname,"myp3",sep="/"))
########################################################################RCP3
i=1:4
#########################
#setwd(new_dir)
#rm(list=ls(all=TRUE))
#################
i=1:4
#########################
#new_dir=getwd()
source("conf1.R",echo=F)
rcp="RCP3PD"

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res1=my_crop
#predd1=my_pred

#################
i=1:4
#########################
new_dir=getwd()
source("conf2.R",echo=F)
rcp="RCP3PD"

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res2=my_crop
#predd2=my_pred

#############
i=1:4
#########################
new_dir=getwd()
source("conf3.R",echo=F)
rcp="RCP3PD"

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res3=my_crop
#predd3=my_pred

################
i=1:4
#########################
new_dir=getwd()
source("conf4.R",echo=F)
rcp="RCP3PD"

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res4=my_crop
#predd4=my_pred

#############
i=1:4
#########################
new_dir=getwd()
source("conf5.R",echo=F)
rcp="RCP3PD"

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res5=my_crop

#######################
i=1:4
#########################
new_dir=getwd()
source("conf6.R",echo=F)
rcp="RCP3PD"

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res6=my_crop
#predd6=my_pred

##############
i=1:4
#########################
new_dir=getwd()
source("conf7.R",echo=F)
rcp="RCP3PD"
if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res7=my_crop
##predd7=my_pred

my_result1=rbind(res1,res2,res3,res4,res5,res6,res7)
save(my_result1,file=paste(dirname,"my_result1",sep="/"))
#=rbind(#predd1,#predd2,#predd3,#predd4,#predd5,#predd6,#predd7)
#save(myp1,file=paste(dirname,"myp1",sep="/"))


################################################################################RCP45
i=1:4
#########################
#setwd(new_dir)
#rm(list=ls(all=TRUE))
#################
i=1:4
#########################
#new_dir=getwd()
source("conf1.R",echo=F)
rcp="RCP45"

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res1=my_crop
#predd1=my_pred
#################
i=1:4
#########################
new_dir=getwd()
source("conf2.R",echo=F)
rcp="RCP45"

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res2=my_crop
#predd2=my_pred
#############
i=1:4
#########################
new_dir=getwd()
source("conf3.R",echo=F)
rcp="RCP45"

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res3=my_crop
#predd3=my_pred
################
i=1:4
#########################
new_dir=getwd()
source("conf4.R",echo=F)
rcp="RCP45"

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res4=my_crop
#predd4=my_pred
#############
i=1:4
#########################
new_dir=getwd()
source("conf5.R",echo=F)
rcp="RCP45"

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res5=my_crop
#predd5=my_pred
#######################
i=1:4
#########################
new_dir=getwd()
source("conf6.R",echo=F)
rcp="RCP45"

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res6=my_crop
#predd6=my_pred
##############
i=1:4
#########################
new_dir=getwd()
source("conf7.R",echo=F)
rcp="RCP45"
if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res7=my_crop
#predd7=my_pred
my_result2=rbind(res1,res2,res3,res4,res5,res6,res7)
save(my_result2,file=paste(dirname,"my_result2",sep="/"))

#myp2=rbind(#predd1,#predd2,#predd3,#predd4,#predd5,#predd6,#predd7)
#save(myp2,file=paste(dirname,"my2",sep="/"))


################################################################################RCP85
i=1:4
#########################
#setwd(new_dir)
#rm(list=ls(all=TRUE))
#################
i=1:4
#########################
#new_dir=getwd()
source("conf1.R",echo=F)
rcp="RCP85"

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res1=my_crop
#predd1=my_pred
#################
i=1:4
#########################
new_dir=getwd()
source("conf2.R",echo=F)
rcp="RCP85"

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res2=my_crop
#predd2=my_pred
#############
i=1:4
#########################
new_dir=getwd()
source("conf3.R",echo=F)
rcp="RCP85"

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res3=my_crop
#predd3=my_pred
################
i=1:4
#########################
new_dir=getwd()
source("conf4.R",echo=F)
rcp="RCP85"

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res4=my_crop
#predd4=my_pred
#############
i=1:4
#########################
new_dir=getwd()
source("conf5.R",echo=F)
rcp="RCP85"

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res5=my_crop
#predd5=my_pred
#######################
i=1:4
#########################
new_dir=getwd()
source("conf6.R",echo=F)
rcp="RCP85"

if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res6=my_crop
#predd6=my_pred
##############
i=1:4
#########################
new_dir=getwd()
source("conf7.R",echo=F)
rcp="RCP85"
if(rcp=="RCP3PD") {
 RC=1
} else {
if(rcp=="RCP45") {
RC=2
} else {
if(rcp=="RCP6") {
RC=3
} else {
if(rcp=="RCP85") {
RC=4
}else{
RC=5
}}}}
########################
#########################
if(j==0) {
   w=1:8
} else {
   w=j
}
########
if(type=="rainfed"){
source("rainfed2/CV.R",echo=F)
} else {
source("irrigated/CV.R",echo=F)
}
setwd(new_dir)
#rm(list=ls(all=TRUE))
res7=my_crop
#predd7=my_pred
my_result4=rbind(res1,res2,res3,res4,res5,res6,res7)

################################################################
m1=array(my_result1,c(186,8,7,4,1))###grid X decade X man X crop X rcp1###
m2=array(my_result2,c(186,8,7,4,1))
m3=array(my_result3,c(186,8,7,4,1))
m4=array(my_result4,c(186,8,7,4,1))
library(abind)
#myp4=rbind(#predd1,#predd2,#predd3,#predd4,#predd5,#predd6,#predd7)

m5=abind(m1,m2,m3,m4,along=5)
h1=m5[,,,1,]
h1=aperm(h1,c(1,2,4,3))
h2=m5[,,,2,]
h2=aperm(h2,c(1,2,4,3))
h3=m5[,,,3,]
h3=aperm(h3,c(1,2,4,3))
h4=m5[,,,4,]
h4=aperm(h4,c(1,2,4,3))
had_on_rain=list(h1,h2,h3,h4)
#had_on_rain=aperm(m5,c(1,2,5,4,3))##grid by decade by rcp by 
save(had_on_rain,file=paste(dirname,"had_on_rain",sep="/"))
save(my_result4,file=paste(dirname,"my_result4",sep="/"))
#save(myp4,file=paste(dirname,"myp4",sep="/"))

######################################
