###################crops  pft_harvest.pft.bin

dirname="/padata/alpha/users/ooyebamiji/NEW"
#setwd("/padata/alpha/users/ooyebamiji/management_and_CO2/")


library(doParallel)
cl <- makeCluster(7)
registerDoParallel(cl)
my_list=cbind("LAI_1","LAI_2","LAI_3","LAI_4","LAI_5","LAI_6","LAI_7")
#my_list=cbind("path/LAI_1","path/LAI_2","path/LAI_3","path/LAI_4","path/LAI_5","path/LAI_6","path/LAI_7")
pft3<-{foreach(i=1:ncol(my_list),.combine =cbind,.verbose=T,.errorhandling="stop") %dopar% {
L=my_list[,i]
L=paste("management_and_CO2",L,sep="/")
GCM="cccma-cgcm31"
#crop.yield=sub("L",L,"L/RCP3PD/GCM/pft_harvest.pft.bin")
crop.yield=sub("GCM",GCM,(sub("L",L,"L/RCP3PD/GCM/pft_harvest.pft.bin")))
library(raster)
sizeof.data=4
j=1
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz1=t(apply(cp1[,,1:nyears],1,rowMeans))

j=2
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz2=t(apply(cp1[,,1:nyears],1,rowMeans))

j=3
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz3=t(apply(cp1[,,1:nyears],1,rowMeans))

j=4
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz4=t(apply(cp1[,,1:nyears],1,rowMeans))

j=5
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz5=t(apply(cp1[,,1:nyears],1,rowMeans))

j=6
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz6=t(apply(cp1[,,1:nyears],1,rowMeans))

j=7
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz7=t(apply(cp1[,,1:nyears],1,rowMeans))

j=8
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz8=t(apply(cp1[,,1:nyears],1,rowMeans))

j=9
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz9=t(apply(cp1[,,1:nyears],1,rowMeans))

g1=rbind(zz2,zz3,zz4,zz5,zz6,zz7,zz8,zz9)
g2=rbind(zz1,zz1,zz1,zz1,zz1,zz1,zz1,zz1)
g3=g1-g2
rain=cbind(g3[,1:3],rowSums(g3[,9:12]),g2[,1:3],rowSums(g2[,9:12]))
irri=cbind(g3[,17:19],rowSums(g3[,25:28]),g2[,17:19],rowSums(g2[,25:28]))
dat=cbind(rain,irri)
}}
#save(CFT3,file=paste(dirname,"CFT3",sep="/"))
stopCluster(cl)
#####################

#####################END
###################crops  pft_harvest.pft.bin
# cd "/data2/physics/paftp/LPJML/managment_and_CO2/LAI1"
#file.path("/data2/physics/paf

library(doParallel)
cl <- makeCluster(7)
registerDoParallel(cl)
my_list=cbind("LAI_1","LAI_2","LAI_3","LAI_4","LAI_5","LAI_6","LAI_7")
pft45<-{foreach(i=1:ncol(my_list),.combine =cbind,.verbose=T,.errorhandling="stop") %dopar% {
L=my_list[,i]
L=paste("management_and_CO2",L,sep="/")
GCM="cccma-cgcm31"
#crop.yield=sub("L",L,"L/RCP3PD/GCM/pft_harvest.pft.bin")
crop.yield=sub("GCM",GCM,(sub("L",L,"L/RCP45/GCM/pft_harvest.pft.bin")))
library(raster)
sizeof.data=4
j=1
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz1=t(apply(cp1[,,1:nyears],1,rowMeans))

j=2
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz2=t(apply(cp1[,,1:nyears],1,rowMeans))

j=3
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz3=t(apply(cp1[,,1:nyears],1,rowMeans))

j=4
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz4=t(apply(cp1[,,1:nyears],1,rowMeans))

j=5
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz5=t(apply(cp1[,,1:nyears],1,rowMeans))

j=6
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz6=t(apply(cp1[,,1:nyears],1,rowMeans))

j=7
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz7=t(apply(cp1[,,1:nyears],1,rowMeans))

j=8
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz8=t(apply(cp1[,,1:nyears],1,rowMeans))

j=9
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz9=t(apply(cp1[,,1:nyears],1,rowMeans))

g1=rbind(zz2,zz3,zz4,zz5,zz6,zz7,zz8,zz9)
g2=rbind(zz1,zz1,zz1,zz1,zz1,zz1,zz1,zz1)
g3=g1-g2
rain=cbind(g3[,1:3],rowSums(g3[,9:12]),g2[,1:3],rowSums(g2[,9:12]))
irri=cbind(g3[,17:19],rowSums(g3[,25:28]),g2[,17:19],rowSums(g2[,25:28]))
dat=cbind(rain,irri)
}}
#save(CFT45,file=paste(dirname,"CFT45",sep="/"))
stopCluster(cl)
#####################

#####################END
###################crops  pft_harvest.pft.bin
# cd "/data2/physics/paftp/LPJML/managment_and_CO2/LAI1"
#file.path("/data2/physics/paf

library(doParallel)
cl <- makeCluster(7)
registerDoParallel(cl)
my_list=cbind("LAI_1","LAI_2","LAI_3","LAI_4","LAI_5","LAI_6","LAI_7")
pft6<-{foreach(i=1:ncol(my_list),.combine =cbind,.verbose=T,.errorhandling="stop") %dopar% {
L=my_list[,i]
L=paste("management_and_CO2",L,sep="/")
GCM="cccma-cgcm31"
#crop.yield=sub("L",L,"L/RCP3PD/GCM/pft_harvest.pft.bin")
crop.yield=sub("GCM",GCM,(sub("L",L,"L/RCP6/GCM/pft_harvest.pft.bin")))
library(raster)
sizeof.data=4
j=1
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz1=t(apply(cp1[,,1:nyears],1,rowMeans))

j=2
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz2=t(apply(cp1[,,1:nyears],1,rowMeans))

j=3
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz3=t(apply(cp1[,,1:nyears],1,rowMeans))

j=4
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz4=t(apply(cp1[,,1:nyears],1,rowMeans))

j=5
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz5=t(apply(cp1[,,1:nyears],1,rowMeans))

j=6
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz6=t(apply(cp1[,,1:nyears],1,rowMeans))

j=7
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz7=t(apply(cp1[,,1:nyears],1,rowMeans))

j=8
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz8=t(apply(cp1[,,1:nyears],1,rowMeans))

j=9
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz9=t(apply(cp1[,,1:nyears],1,rowMeans))

g1=rbind(zz2,zz3,zz4,zz5,zz6,zz7,zz8,zz9)
g2=rbind(zz1,zz1,zz1,zz1,zz1,zz1,zz1,zz1)
g3=g1-g2
rain=cbind(g3[,1:3],rowSums(g3[,9:12]),g2[,1:3],rowSums(g2[,9:12]))
irri=cbind(g3[,17:19],rowSums(g3[,25:28]),g2[,17:19],rowSums(g2[,25:28]))
dat=cbind(rain,irri)
}}
#save(CFT6,file=paste(dirname,"CFT6",sep="/"))
stopCluster(cl)
#####################

#####################END
###################crops  pft_harvest.pft.bin
# cd "/data2/physics/paftp/LPJML/managment_and_CO2/LAI1"
#file.path("/data2/physics/paf

library(doParallel)
cl <- makeCluster(7)
registerDoParallel(cl)
my_list=cbind("LAI_1","LAI_2","LAI_3","LAI_4","LAI_5","LAI_6","LAI_7")
pft85<-{foreach(i=1:ncol(my_list),.combine =cbind,.verbose=T,.errorhandling="stop") %dopar% {
L=my_list[,i]
L=paste("management_and_CO2",L,sep="/")
GCM="cccma-cgcm31"
#crop.yield=sub("L",L,"L/RCP3PD/GCM/pft_harvest.pft.bin")
crop.yield=sub("GCM",GCM,(sub("L",L,"L/RCP85/GCM/pft_harvest.pft.bin")))
library(raster)
sizeof.data=4
j=1
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz1=t(apply(cp1[,,1:nyears],1,rowMeans))

j=2
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz2=t(apply(cp1[,,1:nyears],1,rowMeans))

j=3
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz3=t(apply(cp1[,,1:nyears],1,rowMeans))

j=4
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz4=t(apply(cp1[,,1:nyears],1,rowMeans))

j=5
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz5=t(apply(cp1[,,1:nyears],1,rowMeans))

j=6
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz6=t(apply(cp1[,,1:nyears],1,rowMeans))

j=7
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz7=t(apply(cp1[,,1:nyears],1,rowMeans))

j=8
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz8=t(apply(cp1[,,1:nyears],1,rowMeans))

j=9
library(raster)
sizeof.data=4
nyear.first=2005+(10*(j-1))
first.year_simulation<-2001
nyear_first<- nyear.first - first.year_simulation          
nbands <- 32                
npixel <- 59199
nyears=10
file.data <-file(crop.yield,"rb")
skip <- nyear_first * nbands * npixel * sizeof.data
seek(file.data, where = skip, origin = "start", rw = "read")
# cp=matrix(0,nrow=(npixel*nbands*nyears ),ncol=2)
cp=readBin(file.data, numeric(), n=(npixel*nbands*nyears ), size=4)
close(file.data)
cp1=array(cp,dim=c(npixel,nbands,nyears))
zz9=t(apply(cp1[,,1:nyears],1,rowMeans))

g1=rbind(zz2,zz3,zz4,zz5,zz6,zz7,zz8,zz9)
g2=rbind(zz1,zz1,zz1,zz1,zz1,zz1,zz1,zz1)
g3=g1-g2
rain=cbind(g3[,1:3],rowSums(g3[,9:12]),g2[,1:3],rowSums(g2[,9:12]))
irri=cbind(g3[,17:19],rowSums(g3[,25:28]),g2[,17:19],rowSums(g2[,25:28]))
dat=cbind(rain,irri)
}}
#save(CFT85,file=paste(dirname,"CFT85",sep="/"))
stopCluster(cl)
#####################

#####################END########################
#cft=rbind(CFT3,CFT45,CFT6,CFT85)
c1=array(pft3,dim=c(473592,16,1,7))
c2=array(pft45,dim=c(473592,16,1,7))
c3=array(pft6,dim=c(473592,16,1,7))
c4=array(pft85,dim=c(473592,16,1,7))
library(abind)
cft_cccma=abind(c1,c2,c3,c4,along=3)
save(cft_cccma,file=paste(dirname,"cft_cccma",sep="/"))
########################################END

