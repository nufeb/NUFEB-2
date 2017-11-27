#########aggregate lpj hadgem, masked
oil_rain=lpj_oil_rain[,,,,5]
#oil_irri=lpj_oil_irri[,,,,5]
lpj_had_rain[[4]]=oil_rain
#lpj_had_irri[[4]]=oil_irri

library(rworldmap)
library(raster)
load("irrigated2/masked/crop_N")
load("rainfed2/r3a")
load("rainfed2/grid_out")
grid_len=nrow(grid_out)
#load("rainfed2/Country2")
lon2=seq(-179.75,179.75,.5)
lat2=seq(-89.75,89.75,.5)
lat2a=seq(89.75,-89.75,-.5)
s3a=expand.grid(lon2,lat2a)
np=8*4*7
#crop=array(dat,c(59199,np,1))
#load("rainfed2/crop_lpj4b")
#load("C:\\Users\\olu\\Desktop\\ALL_EMULATORS/error_rain")
i=1
E=lpj_had_rain[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
w=1

lpj.raster <- raster(nrows=360, ncols=720, xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat +datum=WGS84")
lpj.ras <- lpj.raster
yaya=xyFromCell(lpj.ras, cell=r3a, spatial=FALSE)
grid.area=area(lpj.ras)
lpj.area=grid.area[r3a]
lpj.area=lpj.area*10^6
mat=matrix(0,nrow=259200,ncol=(np+2))
s3a=as.matrix(s3a)
mat[,1:2]=s3a
mat[r3a,3:(np+2)]=mydata[,1:np]
mat[r3a,3:(np+2)]=mat[r3a,3:(np+2)]*lpj.area
#mat[r3a,3:30]=mat[r3a,3:(np+2)]
mat=as.data.frame(mat)
coordinates(mat) = c("V1", "V2")
gridded(mat) <- TRUE
mat1 = as(mat, "SpatialGridDataFrame")
crp1=matrix(0,nrow=186,ncol=np)
for(t in 1:np){
crp1[t]=aggregateHalfDegreeGridToCountries(mat1[t], aggregateOption = "sum")[2]
}
residual <- as.matrix(data.frame(matrix(unlist(crp1[1:np]), nrow=186)))
S1=array(residual,c(186,8,4,7))
 
i=2
E=lpj_had_rain[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
w=1

lpj.raster <- raster(nrows=360, ncols=720, xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat +datum=WGS84")
lpj.ras <- lpj.raster
yaya=xyFromCell(lpj.ras, cell=r3a, spatial=FALSE)
grid.area=area(lpj.ras)
lpj.area=grid.area[r3a]
lpj.area=lpj.area*10^6
mat=matrix(0,nrow=259200,ncol=(np+2))
s3a=as.matrix(s3a)
mat[,1:2]=s3a
mat[r3a,3:(np+2)]=mydata[,1:np]
mat[r3a,3:(np+2)]=mat[r3a,3:(np+2)]*lpj.area
#mat[r3a,3:30]=mat[r3a,3:(np+2)]
mat=as.data.frame(mat)
coordinates(mat) = c("V1", "V2")
gridded(mat) <- TRUE
mat1 = as(mat, "SpatialGridDataFrame")
crp1=matrix(0,nrow=186,ncol=np)
for(t in 1:np){
crp1[t]=aggregateHalfDegreeGridToCountries(mat1[t], aggregateOption = "sum")[2]
}
residual <- as.matrix(data.frame(matrix(unlist(crp1[1:np]), nrow=186)))
S2=array(residual,c(186,8,4,7))

i=3
E=lpj_had_rain[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
w=1

lpj.raster <- raster(nrows=360, ncols=720, xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat +datum=WGS84")
lpj.ras <- lpj.raster
yaya=xyFromCell(lpj.ras, cell=r3a, spatial=FALSE)
grid.area=area(lpj.ras)
lpj.area=grid.area[r3a]
lpj.area=lpj.area*10^6
mat=matrix(0,nrow=259200,ncol=(np+2))
s3a=as.matrix(s3a)
mat[,1:2]=s3a
mat[r3a,3:(np+2)]=mydata[,1:np]
mat[r3a,3:(np+2)]=mat[r3a,3:(np+2)]*lpj.area
#mat[r3a,3:30]=mat[r3a,3:(np+2)]
mat=as.data.frame(mat)
coordinates(mat) = c("V1", "V2")
gridded(mat) <- TRUE
mat1 = as(mat, "SpatialGridDataFrame")
crp1=matrix(0,nrow=186,ncol=np)
for(t in 1:np){
crp1[t]=aggregateHalfDegreeGridToCountries(mat1[t], aggregateOption = "sum")[2]
}
residual <- as.matrix(data.frame(matrix(unlist(crp1[1:np]), nrow=186)))
S3=array(residual,c(186,8,4,7))

i=4
E=lpj_had_rain[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
w=1

lpj.raster <- raster(nrows=360, ncols=720, xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat +datum=WGS84")
lpj.ras <- lpj.raster
yaya=xyFromCell(lpj.ras, cell=r3a, spatial=FALSE)
grid.area=area(lpj.ras)
lpj.area=grid.area[r3a]
lpj.area=lpj.area*10^6
mat=matrix(0,nrow=259200,ncol=(np+2))
s3a=as.matrix(s3a)
mat[,1:2]=s3a
mat[r3a,3:(np+2)]=mydata[,1:np]
mat[r3a,3:(np+2)]=mat[r3a,3:(np+2)]*lpj.area
#mat[r3a,3:30]=mat[r3a,3:(np+2)]
mat=as.data.frame(mat)
coordinates(mat) = c("V1", "V2")
gridded(mat) <- TRUE
mat1 = as(mat, "SpatialGridDataFrame")
crp1=matrix(0,nrow=186,ncol=np)
for(t in 1:np){
crp1[t]=aggregateHalfDegreeGridToCountries(mat1[t], aggregateOption = "sum")[2]
}
residual <- as.matrix(data.frame(matrix(unlist(crp1[1:np]), nrow=186)))
S4=array(residual,c(186,8,4,7))
S_rain_lpj=list(S1,S2,S3,S4)

dirname="E:\\GPemulator\\masked\\irrigated2\\residual"

save(S_rain_lpj,file=paste(dirname,"S_rain_lpj",sep="/"))
rm(S1,S2,S3,S4,S_rain_lpj)
#############################################RAIN
#########aggregate residual and prediction for masked
library(rworldmap)
library(raster)
#load("rainfed2/masked/crop_N")
#load("rainfed2/r3a")
#load("rainfed2/grid_out")
grid_len=nrow(grid_out)
#load("rainfed2/Country2")
lon2=seq(-179.75,179.75,.5)
lat2=seq(-89.75,89.75,.5)
lat2a=seq(89.75,-89.75,-.5)
s3a=expand.grid(lon2,lat2a)

#crop=array(dat,c(59199,np,1))
#load("rainfed2/crop_lpj4b")
#load("C:\\Users\\olu\\Desktop\\ALL_EMULATORS/error_rain")
i=1
E=lpj_had_irri[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
w=1

lpj.raster <- raster(nrows=360, ncols=720, xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat +datum=WGS84")
lpj.ras <- lpj.raster
yaya=xyFromCell(lpj.ras, cell=r3a, spatial=FALSE)
grid.area=area(lpj.ras)
lpj.area=grid.area[r3a]
lpj.area=lpj.area*10^6
mat=matrix(0,nrow=259200,ncol=(np+2))
s3a=as.matrix(s3a)
mat[,1:2]=s3a
mat[r3a,3:(np+2)]=mydata[,1:np]
mat[r3a,3:(np+2)]=mat[r3a,3:(np+2)]*lpj.area
#mat[r3a,3:30]=mat[r3a,3:(np+2)]
mat=as.data.frame(mat)
coordinates(mat) = c("V1", "V2")
gridded(mat) <- TRUE
mat1 = as(mat, "SpatialGridDataFrame")
crp1=matrix(0,nrow=186,ncol=np)
for(t in 1:np){
crp1[t]=aggregateHalfDegreeGridToCountries(mat1[t], aggregateOption = "sum")[2]
}
residual <- as.matrix(data.frame(matrix(unlist(crp1[1:np]), nrow=186)))
S1=array(residual,c(186,8,4,7))
 
i=2
E=lpj_had_irri[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
w=1

lpj.raster <- raster(nrows=360, ncols=720, xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat +datum=WGS84")
lpj.ras <- lpj.raster
yaya=xyFromCell(lpj.ras, cell=r3a, spatial=FALSE)
grid.area=area(lpj.ras)
lpj.area=grid.area[r3a]
lpj.area=lpj.area*10^6
mat=matrix(0,nrow=259200,ncol=(np+2))
s3a=as.matrix(s3a)
mat[,1:2]=s3a
mat[r3a,3:(np+2)]=mydata[,1:np]
mat[r3a,3:(np+2)]=mat[r3a,3:(np+2)]*lpj.area
#mat[r3a,3:30]=mat[r3a,3:(np+2)]
mat=as.data.frame(mat)
coordinates(mat) = c("V1", "V2")
gridded(mat) <- TRUE
mat1 = as(mat, "SpatialGridDataFrame")
crp1=matrix(0,nrow=186,ncol=np)
for(t in 1:np){
crp1[t]=aggregateHalfDegreeGridToCountries(mat1[t], aggregateOption = "sum")[2]
}
residual <- as.matrix(data.frame(matrix(unlist(crp1[1:np]), nrow=186)))
S2=array(residual,c(186,8,4,7))
i=3
E=lpj_had_irri[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
w=1

lpj.raster <- raster(nrows=360, ncols=720, xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat +datum=WGS84")
lpj.ras <- lpj.raster
yaya=xyFromCell(lpj.ras, cell=r3a, spatial=FALSE)
grid.area=area(lpj.ras)
lpj.area=grid.area[r3a]
lpj.area=lpj.area*10^6
mat=matrix(0,nrow=259200,ncol=(np+2))
s3a=as.matrix(s3a)
mat[,1:2]=s3a
mat[r3a,3:(np+2)]=mydata[,1:np]
mat[r3a,3:(np+2)]=mat[r3a,3:(np+2)]*lpj.area
#mat[r3a,3:30]=mat[r3a,3:(np+2)]
mat=as.data.frame(mat)
coordinates(mat) = c("V1", "V2")
gridded(mat) <- TRUE
mat1 = as(mat, "SpatialGridDataFrame")
crp1=matrix(0,nrow=186,ncol=np)
for(t in 1:np){
crp1[t]=aggregateHalfDegreeGridToCountries(mat1[t], aggregateOption = "sum")[2]
}
residual <- as.matrix(data.frame(matrix(unlist(crp1[1:np]), nrow=186)))
S3=array(residual,c(186,8,4,7))
i=4
E=lpj_had_irri[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
w=1

lpj.raster <- raster(nrows=360, ncols=720, xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat +datum=WGS84")
lpj.ras <- lpj.raster
yaya=xyFromCell(lpj.ras, cell=r3a, spatial=FALSE)
grid.area=area(lpj.ras)
lpj.area=grid.area[r3a]
lpj.area=lpj.area*10^6
mat=matrix(0,nrow=259200,ncol=(np+2))
s3a=as.matrix(s3a)
mat[,1:2]=s3a
mat[r3a,3:(np+2)]=mydata[,1:np]
mat[r3a,3:(np+2)]=mat[r3a,3:(np+2)]*lpj.area
#mat[r3a,3:30]=mat[r3a,3:(np+2)]
mat=as.data.frame(mat)
coordinates(mat) = c("V1", "V2")
gridded(mat) <- TRUE
mat1 = as(mat, "SpatialGridDataFrame")
crp1=matrix(0,nrow=186,ncol=np)
for(t in 1:np){
crp1[t]=aggregateHalfDegreeGridToCountries(mat1[t], aggregateOption = "sum")[2]
}
residual <- as.matrix(data.frame(matrix(unlist(crp1[1:np]), nrow=186)))
S4=array(residual,c(186,8,4,7))
S_irri_lpj=list(S1,S2,S3,S4)
 dirname="E:\\GPemulator\\masked\\irrigated2\\residual"
save(S_irri_lpj,file=paste(dirname,"S_irri_lpj",sep="/"))
