#########aggregate lpj hadgem, for masked
load("rainfed2/crop_err4b")
load("rainfed2/crop_lpj4b")
load("rainfed2/masked/N1")
library(rworldmap)
library(raster)
load("irrigated2/masked/crop_N")
#load("rainfed2/r3a")
#load("rainfed2/grid_out")
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
E=crop_err4b[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
mydata[-N1[[i]],]=0
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
S1=array(residual,c(186,8,2,7,2))
 
i=2
E=crop_err4b[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
mydata[-N1[[i]],]=0
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
S2=array(residual,c(186,8,2,7,2))

i=3
E=crop_err4b[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
mydata[-N1[[i]],]=0
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
S3=array(residual,c(186,8,2,7,2))

i=4
E=crop_err4b[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
mydata[-N1[[i]],]=0
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
S4=array(residual,c(186,8,2,7,2))
residual_rain_gfdl=list(S1,S2,S3,S4)

dirname="E:\\GPemulator\\masked\\rainfed2\\residual"

save(residual_rain_gfdl,file=paste(dirname,"residual_rain_gfdl",sep="/"))
rm(S1,S2,S3,S4,residual_rain_gfdl)


#####rain lpj
library(rworldmap)
library(raster)
load("irrigated2/masked/crop_N")
#load("rainfed2/r3a")
#load("rainfed2/grid_out")
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
E=crop_lpj4b[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
mydata[-N1[[i]],]=0
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
S1=array(residual,c(186,8,2,7,2))
 
i=2
E=crop_lpj4b[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
mydata[-N1[[i]],]=0
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
S2=array(residual,c(186,8,2,7,2))

i=3
E=crop_lpj4b[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
mydata[-N1[[i]],]=0
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
S3=array(residual,c(186,8,2,7,2))

i=4
E=crop_lpj4b[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
mydata[-N1[[i]],]=0
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
S4=array(residual,c(186,8,2,7,2))
lpj_rain_gfdl=list(S1,S2,S3,S4)

dirname="E:\\GPemulator\\masked\\rainfed2\\residual"

save(lpj_rain_gfdl,file=paste(dirname,"lpj_rain_gfdl",sep="/"))
rm(S1,S2,S3,S4,lpj_rain_gfdl,crop_err4b,crop_lpj4b)






#############################################IRRIGATED#############################################
#########aggregate residual and prediction for masked
load("irrigated2/crop_err4b")
load("irrigated2/crop_lpj4b")
load("irrigated2/masked/N1")
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
E=crop_err4b[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
mydata[-N1[[i]],]=0
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
S1=array(residual,c(186,8,2,7,2))
 
i=2
E=crop_err4b[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
mydata[-N1[[i]],]=0
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
S2=array(residual,c(186,8,2,7,2))
i=3
E=crop_err4b[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
mydata[-N1[[i]],]=0
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
S3=array(residual,c(186,8,2,7,2))
i=4
E=crop_err4b[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
mydata[-N1[[i]],]=0
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
S4=array(residual,c(186,8,2,7,2))
residual_irri_gfdl=list(S1,S2,S3,S4)
 dirname="E:\\GPemulator\\masked\\irrigated2\\residual"
save(residual_irri_gfdl,file=paste(dirname,"residual_irri_gfdl",sep="/"))
rm(S1,S2,S3,S4,residual_irri_gfdl,crop_err4b)


####
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
E=crop_lpj4b[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
mydata[-N1[[i]],]=0
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
S1=array(residual,c(186,8,2,7,2))
 
i=2
E=crop_lpj4b[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
mydata[-N1[[i]],]=0
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
S2=array(residual,c(186,8,2,7,2))
i=3
E=crop_lpj4b[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
mydata[-N1[[i]],]=0
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
S3=array(residual,c(186,8,2,7,2))
i=4
E=crop_lpj4b[[i]]
crop=array(E,c(59199,np,1))
mydata=crop[,,1]
mydata[-N1[[i]],]=0
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
S4=array(residual,c(186,8,2,7,2))
lpj_irri_gfdl=list(S1,S2,S3,S4)
 dirname="E:\\GPemulator\\masked\\irrigated2\\residual"
save(lpj_irri_gfdl,file=paste(dirname,"lpj_irri_gfdl",sep="/"))
rm(S1,S2,S3,S4,lpj_irri_gfdl,crop_lpj4b)

