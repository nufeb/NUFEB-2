library(rworldmap)
library(raster)
load("irrigated2/r3a")
load("irrigated2/grid_out")
grid_len=nrow(grid_out)
load("irrigated2/Country2")
lon2=seq(-179.75,179.75,.5)
lat2=seq(-89.75,89.75,.5)
lat2a=seq(89.75,-89.75,-.5)
s3a=expand.grid(lon2,lat2a)
np=4
if(out=="current"){
mydata=prediction
mydata[-N1[[1]],1]=0
mydata[-N1[[2]],2]=0
mydata[-N1[[3]],3]=0
mydata[-N1[[4]],4]=0
} else {
mydata=prediction
}

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
prediction2 <- as.matrix(data.frame(matrix(unlist(crp1[1:np]), nrow=186)))
prediction2=array(prediction2,c(186,4))