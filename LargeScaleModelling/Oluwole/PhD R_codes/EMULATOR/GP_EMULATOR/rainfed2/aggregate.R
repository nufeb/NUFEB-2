######### LPJmL crop yield aggregation to country levels
library(rworldmap)
library(raster)
load("rainfed2/r3a")
load("rainfed2/grid_out")
grid_len=nrow(grid_out)
load("rainfed/Country2")
lon2=seq(-179.75,179.75,.5)
lat2=seq(-89.75,89.75,.5)
lat2a=seq(89.75,-89.75,-.5)
s3a=expand.grid(lon2,lat2a)
dirname <- getwd()
############
if(j==0) {
    w=1:8
} else {
   w=j;mydata=crop[,,1]
}
library(doParallel)
cl <- makeCluster(length(w))
registerDoParallel(cl)
result_model <-{foreach(r=1:length(w),.combine =cbind,.verbose=FALSE,.packages=c("raster","maptools","fields","sp","rworldmap")) %dopar% {
mydata=crop[,,r]
#library(rworldmap)
#library(raster)
lpj.raster <- raster(nrows=360, ncols=720, xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat +datum=WGS84")
lpj.ras <- lpj.raster
yaya=xyFromCell(lpj.ras, cell=r3a, spatial=FALSE)
grid.area=area(lpj.ras)
lpj.area=grid.area[r3a]
lpj.area=lpj.area*10^6
mat=matrix(0,nrow=259200,ncol=10)
s3a=as.matrix(s3a)
mat[,1:2]=s3a
mat[r3a,3:6]=mydata[,1:4]
mat[r3a,7:10]=mydata[,5:8]
colnames(mat)=c("Lon","Lat","cereal","rice","maize","oil","ini_cereal","ini_rice","ini_maize","ini_oil")
 mat[r3a,3:10]=mat[r3a,3:10]*lpj.area
 mat=as.data.frame(mat)
coordinates(mat) = c("Lon", "Lat")
gridded(mat) <- TRUE
mat1 = as(mat, "SpatialGridDataFrame")
crp1=matrix(0,nrow=186,ncol=8)
for(t in 1:8){
crp1[t]=aggregateHalfDegreeGridToCountries(mat1[t], aggregateOption = "sum")[2]
}
crp <- data.frame(matrix(unlist(crp1[1:8]), nrow=186))
colnames(crp)=c("cereal","rice","maize","oil_max","ini_cereal","ini_rice","ini_maize","ini_oil_max")
res=cbind(Country2,crp)
}}
stopCluster(cl)
dirname <- paste(getwd(),"output_result",sep="/")
if(j==0) {
write.csv(result_model,file=paste(dirname,sub("v3",v3,(sub("j",1,(sub("j",j,"crop_ALL_8_v3_result.csv"))))),sep="/"),row.names=FALSE)

} else {
write.csv(result_model,file=paste(dirname,sub("v3",v3,(sub("j",1,(sub("j",j,"crop_j_j_v3_result.csv"))))),sep="/"),row.names=FALSE)
}


