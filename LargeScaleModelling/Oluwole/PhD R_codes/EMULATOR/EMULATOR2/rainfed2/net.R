############netcdf output
crp=aperm(crop,c(1,3,2))
load("rainfed2/k3")
if(j==0) {
   w=1:8
} else {
   w=j
}
load("rainfed2/lonlat")
xvals <- Lon
yvals <- Lat
zvals=w
nz=length(zvals)
ny <- length(yvals)
nx <- length(xvals)
new_dir=getwd()
dirname <- paste(getwd(),"output_result",sep="/")
setwd(dirname)

xdim <- dim.def.ncdf("LONGITUDE","degrees",xvals )
ydim <- dim.def.ncdf("LATITUDE", "degrees",yvals )
tdim <- dim.def.ncdf("TIME","decade",zvals)
#---------
# Make var
#---------
mv <- 1.e30 # missing value
var_tropical_cereal <- var.def.ncdf("tropical_cereal","gC/m^2",dim=list(xdim,ydim,tdim),mv,longname="mean decadal change tropical_cereal")
var_pulses <- var.def.ncdf("pulses","gC/m^2",dim=list(xdim,ydim,tdim), mv,longname="mean decadal change pulses")
var_temperate_root <- var.def.ncdf("temperate_root","gC/m^2",dim=list(xdim,ydim,tdim),mv,longname="mean decadal change temperate_root")
var_soybeans <- var.def.ncdf("soybeans","gC/m^2",dim=list(xdim,ydim,tdim),mv,longname="mean decadal change soybeans")
var_ini_tropical_cereal <- var.def.ncdf("ini_tropical_cereal","gC/m^2",dim=list(xdim,ydim,tdim),mv,longname="mean decadal change ini_tropical_cereal")
var_ini_pulses <- var.def.ncdf("ini_pulses","gC/m^2",dim=list(xdim,ydim,tdim), mv,longname="mean decadal change ini_pulses")
var_ini_temperate_root <- var.def.ncdf("ini_temperate_root","gC/m^2",dim=list(xdim,ydim,tdim),mv,longname="mean decadal change ini_temperate_root")
var_ini_soybeans <- var.def.ncdf("ini_soybeans","gC/m^2",dim=list(xdim,ydim,tdim),mv,longname="mean decadal change ini_soybeans")

#---------------------
# Make new output file
#---------------------
crop=create.ncdf(paste(as.name(paste("crop",v3,sep="_")),"nc",sep="."),list(var_tropical_cereal,var_pulses,var_temperate_root,var_soybeans,var_ini_tropical_cereal,var_ini_pulses,var_ini_temperate_root,var_ini_soybeans))


#-------------------------------
dat <- array(NA,dim=c(nx*ny,nz))
dat[k3,]=crp[,,1]  ###crop(59199,8,8)
dat=array(dat,dim=c(280,720,length(w)))
dat=aperm(dat,c(2,1,3))

dat2 <- array(NA,dim=c(nx*ny,nz))
dat2[k3,]=crp[,,2]  ###crop(59199,8,8)
dat2=array(dat2,dim=c(280,720,length(w)))
dat2=aperm(dat2,c(2,1,3))

dat3 <- array(NA,dim=c(nx*ny,nz))
dat3[k3,]=crp[,,3]  ###crop(59199,8,8)
dat3=array(dat3,dim=c(280,720,length(w)))
dat3=aperm(dat3,c(2,1,3))

dat4 <- array(NA,dim=c(nx*ny,nz))
dat4[k3,]=crp[,,4]  ###crop(59199,8,8)
dat4=array(dat4,dim=c(280,720,length(w)))
dat4=aperm(dat4,c(2,1,3))

dat5 <- array(NA,dim=c(nx*ny,nz))
dat5[k3,]=crp[,,5]  ###crop(59199,8,8)
dat5=array(dat5,dim=c(280,720,length(w)))
dat5=aperm(dat5,c(2,1,3))

dat6 <- array(NA,dim=c(nx*ny,nz))
dat6[k3,]=crp[,,6]  ###crop(59199,8,8)
dat6=array(dat6,dim=c(280,720,length(w)))
dat6=aperm(dat6,c(2,1,3))

dat7 <- array(NA,dim=c(nx*ny,nz))
dat7[k3,]=crp[,,7]  ###crop(59199,8,8)
dat7=array(dat7,dim=c(280,720,length(w)))
dat7=aperm(dat7,c(2,1,3))

dat8 <- array(NA,dim=c(nx*ny,nz))
dat8[k3,]=crp[,,8]  ###crop(59199,8,8)
dat8=array(dat8,dim=c(280,720,length(w)))
dat8=aperm(dat8,c(2,1,3))


put.var.ncdf(crop,var_tropical_cereal,dat, start=c(1,1,1), count=c(nx,ny,nz))
put.var.ncdf(crop,var_pulses,dat2, start=c(1,1,1), count=c(nx,ny,nz))
put.var.ncdf(crop,var_temperate_root,dat3, start=c(1,1,1), count=c(nx,ny,nz))
put.var.ncdf(crop,var_soybeans,dat4, start=c(1,1,1), count=c(nx,ny,nz))

put.var.ncdf(crop,var_ini_tropical_cereal,dat5, start=c(1,1,1), count=c(nx,ny,nz))
put.var.ncdf(crop,var_ini_pulses,dat6, start=c(1,1,1), count=c(nx,ny,nz))
put.var.ncdf(crop,var_ini_temperate_root,dat7, start=c(1,1,1), count=c(nx,ny,nz))
put.var.ncdf(crop,var_ini_soybeans,dat8, start=c(1,1,1), count=c(nx,ny,nz))


close.ncdf(crop)
setwd(new_dir)

