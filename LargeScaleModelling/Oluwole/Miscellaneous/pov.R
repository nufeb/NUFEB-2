#########
dirname=getwd()
index=grep("ITEM: ATOMS id type diameter mass",readLines("snapshot.bubblemd"),value=FALSE,fixed=TRUE)
ind=length(readLines("snapshot.bubblemd"))
cnam=c("type","diameter","mass","x","y","z","sub","no2","no3","o2","nh4")
dd=list()
index=c(index,ind+9)
for (i in 1:(length(index)-1)){
dd[i]=list(read.csv("snapshot.bubblemd",sep="",skip=index[i]-1,header=TRUE,colClasses=numeric(),nrows=index[i+1]-index[i]-9)[,2:(length(cnam)+1)]);names(dd[[i]])=cnam
}
save(dd,file=paste(dirname,"dd",sep="/"))
###
source("povray.r")
#gaja="C:\\Users\\olu\\Desktop\\NEW_LAMMPS_RESULTS"
gaja=getwd()
for (i in 1:length(dd)) {
dat=dd[[i]][,c(1:6)]
da=dat[,c(2,4:6)]*10^6
dat[,c(2,4:6)]=da
dat=as.matrix(dat)
###########################
source("povray.r")
sc2=Scene()
# Place some lights and cameras
l1=Light(c(52,53,-1),Colour(1,1,1))###c(+2,+4,-3) unit
#c1=Camera(c(50,51,-1),c(50,49,6))
c1=Camera(c(50,51,-99),c(50,49,6))
sca=scale(c(0.24,0.24,0.24))
# Define two colours, and a simple reflective texture
mycol=c(Colour(.2,.2,.8),Colour(.8,.2,.2),Colour(.2,.8,.2),Colour(.5,.5,.5),Colour(.17,.345,.5656))
txt=Texture("finish{reflection 0.2 specular 0.3 ambient 0.42}")
Tp=dat[,1];X=dat[,4];Y=dat[,5];Z=dat[,6];D=dat[,2]
for(l in 1:length(Tp)){
if(Tp[l]==1){
coll=mycol[[1]]}
if(Tp[l]==2){
coll=mycol[[2]]}
if(Tp[l]==3){
coll=mycol[[3]]}
if(Tp[l]==4){
coll=mycol[[4]]}
if(Tp[l]==5){
coll=mycol[[5]]}
sph= Sphere(centre=c(X[l],Y[l],Z[l]),radius=D[l]/2,col=coll,tex=txt)
sc2$push(sph)
}
sc2$push(l1)
#sc2$push(l2)
sc2$push(c1)
sc2$save()
sc2$render(aa=T)
}
####################################PoVRAY###########################
