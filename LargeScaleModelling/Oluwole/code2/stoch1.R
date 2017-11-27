####in diferent folders
nvar3=5##no of initial condition to be varied
dirname=getwd()
setwd(paste(getwd(),"/",gsub("rep",rep,"inputrep"),sep=""))
#val1=format(c(8.00e-02,10.00e-03,9.00e-02,8.000e-03,1.0000e-05,scientific=TRUE))#c("sub","no2","no3","o2","nh4") 
val1=c(8.00e-02,10.00e-03,9.00e-02,8.000e-03,1.0000e-05)
val2=("8.00e-02 10.00e-03 9.00e-02 8.000e-03 1.0000e-05")

#para1=c("1 1 1.0839e-6 150 4e-05 3.5e-05 4e-05 1.4307e-6 8.00e-02 10.00e-03 9.00e-02 8.000e-03 1.0000e-05","2 2 1.2839e-6 150 4e-05 4e-05 4e-05 1.2839e-6 8.00e-02 10.00e-03 9.00e-02 8.000e-03 1.0000e-05","3 3 1.0283e-6 150 4e-05 4.5e-05 4e-05 1.0283e-6 8.00e-02 10.00e-03 9.00e-02 8.000e-03 1.0000e-05","4 4 1.0539e-6 150 3.5e-05 4e-05 4e-05 1.0539e-6 8.00e-02 10.00e-03 9.00e-02 8.000e-03 1.0000e-05","5 5 1.4839e-6 150 4.5e-05 4e-05 4e-05 1.4839e-6 8.00e-02 10.00e-03 9.00e-02 8.000e-03 1.0000e-05") 
para1=c("1 1 1.0839e-6 150 4e-05 3.5e-05 8.6e-07 1.4307e-6 8.00e-02 10.00e-03 9.00e-02 8.000e-03 1.0000e-05","2 2 1.2839e-6 150 4e-05 4e-05 9.26e-07 1.2839e-6 8.00e-02 10.00e-03 9.00e-02 8.000e-03 1.0000e-05","3 3 1.0283e-6 150 4e-05 4.5e-05 8.6e-07 1.0283e-6 8.00e-02 10.00e-03 9.00e-02 8.000e-03 1.0000e-05","4 4 1.0539e-6 150 3.5e-05 4e-05 9.6e-07 1.0539e-6 8.00e-02 10.00e-03 9.00e-02 8.000e-03 1.0000e-05","5 5 1.4839e-6 150 4.5e-05 4e-05 10.6e-07 1.4839e-6 8.00e-02 10.00e-03 9.00e-02 8.000e-03 1.0000e-05")
set.seed(3)
hyper2=maximinLHS(npoints,nvar3, 2)
design1= matrix(NA,nrow=npoints, ncol=nvar3)
for(j in 1:nvar3){
design1[,j]=qunif(hyper2[,j], min=0.8*val1[j], max=1.2*val1[j])
}
gah=design1[,c(1,4,5,2,3)]
design1=gah
text1 <- readLines(paste("flat_surface", ".in",sep=""),encoding="UTF-8")
save(design1,file=paste(dirname,"design1",sep="/"))

des=rep(NA,npoints)
for(i in 1:npoints){
des[i]=paste(design1[i,1],design1[i,2],design1[i,3],design1[i,4],design1[i,5],sep=" ")
}
#########################################finding index
index=c(rep(NA,length(para1)))
for(i in 1:length(para1)){
index[i]=grep(para1[i],text1,value=FALSE,fixed=TRUE)
}
for(j in 1:nvar3){
#text1[index[j]]=gsub(design1[re3,j] ,design1[(re3)+1,j], text1[index[j]])
text1[index[j]]=gsub(val2,des[rep],text1[index[j]])
}
writeLines(text1,paste("flat_surface", ".in",sep=""))
setwd("..")
#######################end


