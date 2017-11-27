###Read in arbitrary CO2 concentration file for other RCP
new_co2=paste(climgen_path,new_co2,sep="/")
cc=as.matrix(read.table(new_co2,header=TRUE,skip=4,nrows=90,row.names=1))
cc2=array(cc,c(10,9))
cc3=colMeans(cc2)
#if(j==0){
co2=rep(cc3[1:8],each=npixel)
cco2=rep((cc3[2:9]-cc3[1:8]),each=npixel)
#}else{
#co2=rep(cc3[j],npixel)
#cco2=rep((cc3[j+1]-cc3[j]),npixel)
#}
co2_dat=cbind(cco2,co2)
############## correct next line code##############
co2_real=co2_dat
kele=other_inp[[1]][,3:4]
other_inp[[1]][,3:4]=co2_real
co2_dat=kele