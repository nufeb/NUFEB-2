##############Read and process LAMMPS input files
step=2000
text <- readLines("Input/Inputscript.lammps",encoding="UTF-8")
ppara=c("variable KsHET equal","variable Ko2HET equal","variable Kno2HET equal","variable Kno3HET equal","variable Knh4AOB equal","variable Ko2AOB equal","variable Kno2NOB equal","variable Ko2NOB equal","variable MumHET equal","variable MumAOB equal","variable MumNOB equal","variable etaHET equal","variable bHET equal","variable bAOB equal","variable bNOB equal","variable bEPS equal","variable YEPS equal","variable YHET equal","variable EPSratio equal","variable factor equal")


id=c(rep(NA,length(ppara)))
for(i in 1:length(ppara)){
#id[i]=grep(ppara[i],text,value=FALSE,fixed=TRUE)##may be prone to error watch
id[i]=grep(ppara[i],text,value=FALSE,fixed=TRUE)[1]##corrected
}
inp1=text[id]
inp2=c(rep(NA,length(ppara)))
for(i in 1:length(ppara)){
inp2[i]=as.numeric(gsub(ppara[i], "",inp1[i],fixed=TRUE))
}
#olu2=as.numeric(gsub("[^[:digit:].]", "",olu))
abc=grep("run",text,value=FALSE,fixed=TRUE)
runtime=as.numeric(gsub("run", "",text[abc],fixed=TRUE))
######################################################################
para1=c("1 1 1.0839e-6 150 4e-05 3.5e-05 8.6e-07 1.4307e-6 8.00e-02 10.00e-03 9.00e-02 8.000e-03 1.0000e-05","2 2 1.2839e-6 150 4e-05 4e-05 9.26e-07 1.2839e-6 8.00e-02 10.00e-03 9.00e-02 8.000e-03 1.0000e-05","3 3 1.0283e-6 150 4e-05 4.5e-05 8.6e-07 1.0283e-6 8.00e-02 10.00e-03 9.00e-02 8.000e-03 1.0000e-05","4 4 1.0539e-6 150 3.5e-05 4e-05 9.6e-07 1.0539e-6 8.00e-02 10.00e-03 9.00e-02 8.000e-03 1.0000e-05","5 5 1.4839e-6 150 4.5e-05 4e-05 10.6e-07 1.4839e-6 8.00e-02 10.00e-03 9.00e-02 8.000e-03 1.0000e-05")

text1 <- readLines(paste("Input/flat_surface", ".in",sep=""),encoding="UTF-8")
index=c(rep(NA,length(para1)))
for(i in 1:length(para1)){
index[i]=grep(para1[i],text1,value=FALSE,fixed=TRUE)
}
inp3=text1[index[1]]
inp4=strsplit(inp3, " ")[[1]]
inp4=as.numeric(inp4[-(1:13)]);inpp=inp4[c(1,4,5,2,3)]
para=c("KsHET","Ko2HET","Kno2HET","Kno3HET","Knh4AOB","Ko2AOB","Kno2NOB","Ko2NOB","MumHET","MumAOB","MumNOB","etaHET","bHET","bAOB","bNOB","bEPS","YEPS","YHET","EPSratio","factor")
namm=c("sub","no2","no3","o2","nh4")
DM=c(inp2,inpp)
ntime=round(runtime/step)
time=cap=seq(0,(runtime-step),step)
j=length(time)
DM=matrix(rep(DM,ntime),nrow=ntime,byrow=TRUE)
DM_new=as.data.frame(cbind(DM,time[1:j]))
names(DM_new)=nam=c(para,namm,"time")
##############################################end



