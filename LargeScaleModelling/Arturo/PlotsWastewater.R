# reading properties of cells
Bac_x = scan(file="H:/Arturo A-A/Master Project/Code/R/floc2.4PDEO2/Bac_x.txt", what = double(0), n=-1, sep ="", dec = ".", skip=0,na.strings="NA")
Bac_y = scan(file="H:/Arturo A-A/Master Project/Code/R/floc2.4PDEO2/Bac_y.txt", what = double(0), n=-1, sep ="", dec = ".", skip=0,na.strings="NA")
Bac_r = scan(file="H:/Arturo A-A/Master Project/Code/R/floc2.4PDEO2/Bac_r.txt", what = double(0), n=-1, sep ="", dec = ".", skip=0,na.strings="NA")
Bac_s = scan(file="H:/Arturo A-A/Master Project/Code/R/floc2.4PDEO2/Bac_s.txt", what = double(0), n=-1, sep ="", dec = ".", skip=0,na.strings="NA")

# reading number of flocs
nFlocs = scan(file="H:/Arturo A-A/Master Project/Code/R/floc2.4PDEO2/3.0/nFlocs.txt", what = double(0), n=-1, sep ="", dec = ".", skip=0,na.strings="NA")

# reading mass of particulate components in the bioreactor
mass_reactorHET = scan(file="H:/Arturo A-A/Master Project/Code/R/mass_reactorHET.txt", what = double(0), n=-1, sep ="", dec = ".", skip=0,na.strings="NA")
mass_reactorAOB = scan(file="H:/Arturo A-A/Master Project/Code/R/mass_reactorAOB.txt", what = double(0), n=-1, sep ="", dec = ".", skip=0,na.strings="NA")
mass_reactorNOB = scan(file="H:/Arturo A-A/Master Project/Code/R/mass_reactorNOB.txt", what = double(0), n=-1, sep ="", dec = ".", skip=0,na.strings="NA")
mass_reactorEPS = scan(file="H:/Arturo A-A/Master Project/Code/R/mass_reactorEPS.txt", what = double(0), n=-1, sep ="", dec = ".", skip=0,na.strings="NA")
mass_reactorDEAD = scan(file="H:/Arturo A-A/Master Project/Code/R/mass_reactorDEAD.txt", what = double(0), n=-1, sep ="", dec = ".", skip=0,na.strings="NA")

# reading mass of particulate components in the bioreactor
soluble_reactorNH4 = scan(file="H:/Arturo A-A/Master Project/Code/R/soluble_reactorNH4.txt", what = double(0), n=-1, sep ="", dec = ".", skip=0,na.strings="NA")
soluble_reactorNO2 = scan(file="H:/Arturo A-A/Master Project/Code/R/soluble_reactorNO2.txt", what = double(0), n=-1, sep ="", dec = ".", skip=0,na.strings="NA")
soluble_reactorNO3 = scan(file="H:/Arturo A-A/Master Project/Code/R/soluble_reactorNO3.txt", what = double(0), n=-1, sep ="", dec = ".", skip=0,na.strings="NA")
soluble_reactorO2 = scan(file="H:/Arturo A-A/Master Project/Code/R/soluble_reactorO2.txt", what = double(0), n=-1, sep ="", dec = ".", skip=0,na.strings="NA")
soluble_reactorS = scan(file="H:/Arturo A-A/Master Project/Code/R/soluble_reactorS.txt", what = double(0), n=-1, sep ="", dec = ".", skip=0,na.strings="NA")

# calculating the length of the vectors ( the PDE version crashes before 100 iterations)
n=length(mass_reactorHET)/10

# ploting number of flocs
plot(x=seq(0,n,0.1),y=nFlocs,type="l",ylab="Number of flocs", xlab="Days", cex.lab=1.5,cex.axis=2,lwd=2,col="blue")

# ploting particulate components
plot(x=seq(0,n-0.1,0.1),y=mass_reactorHET,type="l",ylab="Concentration", xlab="Days", cex.lab=2,cex.axis=2,lwd=2,col="black", main="Particulate concentration 3.5 days PDE",cex.main=2)
lines(x=seq(0,n-0.1,0.1),y=mass_reactorAOB,col="red",lwd=2)
lines(x=seq(0,n-0.1,0.1),y=mass_reactorNOB,col="green",lwd=2)
lines(x=seq(0,n-0.1,0.1),y=mass_reactorEPS,col="blue",lwd=2)
lines(x=seq(0,n-0.1,0.1),y=mass_reactorDEAD,col="lightblue",lwd=2)

legend("topleft",c("HET","AOB","NOB","EPS","DEAD"),lty=c(1,1,1,1,1),lwd=c(4,4,4,4,4),col=c("black","red","green","blue","lightblue"),pt.cex=1, cex=2)

# ploting soluble components
plot(x=seq(0,n,0.1),y=soluble_reactorO2,ylim=c(0,2),type="l",ylab="Concentration", xlab="Days", cex.lab=2,cex.axis=2,lwd=2,col="green", main="Soluble concentration 3.5 days PDE", cex.main=2)
lines(x=seq(0,n,0.1),y=soluble_reactorNO2,col="black",lwd=2)
lines(x=seq(0,n,0.1),y=soluble_reactorNO3,col="red",lwd=2)
lines(x=seq(0,n,0.1),y=soluble_reactorNH4,col="blue",lwd=2)
lines(x=seq(0,n,0.1),y=soluble_reactorS,col="lightblue",lwd=2)

legend("topleft",c("NO2","NO3","O2","NH4","S"),lty=c(1,1,1,1,1),lwd=c(4,4,4,4,4),col=c("black","red","green","blue","lightblue"),pt.cex=1, cex=2)


# we add 1 so HET is not white color.
Bac_s = Bac_s + 1

# for a 50x50 floc (with PDE)
 size_n_color(Bac_x,Bac_y,Bac_r,NA,
  col=Bac_s,
  xlab="X",ylab="Y",
  main="Cells floc 2.4 days PDE",
  xat=seq(-2.5e-05,6.5e-05,by=0.00001),yat=seq(0,9.5e-05,by=0.00001),cex.main=2) 

# for a 450x450 floc
 size_n_color(Bac_x,Bac_y,Bac_r,NA,xlim=c(0.00015,0.00035),
  col=Bac_s,
  xlab="X",ylab="Y",
  main="Cells floc 10 days",
  xat=seq(0.00015,0.00035,by=0.00001),yat=seq(0.00005,0.00033,by=0.00001)) 

