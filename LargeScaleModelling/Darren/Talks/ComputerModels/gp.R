# gp.R
# script to demo simple GP emulation for 1 input and 1 output

# setup
Grid=seq(0,1,0.01)
K=function(d,scale=20,lengthScale=0.4)
{
  (scale^2)*exp(-(d/lengthScale)^2)
}

# function definitions
GPAdjust=function(Grid,ObsLoc,Obs) {
  lG=length(Grid)
  lO=length(ObsLoc)
  All=c(Grid,ObsLoc)
  Dist=as.matrix(dist(All))
  Var=K(Dist)
  VG=Var[1:lG,1:lG]
  VO=Var[(lG+1):(lG+lO),(lG+1):(lG+lO)]
  CovGO=Var[1:lG,(lG+1):(lG+lO)]
  CovOG=t(CovGO)
  ViC=solve(VO,CovOG)
  AdjEx=t(ViC)%*%Obs
  AdjVar=VG-CovGO%*%ViC
  list(Ex=as.vector(AdjEx),Var=diag(AdjVar))
}

plotAdjust=function(Grid,ObsLoc,Obs,...)
{
  plot(ObsLoc,Obs,xlim=range(Grid),pch=19,...)
  Adj=GPAdjust(Grid,ObsLoc,Obs)
  Adj$Var[Adj$Var<0]=0 # numerical stability fudge
  Upper=Adj$Ex+2*sqrt(Adj$Var)
  Lower=Adj$Ex-2*sqrt(Adj$Var)
  polygon(c(Grid,rev(Grid),Grid[1]),c(Upper,rev(Lower),Upper[1]),col="gray",border=NA)
  lines(Grid,Adj$Ex,lwd=2)
  points(ObsLoc,Obs,pch=19)
}

# calculations
ObsLoc=c(0.3,0.6,0.4,0.9,0.1)
Obs=c(10,11,8,12,12)
for (i in 1:length(Obs)) {
  plotAdjust(Grid,ObsLoc[1:i],Obs[1:i],ylim=c(0,20))
}



# eof


