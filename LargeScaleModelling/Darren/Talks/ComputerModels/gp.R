# gp.R
# script to demo simple GP emulation for 1 input and 1 output

# setup
Grid=seq(0,1,0.01)

# GP covariance kernel
K=function(d,scale=50,lengthScale=0.4)
{
  (scale^2)*exp(-(d/lengthScale)^2)
}

# fit a mean zero GP
GPAdjust=function(Grid,ObsLoc,Obs) 
{
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
  AdjVar=AdjVar+diag(rep((1e-6)*max(AdjVar),lG)) # numerical stability fudge
  list(Ex=as.vector(AdjEx),Var=AdjVar)
}

# fit a GP to the residuals from a linear regression fit
GPAdjustLR=function(Grid,ObsLoc,Obs) 
{
  Mod=lm(Obs~ObsLoc)
  Pred=predict(Mod,newdata=data.frame(ObsLoc=Grid))
  Raw=GPAdjust(Grid,ObsLoc,Mod$residuals)
  list(Ex=Raw$Ex+Pred,Var=Raw$Var)
}

# produce a plot illustrating the GP fit
plotAdjust=function(Grid,ObsLoc,Obs,...)
{
  plot(ObsLoc,Obs,xlim=range(Grid),pch=19,...)
  Adj=GPAdjustLR(Grid,ObsLoc,Obs)
  Var=diag(Adj$Var)
  Upper=Adj$Ex+2*sqrt(Var)
  Lower=Adj$Ex-2*sqrt(Var)
  polygon(c(Grid,rev(Grid),Grid[1]),c(Upper,rev(Lower),Upper[1]),col="gray",border=NA)
  lines(Grid,Adj$Ex,lwd=2)
  points(ObsLoc,Obs,pch=19)
}

# generate a sample from an adjusted object
sampleGP=function(Adj)
{
  L=t(chol(Adj$Var))
  Adj$Ex + L%*%rnorm(length(Adj$Ex))
}

# Now actually do some stuff...
ObsLoc=c(0.3,0.6,0.4,0.9,0.1,0.8)
Obs=c(10,11,8,12,12,13)

# Show data progressively
for (i in 1:length(Obs))
  plot(ObsLoc[1:i],Obs[1:i],xlim=range(Grid),pch=19,ylim=c(0,20))

# Show a progressive GP fit
for (i in 1:length(Obs))
  plotAdjust(Grid,ObsLoc[1:i],Obs[1:i],ylim=c(0,20))

# Now show some samples from a GP

for (i in 10^(0:2)) {
  plot(ObsLoc,Obs,xlim=range(Grid),pch=19,ylim=c(0,20))
  Adj=GPAdjustLR(Grid,ObsLoc,Obs)
  for (j in 1:i)
    lines(Grid,sampleGP(Adj),col="gray")
  points(ObsLoc,Obs,pch=19)
  }

# eof


