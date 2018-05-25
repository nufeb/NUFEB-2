###################Relevant functions for GP regression
##  DM==Design matrix (input)
##  RB== Regressor basis
##  Obs== Output vectors
##  RM== Regressor function for matrix
##  corr_matrix== Function to calculate correlation matrix
##  obj==objective function for the optimal scale
##  opts== Function to find optimal scales (Roughness length)
##  slik== scales likelihood
##  sigma== sigmahat2
##  betahat== coefficients
##  DM_new== new Design point (crossvalidation)

Obs=Obs
DM=as.matrix(DM)
DM_new=as.matrix(DM_new)
scales_start=rep(0.5,ncol(DM))

###############
RB=function (x) 
{x <- c(1)
#x <- c(1, x)
#x=c(1, x, x[1] * x[2])
#x=c(1, x,log(x+2))
    names(x)[1] <- "const"
    return(x)}
#####Gives regressor function for matrix
RM=function (DM,RB) 
{
    out <- t(apply(DM, 1,RB))
    if (nrow(out) == 1) {
        return(t(out))
    }
    else {
        return(out)
    }
}
func=RB
############################ginv
ginv=function (X, tol =.Machine$double.eps) 
{
    if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
        stop("'X' must be a numeric or complex matrix")
    if (!is.matrix(X)) 
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X)) 
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    if (all(Positive)) 
        Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive)) 
        array(0, dim(X)[2L:1L])
    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
        t(Xsvd$u[, Positive, drop = FALSE]))
}
################################################
corr_matrix=function(DM,Obs = NULL, method = 1) 
{scales=scales_start
PDM <- diag(scales, nrow = length(scales))
    if (is.null(Obs)) {
        nully <- TRUE
        Obs <- DM
    }
    else {
        nully <- FALSE
    }	
                         jj <- function(x) {
            as.matrix(diag(as.matrix(x)))
        }
        if (nully) {
            R <- DM%*%PDM%*%t(DM)  
            S <- kronecker(jj(R), t(rep(1, nrow(DM))))
            return(exp(Re(R + t(R) - S - t(S))))
        }
        tDM <- t(DM)
        tObs <- t(Obs)
        R1 <- DM %*% PDM %*% tDM
        R2 <- Obs %*% PDM %*% t(Obs)
        S1 <- kronecker(t(jj(R1)), rep(1, nrow(Obs)))
        S2 <- kronecker(jj(R2), t(rep(1, nrow(DM))))
        a1 <- (t(tDM)%*%PDM)%*%tObs 
        a2 <- (t(tObs)%*%PDM)%*%tDM
        return(exp(Re(t(a1) + a2 - S1 - S2)))
    }  
#corr_matrix(DM=DM,Obs=DM)
###############################################
slik=function(scales,DM,Obs,RB) 
{scales=scales_start
sigma=sigmahat2(H, C, Obs)
 PDM<- diag(scales, nrow = length(scales))
    H <- RM(DM,RB)
    q <- ncol(H)
    n <- nrow(H)
    C <- corr_matrix(DM,DM)
    g <- function(C) {
        (-0.5) * sum(log(eigen(C, TRUE, TRUE)$values))
    }
    k2 <- g(C)
    k1 <- log(sigma)* (-(n - q)/2)
    k3 <- g(t(H) %*%ginv(C) %*%H)
    out <- drop(k1 + k2 + k3)
    return(out)
    }

##########################################
opts=function (DM, scales,Obs,RB,control=list(trace=1000,maxit=3)) 
{
obj <- function(scales,DM, Obs) {
-slik(scales = exp(scales),DM,Obs,RB)}
#jj <- optim(par=log(scales=scales_start),fn=obj,gr=NULL,DM,Obs)
#method=c("Nelder-Mead","Brent")
if(length(scales)==1){
method="Brent"}
else{method="Nelder-Mead"}
jj= optim(par=(scales),fn=obj,gr=NULL,DM,Obs,method=method,lower=scales_star-.5,upper=scales_star+.5)
           return(jj)}
################################################

sigmahat2=function (H, C, Obs) 
{
    C <- corr_matrix(DM=DM,Obs=DM)
    n <- nrow(C)
    q <- ncol(H)
	w=(n - q - 2)
	V=t(H)%*%ginv(C)%*%H
	W=t(Obs)%*%ginv(C)%*%Obs
	R=t(ginv(C)%*%H)
	T=t(R)%*%ginv(V)%*%R
	S=t(Obs)%*%T%*%Obs
    	(W - S)/w
}
########################################
###################################################
##################################################prediction and error of prediction
prediction=function(DM,Obs,DM_new){
scales=opts(DM, scales=scales_start,Obs,RB)$par
C <- corr_matrix(DM,DM)
PDM=diag(scales, nrow = length(scales))
DM=as.matrix(DM)
sigma=sigmahat2(H, C, Obs)
betahat= as.vector(ginv(t(H)%*%ginv(C)%*%H) %*% (t(H)%*%ginv(C)%*% Obs))
posterior_mean <- rep(NA, nrow(DM_new))
prior <- rep(NA, nrow(DM_new))
Z <- rep(NA, nrow(DM_new))
cstarx=rep(NA, nrow(DM_new))
posterior_var=rep(NA, nrow(DM_new))
for (i in 1:nrow(DM_new)) {
hh <- func(DM_new[i,])  
tx <- as.vector(corr_matrix(DM=DM, Obs= DM_new[i,, drop = FALSE], method = 1))
#prior[i] <- t(as.matrix(hh))%*%betahat 
prior[i] <- t(hh)%*%betahat
D=Obs - H %*% betahat
posterior_mean[i] = prior[i] + tx%*%ginv(C)%*%D
cstarx[i] <- 1 - t(tx)%*%ginv(C)%*%tx
posterior_var[i]=cstarx[i] + (t(hh) - t(tx)%*%ginv(C)%*%H) %*% ginv(t(H)%*%ginv(C)%*%H) %*% t(t(hh)-t(tx)%*%ginv(C)%*%H)
Z[i] <- sqrt(abs(sigma* posterior_var[i]))
result=cbind(posterior_mean,Z)}
return(result)
}
H <- RM(DM,RB)
######################################
mod=lm(Y~log(V1),data=as.data.frame(DM))
X=as.data.frame(DM_new);names(X)="V1"
lm_pred=predict(mod,newdata=X)

