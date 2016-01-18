############COVMAT1MAT2### cross covariance using C-functions
dyn.load("access/CovFuns.so")###load .c function into R
#model=model
##
covparam2vect<- function(model){
  if (model@paramset.n==1) {
    param <- model@range.val
  } else param <- c(model@range.val, model@shape.val)
  return(as.numeric(param))
}

################CovMat
source("access/atta.R")
covMatrix<- function(model, X, noise.var=NULL) {
  d <- ncol(X)
  n <- nrow(X)
param <- covparam2vect(model)
out <- .C("C_covMatrix", 
            as.double(X), as.integer(n), as.integer(d), 
            as.double(param), as.double(model@sd2), as.character(model@name), 
            ans = double(n * n))
  C <- matrix(out$ans, n, n)   # covariance matrix when there is no nugget effect
 if (model@nugget.flag) {
vn <- rep(model@nugget, n)
C <- C + diag(vn, nrow = n)
  } else if (length(noise.var)>0) {
vn <- noise.var
C <- C + diag(noise.var, nrow = n)
  } else {
vn <- rep(0, n)
  }
return(list(C=C, vn=vn))
}
################
covMat1Mat2<- function(model, X1, X2, nugget.flag=FALSE) {
  
  # X1 : matrix n1 x d - containing training points
  # X2 : matrix n2 x d - containing test points
  
  # X1 <- as.matrix(X1)
  # X2 <- checkNames(X1, X2)
  
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  d <- ncol(X1)
  
  param <- covparam2vect(model)
  
  out <- .C("C_covMat1Mat2", 
            as.double(X1), as.integer(n1),
            as.double(X2), as.integer(n2), 
            as.integer(d),
            as.double(param), as.double(model@sd2), as.character(model@name),
            ans = double(n1 * n2))
  
  M <- matrix(out$ans, n1, n2)
  
  if ((!nugget.flag) | (!model@nugget.flag)) {
    return(M)
  } else {
    out <- .C("C_covMat1Mat2", 
              as.double(X1), as.integer(n1),
              as.double(X2), as.integer(n2), 
              as.integer(d),
              as.double(param), as.double(model@nugget), "whitenoise",
              ans = double(n1 * n2))
    N <- matrix(out$ans, n1, n2)
    return(M+N)
  }
 }
