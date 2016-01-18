#source("CovFuns.c")

dyn.load("CovFuns.so")
#######################################
#source("covStruct_Base_TensorProduct.R")
## -----------------------------------------
## Useful METHOD for prediction: covMat1Mat2
## -----------------------------------------
## --------------------
## Tensor product class
## --------------------

## covTensorProduct : separable (or tensor product) covariances, depending on 1 set of parameters
## Examples - Gaussian, Exponential, Matern with fixed nu=p+1/2, Power-Exponential
##

setClass("covTensorProduct", 		
         representation(
                        d = "integer",            	## (spatial) dimension
                        name = "character",             ## "powexp"
                        paramset.n = "integer",         ## number of parameters sets 
                        ##   gauss, exp : 1;  powexp : 2
                        var.names = "character",  	## e.g.  c("Lat", "Long") length d
			## s.d. of dor the non-nugget part of error
                        sd2 = "numeric",       		## variance (stationarity)
			## nugget part
                        known.covparam = "character",   ## known covariance parameters (except nugget): "All" or "Known"
                        nugget.flag = "logical",  	## logical : is there a nugget effect ?
                        nugget.estim = "logical", 	## logical : is it estimated (TRUE) or known ?
                        nugget = "numeric",    		## nugget (variance)
  			## total number of parameters (except sigma and nugget)
                        param.n = "integer",            ## range.n + shape.n
  			## range part 
                        range.n = "integer",            ## number of distinct range parms
                        range.names = "character",	## their name (usually "theta")
                        range.val = "numeric",          ## their values
  			## shape part, if any 
                        shape.n = "integer",            ## number of distinct shape parms
                        shape.names = "character",      ## their name ("p", "nu", "alpha", etc.)
                        shape.val = "numeric"           ## their values
                        ),
         validity = function(object) {
           
           covset <- c("gauss", "exp", "matern3_2", "matern5_2", "powexp")
           if (!is.element(object@name, covset)) {
             cat("The list of available covariance functions is:\n", covset, "\n")
             return("invalid character string for 'covtype' argument")
           }
           
           if (!identical(object@sd2, numeric(0))) {
             if (object@sd2 < 0) {
               return("The model variance should be non negative")
             }
           }
           
           if (length(object@nugget) > 1L) {
             return("'nugget' must be a single non-negative number. For heteroskedastic noisy observations, use 'noise.var' instead.")
           }
           
           if (!identical(object@nugget, numeric(0))) {
             if (object@nugget < 0) {
               return("The nugget effect should be non negative")
             }
           }
           
           if (!identical(object@range.val, numeric(0))) {
             if (min(object@range.val) < 0) {
               return("The range parameters must have positive values")
             }
             if (length(object@range.val) != object@d) {
               return("Incorrect number of range parameters")
             }
           }
         
           if (!identical(object@shape.val, numeric(0))) {
             if (min(object@shape.val) < 0) {
               return("The shape parameters must have positive values")
             }
             if (length(object@shape.val) != object@d) {
               return("Incorrect number of shape parameters")
             }
             if (identical(object@name, "powexp") && (max(object@shape.val) > 2)) {
               return("The exponents must be <= 2 for a Power-Exponential covariance")
             }
           }
           TRUE
         }
         )
#################
covparam2vect<- function(object){
  if (object@paramset.n==1) {
    param <- object@range.val
  } else param <- c(object@range.val, object@shape.val)
  return(as.numeric(param))
}
########################covMatrix

covMatrix<- function(object, X, noise.var=NULL) {
  d <- ncol(X)
  n <- nrow(X)
param <- covparam2vect(object)
out <- .C("C_covMatrix", 
            as.double(X), as.integer(n), as.integer(d), 
            as.double(param), as.double(object@sd2), as.character(object@name), 
            ans = double(n * n),
            PACKAGE="DiceKriging")
  C <- matrix(out$ans, n, n)   # covariance matrix when there is no nugget effect
 if (object@nugget.flag) {
vn <- rep(object@nugget, n)
C <- C + diag(vn, nrow = n)
  } else if (length(noise.var)>0) {
vn <- noise.var
C <- C + diag(noise.var, nrow = n)
  } else {
vn <- rep(0, n)
  }
return(list(C=C, vn=vn))	
}


###############

covMat1Mat2<- function(object, X1, X2, nugget.flag=FALSE) {
  
  # X1 : matrix n1 x d - containing training points
  # X2 : matrix n2 x d - containing test points
  
  # X1 <- as.matrix(X1)
  # X2 <- checkNames(X1, X2)
  
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  d <- ncol(X1)
  
  param <- covparam2vect(object)
  
  out <- .C("C_covMat1Mat2", 
            as.double(X1), as.integer(n1),
            as.double(X2), as.integer(n2), 
            as.integer(d),
            as.double(param), as.double(object@sd2), as.character(object@name),
            ans = double(n1 * n2), PACKAGE="DiceKriging")
  
  M <- matrix(out$ans, n1, n2)
  
  if ((!nugget.flag) | (!object@nugget.flag)) {
    return(M)
  } else {
    out <- .C("C_covMat1Mat2", 
              as.double(X1), as.integer(n1),
              as.double(X2), as.integer(n2), 
              as.integer(d),
              as.double(param), as.double(object@nugget), "whitenoise",
              ans = double(n1 * n2), PACKAGE="DiceKriging")
    N <- matrix(out$ans, n1, n2)
    return(M+N)
  }
}

##*****************************************************************************
##                        P R E D I C T  METHOD
##*****************************************************************************


emu_pred<- function(object, newdata, type,
                       se.compute = TRUE, cov.compute = FALSE, light.return = FALSE,
                       bias.correct = FALSE, ...) {
  ## newdata : n x d
  
  nugget.flag <- object@covariance@nugget.flag 
  
  X <- object@X
  y <- object@y
 newdata <- as.matrix(newdata)
 d.newdata <- ncol(newdata)
 if (!identical(d.newdata, object@d)) {
 stop("newdata must have the same numbers of columns than the experimental design")
    #}
    if (!identical(colnames(newdata), colnames(X))) {
colnames(newdata) <- colnames(X)
    }
  }
  
  T <- object@T
  z <- object@z
  M <- object@M
  
  beta <- object@trend.coef
    
  F.newdata <- model.matrix(object@trend.formula, data = data.frame(newdata))
  y.predict.trend <- F.newdata%*%beta
  
  c.newdata <- covMat1Mat2(object@covariance, X1 =as.matrix(X), X2 = as.matrix(newdata),nugget.flag = object@covariance@nugget.flag)
  
  Tinv.c.newdata <- backsolve(t(T), c.newdata, upper.tri=FALSE)
  y.predict.complement <- t(Tinv.c.newdata)%*%z
  y.predict <- y.predict.trend + y.predict.complement
  y.predict <- as.numeric(y.predict)
  
  output.list <- list()
  output.list$trend <- y.predict.trend
  output.list$mean <- y.predict
  ################Standard error of predictions
  #if ((se.compute) || (cov.compute)) {
      total.sd2 <- object@covariance@sd2
        #if (object@covariance@nugget.flag) {
        total.sd2 <- total.sd2 + object@covariance@nugget
    #}
  #}
if (se.compute) {		
 s2.predict.1 <- apply(Tinv.c.newdata, 2,crossprod)#compute c(x)'*C^(-1)*c(x)for x = newdata
        
    if (type == "SK") {
      s2.predict <- pmax(total.sd2 - s2.predict.1, 0)
      s2.predict <- as.numeric(s2.predict)
      q95 <- qnorm(0.975)
    }
    else if (type == "UK") {
      T.M <- chol(t(M)%*%M)   # equivalently : qrR <- qr.R(qr(M))
      s2.predict.mat <- backsolve(t(T.M), t(F.newdata - t(Tinv.c.newdata)%*%M) , upper.tri = FALSE)
      s2.predict.2 <- apply(s2.predict.mat, 2, crossprod)
      s2.predict <- pmax(total.sd2 - s2.predict.1 + s2.predict.2, 0)
      s2.predict <- as.numeric(s2.predict)
      #if (bias.correct) s2.predict <- s2.predict * object@n/(object@n - object@p)
      q95 <- qt(0.975, object@n - object@p)
    }
    
    lower95 <- y.predict - q95*sqrt(s2.predict)
    upper95 <- y.predict + q95*sqrt(s2.predict)
    
    output.list$sd <- sqrt(s2.predict)
    output.list$lower95 <- lower95
    output.list$upper95 <- upper95
  }
 return(output.list)
 }

####################END



