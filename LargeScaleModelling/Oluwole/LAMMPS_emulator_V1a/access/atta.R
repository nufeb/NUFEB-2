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
