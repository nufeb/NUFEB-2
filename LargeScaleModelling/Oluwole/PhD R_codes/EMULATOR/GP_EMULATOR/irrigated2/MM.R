GP=function (X, Z, constantMean = 1, nugget = NULL, nugget.known = 0, 
    min.nugget = 0, param.names = NULL, gp.names = NULL, PC.UD = NULL, 
    PC.num = NULL, PC.percent = NULL, simplex.ntries = 5, simplex.maxiter = 500, 
    simplex.reltol = 1e-08, BFGS.maxiter = 500, BFGS.tol = 0.01, 
    BFGS.h = 1e-10, seed = 0, verbose = 1, parallel = FALSE) 
{
    X = as.matrix(X)
    orig.X = X
    Z = as.matrix(Z)
    reps = NULL
    un.sum = NULL
    if (nrow(X) == 1) {
        stop("error: there must be at least two input values")
    }
    if (!is.null(param.names) && length(param.names) != dim(X)[2]) {
        stop("length of param.names must match number of columns of X")
    }
    if (!is.null(PC.percent)) 
        PC.num = numSingularValues(Z, PC.percent)
    if (!is.null(PC.num)) {
        weights = pcweights(Z, weights.num = PC.num)
        Z = t(weights$Vprime)
        if (PC.num == 1) {
            stop("error: currently not implemented for PC.num = 1 (set PC.num directly)")
        }
        PC.UD = weights$UD
    }
    if (!is.null(dim(Z))) {
        numGPs = dim(Z)[2]
    }
    else {
        numGPs = 1
        Z = matrix(Z)
    }
    if (length(seed) != 1 && length(seed) != ncol(Z)) {
        stop("error: seed length must equal dimensions of Z, or 1")
    }
    if (is.null(gp.names)) {
        gp.names = paste("gp #", 1:numGPs)
    }
    if (length(gp.names) != numGPs) {
        stop("length of gp.names must match number of GPs")
    }
    if (dim(X)[1] != dim(Z)[1]) {
        stop("error: X matrix and output matrix must have same number of rows; \n\tfor PC weights, make sure each column of Z is a computer model run")
    }
    if (min(diag(var(Z))) <= 0) {
        stop("error: Z does not vary for at least one column")
    }
    if (dim(Z)[2] > 1 && !is.null(nugget) && length(nugget) > 
        1) {
        stop("error: cannot handle nugget matrix with multiple field observations; try fitting 1 gp at a time")
    }
    if (!is.null(nugget) && length(nugget) > dim(Z)[1]) {
        stop("length of nugget matrix must match number of observations")
    }
    if (min(apply(X, 2, max) - apply(X, 2, min)) == 0) {
        stop("error: at least one parameter does not vary in design")
    }
    if (constantMean != 1) {
        ones = rep(1, dim(X)[1])
        dX = cbind(ones, X)
		#library(MASS)
        t = try(ginv(t(dX) %*% dX), TRUE)
        if (class(t) == "try-error") {
            stop("error: design matrix with intercept is not full rank; set constantMean = 1")
        }
    }
    if (nugget.known == 1) {
        if (!is.null(nugget) && length(nugget) > 1 && length(nugget) != 
            nrow(X)) {
            stop("error: for nugget.known = 1, length of nugget matrix must be nrow(X) or 1")
        }
        un.sum = uniqueSummary(X, Z)
        orig.X = X
        X = un.sum$uniqueX
        Z = un.sum$uniqueMeans
        reps = un.sum$reps
    }
    if (!is.null(nugget) && length(nugget) == 1 && nugget <= 
        0) {
        if (anyReps(X)) 
            stop("at least 2 or more inputs are identical...must use nugget!")
    }
    if (is.null(nugget)) {
        if (!anyReps(X)) {
            if (verbose > 0) 
                cat("no reps detected - nugget will not be estimated\n")
        }
        else {
            if (verbose > 0 && nugget.known == 0) 
                cat("reps detected - nugget will be estimated\n")
        }
    }
    numEstimates = dim(X)[2] + 1
    if (anyReps(X) || !is.null(nugget) || length(nugget) > 1) 
        numEstimates = numEstimates + 1
    if (constantMean == 1) {
        numEstimates = numEstimates + 1
    }
    else {
        numEstimates = numEstimates + dim(X)[2] + 1
    }
    estimates = rep(0, numEstimates)
    if (verbose > 0) 
        cat("\n")
    l = NULL
    Zlist = list()
    if (length(seed) == 1) 
        seed = rep(seed, ncol(Z))
    for (i in 1:ncol(Z)) {
        Zlist[[i]] = list()
        Zlist[[i]]$Z = matrix(Z[, i])
        Zlist[[i]]$i = i
        Zlist[[i]]$seed = seed[i]
    }
    if (!parallel) {
        l = lapply(Zlist,GP2, XX = X, orig.XX = orig.X, nugget = nugget, 
            nugget.known = nugget.known, reps = reps, un.sum = un.sum, 
            numEstimates = numEstimates, constantMean = constantMean, 
            simplex.ntries = simplex.ntries, simplex.maxiter = simplex.maxiter, 
            simplex.reltol = simplex.reltol, BFGS.maxiter = BFGS.maxiter, 
            BFGS.tol = BFGS.tol, BFGS.h = BFGS.h, min.nugget = min.nugget, 
            verbose = verbose, param.names = param.names)
    }
    if (parallel) {
        cat("fitting ")
        cat(length(Zlist))
        cat(" GPs in parallel mode...\n")
        l = sfLapply(Zlist,GP2, XX = X, orig.XX = orig.X, 
            nugget = nugget, nugget.known = nugget.known, reps = reps, 
            un.sum = un.sum, numEstimates = numEstimates, constantMean = constantMean, 
            simplex.ntries = simplex.ntries, simplex.maxiter = simplex.maxiter, 
            simplex.reltol = simplex.reltol, BFGS.maxiter = BFGS.maxiter, 
            BFGS.tol = BFGS.tol, BFGS.h = BFGS.h, min.nugget = min.nugget, 
            verbose = verbose, param.names = param.names)
        cat("...done\n")
    }
    if (numGPs == 1) {
        return(l[[1]])
    }
    return(gp.list(l, UD = PC.UD, param.names = param.names, 
        gp.names = gp.names))
	#detach("package:MASS", unload=TRUE)
}
