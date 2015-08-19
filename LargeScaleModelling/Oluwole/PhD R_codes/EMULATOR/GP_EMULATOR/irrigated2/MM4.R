CV2=function (gp, predictObserved = TRUE, verbose = FALSE) 
{
    Var =ginv(gp$invVarMatrix)
    Z.pred = matrix(0, ncol = 1, nrow = gp$numObs)
    pred.error = matrix(1, ncol = 1, nrow = gp$numObs)
    complete = matrix(FALSE, ncol = 1, nrow = gp$numObs)
    i = 1
    for (i in 1:gp$numObs) {
        if (verbose) {
            s = paste("cross validate at index # ", i)
            cat(s)
            cat("\n")
        }
        if (complete[i]) 
            next
        index = matrix(TRUE, gp$numObs)
        for (j in 1:gp$numObs) {
            index[j] = !all(gp$X[i, ] == gp$X[j, ])
        }
        newX = gp$X[index, ]
        if (!is.matrix(newX)) 
            newX = as.matrix(newX)
        p = predictNewYCV2(gp, newX, gp$Z[index], gp$mu[index], 
            gp$X[i, ], Var[index, index])
        n = 0
        if (predictObserved) {
            n = gp$nugget
            if (length(gp$nugget) > 1) {
                n = gp$nugget[i]
            }
        }
        v2 = calcPredictionErrorCV2(gp, newX, gp$X[i, ], Var[index, 
            index], n)
        Z.pred[!index] = p
        pred.error[!index] = v2
        complete = complete | !index
    }
    pred.error[pred.error < 0] = 0
    return(cbind(Z.pred, pred.error))
}
########
predictNewYCV2=function (gp, X, Z, muMatrix, newTheta, varMatrix) 
{
    r = calcCorOneObs(X, gp$beta, gp$a, newTheta)
    newY = predictMu(gp, newTheta) + gp$sig2 * r %*% ginv(varMatrix) %*% 
        (Z - muMatrix)
    return(newY)
}
#######
calcLogLikeManual=function (gp) 
{
    N = dim(as.matrix(gp$Z))[1]
    ans = -N/2 * log(2 * pi) - 0.5 * determinant(ginv(gp$invVarMatrix), 
        logarithm = TRUE)$modulus[1]
    ans = ans - 0.5 * t(gp$Z - gp$mu) %*% gp$invVarMatrix %*% 
        (gp$Z - gp$mu)
    return(ans)
}
##########
calcPredictionErrorCV2=function (gp, X, newTheta, varMatrix, nugget) 
{
    r = calcCorOneObs(X, gp$beta, gp$a, newTheta)
    return((gp$sig2 + nugget) - gp$sig2 * r %*% ginv(varMatrix) %*% 
        t(r) * gp$sig2)
}