### ======================================================================== ###
### log-likelihood assuming zeta_t follows a Weibull distribution [unit mean]
### ======================================================================== ###
LL_Weibull <- function(theta, n, zVec) {
    ### zVec: l* - re-scaled illiquidity process
    ### theta=[beta, gamma, shape] 
    shape <- tail(theta, 1)
    lamVec <- lambda_series(theta, n, zVec)
    lik <- (log(shape/zVec)+shape*log(gamma(1+1/shape)*zVec/lamVec)-
                (gamma(1+1/shape)*zVec/lamVec)^shape)
    return(-sum(lik))
}
### ======================================================================== ###


### ======================================================================== ###
### one-step update estimate using efficient score 
### ======================================================================== ###
MLE_Weibull_Eff <- function(theta, shape, n, zVec) {
    lamVec <- lambda_series(theta=theta, n=n, zVec=zVec)
    zetaVec <- zVec/lamVec
    
    lmdWei <- 1/gamma(1+1/shape)
    fLnDzeta <- function(x, shape) {
        return((shape-1-shape*(x/lmdWei)^shape)/x)
    } 
    fLnDshape <- function(x, shape) {
        dshape <- 1/shape+log(x/lmdWei)-log(x/lmdWei)*(x/lmdWei)^shape
        ### D lambda d shape as lambda is a function of shape [unit mean]
        dlmdshape <- (shape*((x/lmdWei)^shape-1)/lmdWei*(psigamma(x=1+1/shape,deriv=0)
                                                         /shape^2/gamma(1+1/shape)))
        return(dshape + dlmdshape)
    } 
    
    s2zeta <- -(1+zetaVec*fLnDzeta(x=zetaVec, shape=shape))
  
    dfShapeDiv <- sapply(zetaVec, function(zeta){
        fLnDshape(x=zeta, shape=shape)
    })
    
    dLambdaTheta <- lambda_dev_theta(theta=theta, n=n, zVec=zVec)
    dLambdaThetaDiv <- sapply(1:ncol(dLambdaTheta), function(j) {
        dLambdaTheta[,j]/lamVec^2
    })
    subTerm <- colMeans(dLambdaThetaDiv)/mean(1/lamVec^2)
    ltheta <- sapply(1:ncol(dLambdaTheta), function(j) {
        sweep(x=dLambdaTheta, MARGIN=2, STATS=subTerm, FUN="-")[,j]*s2zeta/lamVec
    })
    
    I2 <- mean(s2zeta^2)
    lshape <- dfShapeDiv-(mean(dfShapeDiv*s2zeta)*mean(1/lamVec)/I2/mean(1/lamVec^2)*s2zeta/lamVec)
    leta <- cbind(ltheta, matrix(lshape, ncol=1))
    dim(leta)
    Ilst <- lapply(1:nrow(leta), function(i) {
        leta[i, ] %*% t(leta[i, ])
    })
    II <- Reduce('+', Ilst)/n
    
    estEff <- c(theta, shape) + drop(solve(II)%*%colMeans(leta))
    seEff <- diag(solve(II)/sqrt(n))
    tStatsEff <- estEff/seEff
    return(list(theta=estEff, se=seEff, tStats=tStatsEff))
}
### ======================================================================== ###