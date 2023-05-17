### ======================================================================== ###
### GMM for lambda parameter estimation
### based on the first moment restriction
### ======================================================================== ###


### ======================================================================== ###
### instruments and conditions series construction for GMM
### ======================================================================== ###
instruments_lags <- function(zVec, lamVec, n) {
    numLag <- 2   ### number of lags as instruments
    zt <- cbind(rep(1, n-numLag),
                zVec[2:(n-1)],
                zVec[1:(n-2)],
                zVec[2:(n-1)]/lamVec[2:(n-1)],
                zVec[1:(n-2)]/lamVec[1:(n-2)])   ### instruments w_t(1, z_t-1)
    return(zt)
}
### ======================================================================== ###
rho_GMM <- function(theta, n, zVec, retE=FALSE) {
    lamVec <- lambda_series(theta, n=n, zVec=zVec)
    zt <- instruments_lags(zVec=zVec, lamVec=lamVec, n=n)
    numLag <- n-nrow(zt)   
    ehat <- zVec[-(1:numLag)]/lamVec[-(1:numLag)]-1   
    rho <- ehat*zt
    
    if (retE) return(ehat)
    
    return(rho)
}
### ======================================================================== ###


### ======================================================================== ###
### Objective function GMM
### ======================================================================== ###
Q_GMM <- function(theta, n, zVec, lmax=5) {
    rho <- rho_GMM(theta=theta, n=n, zVec=zVec)
    W <- diag(ncol(rho))   ### identity matrix
    M <- matrix(colMeans(rho), ncol=1)
    Q <- 1/2*t(M)%*%solve(W)%*%M
    return(as.numeric(Q))
}    
### ======================================================================== ###


### ======================================================================== ###
### Newey-West Covariance matrix estimation    
### ======================================================================== ###
cov_NW <- function(rho, lmax=5) {
    W   <- t(rho) %*% rho
    tau <- 1
    while (tau <= lmax) {
        Wtau <- t(rho[(tau+1):nrow(rho),]) %*% rho[1:(nrow(rho)-tau),]
        W <- W + (1.0-tau/(lmax+1))*(Wtau + t(Wtau))
        tau  <- tau + 1      
    }  
    n <- nrow(rho)
    W <- W/n  
    return(W)
}
### ======================================================================== ###


### ======================================================================== ###
fgrad <- function(theta, n, zVec, retE) {
    rho <- rho_GMM(theta=theta, n=n, zVec=zVec, retE=FALSE)
    return(colMeans(rho))
}
stdError_GMM <- function(theta, n, zVec, lmax=10) {
    dg <- numgrad(fgrad, theta, n=n, zVec=zVec, retE=FALSE)
    rho <- rho_GMM(theta=theta, n=n, zVec=zVec, retE=FALSE)
    Omega <- cov_NW(rho=rho, lmax=lmax)
    WI <- diag(ncol(rho)) 
    V  <- solve(t(dg) %*% WI %*% dg)
    VMid <- t(dg) %*% WI %*% Omega %*% WI %*% dg
    stdError <- sqrt(diag((V %*% VMid  %*% V)/n))
}
### ======================================================================== ###

