### ======================================================================== ###
### updating equation of lambda_t process
### ======================================================================== ###
lambda_series <- function(theta, n, zVec) {
    lamVec <- numeric(n)+1   
    
    for (i in 2:n) {
        lamVec[i] <- 1-theta[1]-theta[2]+theta[1]*lamVec[i-1]+theta[2]*zVec[i-1]
    } 
   
    return(lamVec)
}
### ======================================================================== ###


### ======================================================================== ###
### derivative of lambda wrt theta
### ======================================================================== ###
lambda_dev_theta <- function(theta, n, zVec) {
    lamVec <- numeric(n)+1
    dThetaLst <- list(dBeta=numeric(n), dGamma=numeric(n))
    for (i in 2:n) {
        lamVec[i] <- 1-theta[1]-theta[2]+theta[1]*lamVec[i-1]+theta[2]*zVec[i-1]
        dThetaLst$dBeta[i] <- -1+theta[1]*dThetaLst$dBeta[i-1]+lamVec[i-1]
        dThetaLst$dGamma[i] <- -1+theta[1]*dThetaLst$dGamma[i-1]+zVec[i-1]
    } 
    return(do.call(cbind, dThetaLst))
}
### ======================================================================== ###