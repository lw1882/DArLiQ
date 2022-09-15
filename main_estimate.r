### ======================================================================== ###
rm(list=ls())
### ======================================================================== ###
library("Rsolnp")
library("xts")
### ======================================================================== ###
source("main_data.r")
source("utils.r")
source("g_trend_LL.r")
source("lambda_update.r")
source("lambda_GMM.r")
source("lambda_MLE.r")
### ======================================================================== ###



### ======================================================================== ###
### Step 1: estimate long-run trend function g(t/T) -- local linear
### ======================================================================== ###
tV <- (1:n)/n   # tVec=(1:n)/n
gKS <- g_trend_LL(tVec=tV, lVec=liquidity)   # estimated trend function (LL)
zV <- liquidity/gKS$y   # l*: re-scaled illiquidity
### ======================================================================== ###
### plot illiquidity series and the trend function
liquidityTS <- as.xts(cbind(gKS$y, liquidity), 
                      order.by=dates[2:(n+1)])
colnames(liquidityTS) <- c("Trend", "Illiquidity")
cols <- c("darkred", "darkblue")
lwds <- c(3, 1.5)
plot(liquidityTS*1e10, col=cols,
     main="Illiquidity (X 10^10)", ylab="",
     lwd=lwds, cex.lab=1.25, cex.axis=1.25)
addLegend("topright", legend.names=c("Trend", "Illiquidity"),
          lty=c(1, 1), lwd=lwds, cex=1.5, col=cols)
plot(log(liquidityTS), col=cols,
     main="Log of Illiquidity", ylab="",
     lwd=lwds, cex.lab=1.25, cex.axis=1.25)
addLegend("topright", legend.names=c("Trend", "Illiquidity"),
          lty=c(1, 1), lwd=lwds, cex=1.5, col=cols)
### ======================================================================== ###


### ======================================================================== ###
### Step 2: estimate the dynamic parameters of lambda_t process
### ======================================================================== ###
theta0 <- c(0.9, 0.05, 1.3)
inequal <- function(theta, n, zVec) {
    return(theta[1]+theta[2])
}
## ======================================================================== ###
### GMM
### ======================================================================== ###
### 1-step GMM estimation to obtain consistent estimates
estGMM <- solnp(pars=theta0[1:2], fun=Q_GMM, LB=c(0, 0), UB=c(1, 0.5),
                ineqfun=inequal, ineqLB=0, ineqUB=1, n=n, zVec=zV)
print(paste("The estimated parameters [beta, gamma] are:", 
            paste(round(estGMM$pars, 3), collapse = " ")))
fgrad <- function(theta, n, zVec, retE) {
    rho <- rho_GMM(theta=theta, n=n, zVec=zVec, retE=FALSE)
    return(colMeans(rho))
}
dg <- numgrad(fgrad, estGMM$pars, n=n, zVec=zV, retE=FALSE)
rho <- rho_GMM(theta=estGMM$pars, n=n, zVec=zV, retE=FALSE)
Omega <- cov_NW(rho=rho, lmax=10)
WI <- diag(ncol(rho)) 
V  <- solve(t(dg) %*% WI %*% dg)
VMid <- t(dg) %*% WI %*% Omega %*% WI %*% dg
cov <- (V %*% VMid  %*% V)/n
se <- sqrt(diag(cov))
tStats <- estGMM$pars/se
print(paste("The t-statistics for the estimated parameters are:", 
            paste(round(tStats, 3), collapse = " ")))
### ======================================================================== ###
### MLE
### ======================================================================== ###
res <- solnp(pars=theta0, fun=LL_Weibull, LB=c(0, 0, 0), UB=c(1, 0.5, 10), 
             ineqfun=inequal, ineqLB=0, ineqUB=1, n=n, zVec=zV)
print(paste("The estimated parameters [beta, gamma, shape] are:", 
            paste(round(res$pars, 3), collapse = " ")))
resEff <- MLE_Weibull_Eff(theta=res$pars[1:2], shape=res$pars[3], n=n, zVec=zV)
print(paste("The estimated parameters [beta, gamma, shape] after one-step update are:", 
            paste(round(resEff$theta, 3), collapse = " ")))
print(paste("The t-statistics for the estimated parameters are:", 
            paste(round(resEff$tStats, 3), collapse = " ")))
### ======================================================================== ###



