### ======================================================================== ###
rm(list=ls())
### ======================================================================== ###
library("Rsolnp")
library("xts")
library("KernSmooth")
### ======================================================================== ###
source("main_data.r")
source("func_utils.R")
source("func_utils_lvar.R")
source("func_step1_g_bandwidth.r")
source("func_step1_g_trend.r")
source("func_step2_lambda_update.r")
source("func_step2_lambda_GMM.r")
source("func_step2_lambda_MLE.r")
### ======================================================================== ###



### ======================================================================== ###
### Step 1: estimate long-run trend function g(t/T) -- local linear
### ======================================================================== ###
tVec <- (1:n)/n   # tVec=(1:n)/n
# estimated trend function (Local Linear)
gLL <- g_trend_LL_bwRoT(tVec=tVec, lVec=liquidity, ifTheta=FALSE, 
                        ifUpdate=FALSE, sigmaZeta=NA)
### ======================================================================== ###
### plot illiquidity series and the trend function
liquidityTS <- as.xts(cbind(gLL$y, liquidity), 
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
gVectheta <- g_trend_LL_bwRoT(tVec=tVec, lVec=liquidity, ifTheta=TRUE, 
                              ifUpdate=FALSE, sigmaZeta=NA)
zVec <- liquidity/gVectheta$y   # l*: re-scaled illiquidity
### ======================================================================== ###
theta0 <- c(0.95, 0.03, 1.3)
inequal <- function(theta, n, zVec) {
    return(theta[1]+theta[2])
}
## ======================================================================== ###
### GMM
### ======================================================================== ###
### 1-step GMM estimation to obtain consistent estimates
estGMM <- solnp(pars=theta0[1:2], fun=Q_GMM, LB=c(0, 0), UB=c(1, 0.5),
                ineqfun=inequal, ineqLB=0, ineqUB=1, n=n, zVec=zVec)
estGMM <- solnp(pars=theta0[1:2], fun=Q_GMM, LB=c(0, 0), UB=c(1, 0.5),
                n=n, zVec=zVec)
print(paste("The estimated parameters [beta, gamma] are:", 
            paste(round(estGMM$pars, 3), collapse = " ")))
stdError <- stdError_GMM(theta=estGMM$pars, n=n, zVec=zVec, lmax=10)
print(paste("The standard errors for the estimated parameters are:", 
            paste(round(stdError, 3), collapse = " ")))
lamVec <- lambda_series(theta=estGMM$pars, n=n, zVec=zVec)
### ==================================================================== ###
zetaVec <- zVec/lamVec
### ======================================================================== ###
### update based on the estimated lambda [work with lt/lambda]
gVecthetaU <- g_trend_LL_bwRoT(tVec=tVec, lVec=liquidity/lamVec, ifTheta=TRUE, 
                               ifUpdate=TRUE, sigmaZeta=sd(zetaVec))
zVecU <- liquidity/gVecthetaU$llTheta$y   # l*: re-scaled illiquidity
estGMMU <- solnp(pars=estGMM$pars[1:2], fun=Q_GMM, LB=c(0, 0), UB=c(1, 0.5),
                 ineqfun=inequal, ineqLB=0, ineqUB=1, n=n, zVec=zVecU)
estGMMU <- solnp(pars=estGMM$pars[1:2], fun=Q_GMM, LB=c(0, 0), UB=c(1, 0.5),
                 n=n, zVec=zVecU)
print(paste("The estimated parameters [beta, gamma] are:", 
            paste(round(estGMMU$pars, 3), collapse = " ")))
stdError <- stdError_GMM(theta=estGMMU$pars, n=n, zVec=zVecU, lmax=10)
print(paste("The standard errors for the estimated parameters are:", 
            paste(round(stdError, 3), collapse = " ")))
### ======================================================================== ###


### ======================================================================== ###
### MLE
### ======================================================================== ###
res <- solnp(pars=theta0, fun=LL_Weibull, LB=c(0, 0, 0), UB=c(1, 0.5, 10), 
             ineqfun=inequal, ineqLB=0, ineqUB=1, n=n, zVec=zVec)
res <- solnp(pars=theta0, fun=LL_Weibull, LB=c(0, 0, 0), UB=c(1, 0.5, 10), 
             n=n, zVec=zVec)
print(paste("The estimated parameters [beta, gamma, shape] are:", 
            paste(round(res$pars, 3), collapse = " ")))
resEff <- MLE_Weibull_Eff(theta=res$pars[1:2], shape=res$pars[3], n=n, zVec=zVec)
print(paste("The estimated parameters [beta, gamma, shape] after one-step update are:", 
            paste(round(resEff$theta, 3), collapse = " ")))
print(paste("The standard errors for the estimated parameters are:", 
            paste(round(resEff$se, 3), collapse = " ")))
### ======================================================================== ###



