### ======================================================================== ###
rm(list=ls())
### ======================================================================== ###
library("Rsolnp")
### ======================================================================== ###
source("main_data.r")
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


### ======================================================================== ###
### Step 2: GMM to estimate the dynamic parameters of lambda_t process
### ======================================================================== ###
theta0 <- c(0.9, 0.05, 1.3)
inequal <- function(theta, n, zVec) {
    return(theta[1]+theta[2])
}
### 1-step GMM estimation to obtain consistent estimates 
resGMM <- solnp(pars=theta0[1:2], fun=Q_GMM, LB=c(0, 0), UB=c(1, 0.5),
             ineqfun=inequal, ineqLB=0, ineqUB=1, n=n, zVec=zV)
print(paste("The estimated parameters [beta, gamma] are:", 
            paste(round(resGMM$pars, 3), collapse = " ")))
### ======================================================================== ###


### ======================================================================== ###
lambdaVec <- lambda_series(theta=resGMM$pars, n=n, zVec=zV)
zetaVec <- zV/lambdaVec
### work with improved estimator by smoothing out lt = lt/lambda_t
liquidityS2 <- liquidity/lambdaVec
sigma <- sd(zetaVec)   ### sigma_zeta
### ======================================================================== ###


### ======================================================================== ###
### Test for permanent breaks
### ======================================================================== ###
K_gaussian <- function(x) {
    1/sqrt(2*pi)*exp(-1/2*x^2)
}
K_gaussian_E <- function(x) x*K_gaussian(x)
K_gaussian_Norm <- function(x) K_gaussian(x)^2
integrate(K_gaussian, lower=-Inf, upper=Inf)
integrate(K_gaussian_E, lower=-Inf, upper=Inf)
integrate(K_gaussian_Norm, lower=-Inf, upper=Inf)
KGaussianNorm <- 0.2821
integrate(K_gaussian_Norm, lower=-Inf, upper=Inf)
### ======================================================================== ###


### ======================================================================== ###
### event: stock split 
eventSplit <- na.omit(read_csv(paste("data/", ticker, "_split.csv", sep=""),
                               col_types = cols(`Stock Splits` = col_character())))
eventSplit$Date <- as.Date(eventSplit$Date, format="%d/%m/%Y")

ith <- 1  ### ith event 
idx <- which(dates==eventSplit$Date[ith])
left <- 1:(idx-1)   ### left-hand side index
right <- (idx:n)   ### right-hand side index
h <- KernSmooth::dpill(x=tV, y=liquidityS2)   ### bandwidth
gMinusVec <- KernSmooth::locpoly(x=tV[left], y=liquidityS2[left], bandwidth=h,
                                 degree=1, gridsize=length(left))$y
gPlusVec <- KernSmooth::locpoly(x=tV[right], y=liquidityS2[right], bandwidth=h,
                                degree=1, gridsize=length(right))$y
gMinus <- tail(gMinusVec, 1)
gPlus <- head(gPlusVec, 1)
### ======================================================================== ###
KNormPlus <- KGaussianNorm*2
KNormMinus <- KGaussianNorm*2
asyErrorVar <- sigma^2*(gPlus^2*KNormPlus+gMinus^2*KNormMinus)/(n*h)
asyError <- sqrt(asyErrorVar)   
tauLR <- (gPlus-gMinus)/asyError
print(paste("The test statistics for permanent breaks:", 
            paste(round(tauLR, 3), collapse = " ")))
### ======================================================================== ###



### ======================================================================== ###
### log-likelihood assuming zeta_t follows a Weibull distribution [unit mean]
### ======================================================================== ###
LL_Weibull <- function(theta, n, zVec) {
    ### zVec: l* - re-scaled illiquidity process; rVec: returns
    ### theta=[beta, gamma, (gamma_minus), shape] for symmetric (asymmetric) specification
    shape <- tail(theta, 1)
    lamVec <- lambda_series(theta, n, zVec)
    lik <- (log(shape/zVec)+shape*log(gamma(1+1/shape)*zVec/lamVec)-
                (gamma(1+1/shape)*zVec/lamVec)^shape)
    return(-sum(lik))
}
### ======================================================================== ###
weibullDensity <- function(x, k) {
    ### k=shape
    gmx <- gamma(1+1/k)*x
    denW <- k/x*gmx^k*exp(-gmx^k)
    return(denW)
}
### ======================================================================== ###

### ======================================================================== ###
### Test for temporary effects
### ======================================================================== ###
### j=1:J=5 corresponds to -2, -1, 0, 1, 2 around the stock split date
J <- 5
### ith event; index of ith event
print(c(ith, idx))
### ======================================================================== ###
zVecMP <- numeric(n)
zVecMP[left] <- liquidityS2[left]/gMinus
zVecMP[right] <- liquidityS2[right]/gPlus
theta0ML <- c(0.9, 0.05, 1.3)
inequalML <- function(theta, n, zVec) {
    return(theta[1]+theta[2])
}
### ML estimation of the dynamic parameters [zeta~Weibull]
resML <- solnp(pars=theta0ML, fun=LL_Weibull, 
               LB=c(0, 0, 0), UB=c(1, 0.5, 10), 
               ineqfun=inequalML, ineqLB=0, ineqUB=1,
               n=n, zVec=zVecMP)
print(paste("The estimated parameters [beta, gamma, shape] are:", 
            paste(round(resML$pars, 3), collapse = " ")))
shape <- resML$pars[3]
beta <- resML$pars[1]
lambdaVecMP <- lambda_series(theta=resML$pars, n=n, zVec=zVecMP)
zetaVecMP <- zVecMP/lambdaVecMP
fLnDzeta <- function(x, shape) {
    lmdWei <- 1/gamma(1+1/shape)
    return((shape-1-shape*(x/lmdWei)^shape)/x)
} 
s2zeta <- -(1+zetaVec*fLnDzeta(x=zetaVec, shape=shape))
### ==================================================================== ###
### efficient score function test statistics
tauCAR <- cumsum(sapply(1:J, function(j) {
    tj <- idx-3+j
    tt <- tj:n
    sum(s2zeta[tt]*beta^(tt-tj)/lambdaVecMP[tt])
}))
### ==================================================================== ###
### use data in between two events to estimate the two quantiles of CAR
### not include dates too close to events
Idx0 <- ifelse(ith==1, 3, eventIdx[ith-1]+5*2)
t1Idx <- idx-5
tauCARSim <- t(sapply(Idx0:(t1Idx+J), function(r) {
    cumsum(sapply(1:J, function(j) {
        rj <- r-3+j
        tt <- (rj:n)
        sum(s2zeta[tt]*beta^(tt-rj)/lambdaVecMP[tt]) 
    }))
    
}))
print(paste("The test statistics for temporary breaks:", 
            paste(round(tauCAR[5], 3), collapse = " ")))
print(paste("The 2.5% and 97.5% quantiles for temporary breaks:", 
            paste(round(quantile(tauCARSim[,J], probs=c(0.025, 0.975)), 3), 
                  collapse = " ")))
### ==================================================================== ###



