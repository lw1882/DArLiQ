### ======================================================================== ###
### bandwidth selection using derived rule-of-thumb as derived in DArLiq
### ======================================================================== ###
bwRoT_DArLiq <- function(tVec, lVec, ifUpdate=FALSE, sigmaZeta=NA) {
    n <- length(tVec)
    ols_data <- data.frame(u=tVec, l=log(lVec))
    ols <- lm(l ~ u, data = ols_data) 
    a0 <- as.numeric(ols$coefficients[1])
    a1 <- as.numeric(ols$coefficients[2])

    # plot(tVec, log(lVec), type="l", col="darkblue", lwd=2)
    # lines(tVec, ols$fitted.values, col="darkred", lwd=2)
    
    K_gaussian <- function(x) {
        1/sqrt(2*pi)*exp(-1/2*x^2)
    }
    K_gaussian_mu1 <- function(x) x*K_gaussian(x)
    K_gaussian_mu2 <- function(x) x^2*K_gaussian(x)
    K_gaussian_Norm <- function(x) K_gaussian(x)^2
    integrate(K_gaussian, lower=-Inf, upper=Inf)
    integrate(K_gaussian_mu1, lower=-Inf, upper=Inf)
    integrate(K_gaussian_mu2, lower=-Inf, upper=Inf)
    integrate(K_gaussian_Norm, lower=-Inf, upper=Inf)
    KNorm <- 0.2821
    KMu2 <- 1
    
    if(ifUpdate) {
        ### if in the update step (working with lt/lambda), use sigma_zeta
        sigmaSq <- sigmaZeta^2
    } else {
        ### if in the initial step (working with lt directly), use long-run variance
        vt <- matrix(liquidity/exp(as.vector(ols$fitted.values))-1, ncol=1)
        sigmaSq <- lr_var(u=vt)$omega
    }
    return((KNorm*sigmaSq/(KMu2^2*a1^4))^(1/5)*n^(-1/5))
}
