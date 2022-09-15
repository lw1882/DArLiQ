### ======================================================================== ###
### local linear estimator of the trend function g(t/T)
### ======================================================================== ###
g_trend_LL <- function(tVec, lVec) {
    h <- KernSmooth::dpill(x=tVec, y=lVec)
    ll <- KernSmooth::locpoly(x=tVec, y=lVec, bandwidth=h, degree=1,
                              gridsize=length(tVec))
    return(ll)
}
### ======================================================================== ###

