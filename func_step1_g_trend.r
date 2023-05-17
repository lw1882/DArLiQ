### ======================================================================== ###
### Step 1: estimation of g function via kernel smoothing method
###         using derived rule-of-thumb to select bandwidth
### ======================================================================== ###
g_trend_LL_bwRoT <- function(tVec, lVec, ifTheta=FALSE, ifUpdate=FALSE, 
                             sigmaZeta=NA) {
    n <- length(tVec)
    h <- bwRoT_DArLiq(tVec=tVec, lVec=lVec, ifUpdate=ifUpdate, 
                      sigmaZeta=sigmaZeta)    
    if (ifUpdate) {
        ll <- locpoly(x=tVec, y=lVec, bandwidth=h, degree=1,
                      gridsize=length(tVec))
        llTheta <- locpoly(x=tVec, y=lVec, bandwidth=h/2, degree=1,
                           gridsize=length(tVec))
        return(list(ll=ll, llTheta=llTheta))
    } 
    if (ifTheta) {
        ### undersmoothing if estimating lambda parameters e.g. h/2
        h <- h*0.5
    }
    print(round(h, 3))
    ll <- locpoly(x=tVec, y=lVec, bandwidth=h, degree=1, 
                  gridsize=length(tVec))
    return(ll)
}
### ======================================================================== ###
