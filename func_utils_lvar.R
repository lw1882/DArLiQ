### ====================================================================== ###
### code source: https://www.ssc.wisc.edu/~bhansen/progs/urreg.R
### ====================================================================== ###
kernel  <- 1
band    <- 0
white   <- 1
### ====================================================================== ###
### long-run variance with Bartlett kernel
### ====================================================================== ###
lr_var <- function(u){
    tu <- nrow(u)
    p <- ncol(u)
    if (white == 1){
        te <- tu-1
        au <- qr.solve(as.matrix(u[1:te,]),as.matrix(u[2:tu,]))
        e <- as.matrix(u[2:tu,]) - as.matrix(u[1:te,])%*%au
    }else{
        e <- u
        te <- tu
    }
    if (band == 0){
        eb <- as.matrix(e[1:(te-1),])
        ef <- as.matrix(e[2:te,])
        ae <- as.matrix(colSums(eb*ef)/colSums(eb^2))
        ee <- ef - eb*(matrix(1,nrow(eb),1)%*%t(ae))
        se <- as.matrix(colMeans(ee^2))
        ad <- sum((se/((1-ae)^2))^2)
        a1 <- 4*sum((ae*se/(((1-ae)^3)*(1+ae)))^2)/ad
        a2 <- 4*sum((ae*se/((1-ae)^4))^2)/ad
        if (kernel == 2){                   # Quadratic Spectral #
            bandw <- 1.3221*((a2*te)^.2)
        }   
        if (kernel == 1){                   # Parzen #
            bandw <- 2.6614*((a2*te)^.2)               
            if (bandw > (te-2)) bandw <- te-2
        }
        if (kernel == 3){                   # Bartlett #   
            bandw <- 1.1447*((a1*te)^.333)
            if (bandw > (te-2)) bandw <- te-2
        }
    }else{
        bandw <- band
    }
    
    # Estimate Covariances #
    if (kernel == 2){                             # Quadratic Spectral Kernel #
        tm <- te-1
        jb <- as.matrix(seq(1,tm,1)/bandw)
        jband <- jb*1.2*pi
        kern <- ((sin(jband)/jband - cos(jband))/(jband^2))*3
    }
    if (kernel == 1){                             # Parzen kernel #
        tm <- floor(bandw)
        if (tm > 0){
            jb <- as.matrix(seq(1,tm,1)/bandw)
            kern <- (1 - (jb^2)*6 + (jb^3)*6)*(jb <= .5)
            kern <- kern + ((1-jb)^3)*(jb > .5)*2
        }
    }
    if (kernel == 3){                             # Bartlett kernel #
        tm <- floor(bandw)
        if (tm > 0) kern <- as.matrix(1 - seq(1,tm,1)/bandw)
    }
    
    lam <- matrix(0,p,p)
    for (j in 1:tm){
        kj <- kern[j]
        lam <- lam + (t(as.matrix(e[1:(te-j),]))%*%as.matrix(e[(1+j):te,]))*kj
    }
    omega <- (t(e)%*%e + lam + t(lam))/te
    
    if (white == 1){
        eau <- solve(diag(p) - au)
        omega <- t(eau)%*%omega%*%eau
    }
    list(omega=omega,bandw=bandw)
}

