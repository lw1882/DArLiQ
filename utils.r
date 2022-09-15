### ======================================================================== ###
### numerical gradient
### ======================================================================== ###
numgrad <- function( fun,x, ... ) {
    f0  <-  fun(x, ... )             # n by 1    
    n   <- length( f0 )
    k   <- length( x )
    fdf <- array(0, dim=c(n,k) )
    
    # Compute step size 
    eps <- .Machine$double.eps
    dx  <- sqrt( eps )*( abs( x ) + eps )
    xh  <- x + dx
    dx  <- xh - x
    ind <- dx < sqrt(eps)
    dx[ind] <- sqrt(eps)
    
    # Compute gradient
    xdx <- diag(dx) + x  
    for (i in seq(k)) {
        fdf[,i] <- fun(xdx[,i], ...)    
    }   
    G0 <- kronecker(matrix(1, 1, k), f0 )                       # n by k 
    G1 <- kronecker(matrix(1, n, 1), t(dx) )  
    G  <- ( fdf-G0 ) / G1
    return(G)  
}
### ======================================================================== ###
