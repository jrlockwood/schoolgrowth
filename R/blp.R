blp <- function(Y, lambda, mu, SigmaX, SigmaU, eig.tol, checkMSE=FALSE){
    ## X has E[X] = mu, var[X] = SigmaX
    ## Y = X + U where E[U|X] = 0, var[U] = SigmaU
    ## compute BLP of lambda'X from Y, along with MSE and PRMSE, conditional
    ## on mu being fixed and known.
    ##
    ## algebra:
    ## Q = SigmaX(SigmaX + SigmaU)^{-1}
    ## BLP = lambda'[ (I-Q)mu + QY ]
    ## MSE = tr(lambda lambda' (Q SigmaU Q' + (I-Q)SigmaX(I-Q)')) = lambda' (I-Q)SigmaX lambda
    ##
    ## NOTE: all eigenvalues of (SigmaX + SigmaU) that are less than eig.tol*max eigenvalue
    ## are zeroed, and the Moore-Penrose inverse is used to get coefficients for BLP.

    ## checks
    stopifnot(is.matrix(SigmaX) && is.matrix(SigmaU))    
    stopifnot(is.numeric(Y) && is.numeric(lambda) && is.numeric(mu) && is.numeric(SigmaX) && is.numeric(SigmaU) && is.numeric(eig.tol))
    stopifnot(all(!is.na(c(Y,lambda,mu,SigmaX,SigmaU,eig.tol))))
    p <- length(Y)
    if(p==1){
        dim(SigmaX) <- dim(SigmaU) <- c(1,1)
    }
    stopifnot( (length(lambda) == p) && (length(mu) == p) && (nrow(SigmaX) == p) && (ncol(SigmaX) == p) && (nrow(SigmaU) == p) && (ncol(SigmaU) == p) )
    stopifnot(max(abs(SigmaX - t(SigmaX))) < 1e-8)
    stopifnot(max(abs(SigmaU - t(SigmaU))) < 1e-8)
    dim(Y) <- dim(lambda) <- dim(mu) <-  c(p,1)

    ## adjustment based on eigenvalues
    SigmaXpU <- SigmaX + SigmaU
    e        <- eigen(SigmaXpU)
    tozero   <- which(e$values < max(abs(e$values))*eig.tol)
    if(length(tozero) > 0){
        print("warning in blp: eigenvalues of SigmaX + SigmaU adjusted")
        e$values[tozero] <- 0.0
        SigmaXpU <- e$vectors %*% diag(e$values) %*% t(e$vectors)
    }
    
    ## computations
    est.direct <- as.vector(t(lambda) %*% Y)
    mse.direct <- as.vector(t(lambda) %*% SigmaU %*% lambda)
    Q          <- SigmaX %*% ginv(SigmaXpU, tol = eig.tol)
    if(max(abs( (Q %*% SigmaXpU) - SigmaX)) > 1e-06){
        print("warning in blp: G-inverse solution to BLP coefficients only approximate")
    }
    IminusQ    <- diag(p) - Q
    est.blp    <- as.vector(t(lambda) %*% ( (IminusQ %*% mu) + (Q %*% Y) ))
    mse.blp    <- as.vector(t(lambda) %*% (IminusQ %*% SigmaX) %*% lambda)

    if(checkMSE){
        part       <- (Q %*% SigmaU %*% t(Q)) + (IminusQ %*% SigmaX %*% t(IminusQ))
        mse.blp2    <- sum(diag( (lambda %*% t(lambda)) %*% part))
        stopifnot(abs(mse.blp - mse.blp2) < 1e-10)
    }

    ## PRMSE of BLP relative to assigning everyone the mean lambda'mu
    mse.null   <- as.vector(t(lambda) %*% SigmaX %*% lambda)
    prmse.null <- 1 - mse.blp/mse.null

    ## PRMSE of BLP relative to direct estimator
    prmse.direct <- 1 - mse.blp/mse.direct

    return(list(est     = c(est.direct = est.direct, mse.direct = mse.direct, est.blp = est.blp, mse.blp = mse.blp, prmse.null = prmse.null, prmse.direct = prmse.direct),
                Q       = Q,
                IminusQ = IminusQ))
}

## simple case
##
## library(JRLmisc)
## Y       <- rnorm(6)
## lambda  <- c(1,0,3,2,0,0)
## mu      <- c(0,1,4,0,0,1)
## sdx     <- sqrt(c(0.2, 0.3, 0.5, 0.2, 0.2, 0.3))
## SigmaX  <- diag(sdx) %*% CS(0.3,6) %*% diag(sdx)
## sdu     <- sqrt(c(0.1, 0.1, 0.15, 0.05, 0.10, 0.15))
## SigmaU  <- diag(sdu) %*% CS(0.1,6) %*% diag(sdu)
## blp(Y, lambda, mu, SigmaX, SigmaU, checkMSE=TRUE)
##
## case where eigenvalues get adjusted
##
## e <- eigen(SigmaX)
## e$values[which(e$values < 0.20)] <- -0.10
## SigmaX <- e$vectors %*% diag(e$values) %*% t(e$vectors)
## blp(Y, lambda, mu, SigmaX, SigmaU, checkMSE=TRUE)

