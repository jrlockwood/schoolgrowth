nearPD2	<- function (x,fix0s=FALSE,corr = FALSE, keepDiag = FALSE, do2eigen = TRUE, 
    doSym = FALSE, doDykstra = TRUE, only.values = FALSE, ensureSymmetry = !isSymmetric(x), 
    eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08, maxit = 100, 
    conv.norm.type = "I", trace = FALSE)
{
    ## adaptation of nearPD() from library(Matrix) to enforce fixed zeros via 'fix0s' argument
    if(fix0s && !is(x,"sparseMatrix")){
        stop("fix0s option valid only when x is a sparse matrix")
    }
    if (ensureSymmetry) {
        x <- symmpart(x)
    }
    n <- ncol(x)
    if (keepDiag) 
        diagX0 <- diag(x)
    if (doDykstra) {
        D_S <- x
        D_S[] <- 0
    }
    X <- x
    iter <- 0
    converged <- FALSE
    conv <- Inf
    while (iter < maxit && !converged) {
        Y <- X
        if (doDykstra) 
            R <- Y - D_S
        e <- eigen(if (doDykstra) 
            R
        else Y, symmetric = TRUE)
        Q <- e$vectors
        d <- e$values
        p <- d > eig.tol * d[1]
        if (!any(p)) 
            stop("Matrix seems negative semi-definite")
        Q <- Q[, p, drop = FALSE]
        X <- tcrossprod(Q * rep(d[p], each = nrow(Q)), Q)
        if (doDykstra) 
            D_S <- X - R
        if (doSym) 
            X <- (X + t(X))/2
        if (corr) 
            diag(X) <- 1
        else if (keepDiag) 
            diag(X) <- diagX0
##added
##	X[which0]	<- 0
	if(fix0s==TRUE) {
		vind		<- data.frame(i=x@i)
		vind$Row	<- rep(1:nrow(x),diff(x@p))
		vind$Col	<- vind$i+1
		vind$ind	<- (vind$Col-1)*nrow(x)+vind$Row
		X		<- sparseMatrix(i=vind$Row,j=vind$Col,x=X[vind$ind],symmetric=TRUE,dims=c(nrow(X),ncol(X)))
	}

        conv <- norm(Y - X, conv.norm.type)/norm(Y, conv.norm.type)
        iter <- iter + 1
        if (trace) 
            cat(sprintf("iter %3d : #{p}=%d, ||Y-X|| / ||Y||= %11g\n", 
                iter, sum(p), conv))
        converged <- (conv <= conv.tol)
    }
    if (!converged) {
        warning(gettextf("'nearPD()' did not converge in %d iterations", 
                         iter), domain = NA)
        return(list(converged=FALSE))
    } else {
        if (do2eigen || only.values) {
            e <- eigen(X, symmetric = TRUE)
            d <- e$values
            Eps <- posd.tol * abs(d[1])
            if (d[n] < Eps) {
                d[d < Eps] <- Eps
                if (!only.values) {
                    Q <- e$vectors
                    o.diag <- diag(X)
                    X <- Q %*% (d * t(Q))
                    D <- sqrt(pmax(Eps, o.diag)/diag(X))
                    X[] <- D * X * rep(D, each = n)
                }
            }
            if (only.values) 
                return(d)
            if (corr) 
                diag(X) <- 1
            else if (keepDiag) 
                diag(X) <- diagX0
        }
        structure(list(mat = new("dpoMatrix", x = as.vector(X), Dim = c(n, 
                                                                        n)), 
                       ##	Dimnames = .M.DN(x)), ##commented out 
                       eigenvalues = d, corr = corr, 
                       normF = norm(x - X, "F"), iterations = iter, rel.tol = conv, 
                       converged = converged), class = "nearPD")
    }
}

GRadj <- function(M, nzlocs, stz, adjmethod, eig.tol, eig.min=NULL){
    ## function to iteratively adjust too small, or negative, eigenvalues and set
    ## such eigenvalues to smallest desired value, using either spectral
    ## decomposition or nearPD2, while enforcing fixed zeros, and optionally
    ## maintaining sum-to-zero constraints.  Note that if adjmethod=="nearPD",
    ## there is not option to set minimim desired eigenvalue, we just take
    ## whatever non-negative values are returned.
    stopifnot( (class(M) == "dsCMatrix") && is.logical(stz) && is.matrix(nzlocs) && is.character(adjmethod) && (eig.tol >= 0.0) )
    .B <- nrow(M)

    if(!(adjmethod %in% c("spectral","nearPD"))){
        stop("in GRadj: invalid adjmethod")
    }

    if(adjmethod=="nearPD" && !is.null(eig.min)){
        stop("in GRadj: eig.min applies only with adjmethod=='spectral'")
    }

    if(adjmethod=="spectral"){
        if(is.null(eig.min)){
            eig.min <- 0.0
        }
        if(eig.tol < eig.min){
            stop("in GRadj: eig.tol must be at least as large as eig.min")
        }
    }
    
    if(stz && (max(abs(c(apply(M, 1, sum), apply(M, 2, sum)))) > 1e-10) ){
        stop("in GRadj: M does not appear to satisfy sum-to-zero constraints")
    }

    ## get initial spectral decomposition to decide if adjustment is needed
    e    <- eigen(M)
    lam  <- e$values
    U    <- e$vectors

    if(stz){
        zloc <- which(apply(U, 2, var) < 1e-10)
        stopifnot(length(zloc)==1L)
        adj_needed <- any(lam[-zloc] < (eig.tol * max(lam)))
    } else {
        zloc <- integer(0)
        adj_needed <- any(lam < (eig.tol * max(lam)))
    }

    Madj <- M
    if(adj_needed){
        done <- FALSE
        while(!done){
            if(adjmethod=="spectral"){
                tofix      <- setdiff( which(lam < (eig.tol * max(lam))), zloc )
                lam[tofix] <- eig.min * max(lam)
                tmp        <- U %*% diag(lam) %*% t(U)
            } else {
                .res <- nearPD2(Madj, fix0s=TRUE, do2eigen=FALSE, eig.tol=eig.tol, conv.tol=1e-11, maxit=5000, keepDiag=FALSE)
                if(!.res$converged){
                    stop("in GRadj: nearPD2 did not converge; try spectral adjustment method")
                }
                tmp <- as.matrix(.res$mat)
            }
            
            Madj <- sparseMatrix(i = nzlocs[,1], j = nzlocs[,2], x = tmp[nzlocs], dims=c(.B,.B), symmetric=TRUE)
            if(stz && (max(abs(c(apply(Madj, 1, sum), apply(Madj, 2, sum)))) > 1e-10) ){
                stop("in GRadj: Madj does not satisfy sum-to-zero constraints")
            }
            
            e    <- eigen(Madj)
            print(lam  <- e$values)
            U    <- e$vectors

            if(stz){
                zloc <- which(apply(U, 2, var) < 1e-10)
                stopifnot(length(zloc)==1L)
                done <- all(lam[-zloc] >= (eig.tol * max(lam)))
            } else {
                zloc <- integer(0)
                done <- all(lam >= (eig.tol * max(lam)))
            }
        }
    }
    return(Madj)
}
