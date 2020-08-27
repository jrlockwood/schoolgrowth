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


Radj <- function(M, nzlocs, adj_method, eig.tol, eig.min=NULL){
    ## function for R adjustment to PSD, maintaining fixed zeros
    ## if specified via allowing only the non-zero elements in "nzlocs"
    ##
    ## if adj_method=="nearPD", calls nearPD2 with eig.tol
    ##
    ## if adj_method=="spectral", iteratively adjusts eigenvalues that are
    ## less than (eig.tol * max) and imposing fixed zeros until
    ## all eigenvalues exceed (eig.min * max)
    ##
    ## NOTE: no sum-to-zero constraints in this case
    ## NOTE: when adj_method=="nearPD" there should be no need for iteration,
    ## although for simplificity of code, it is placed within the while() loop
    stopifnot( (class(M) == "dsCMatrix") && is.matrix(nzlocs) && is.character(adj_method) && (eig.tol >= 0.0) )
    .B <- nrow(M)

    if(!(adj_method %in% c("spectral","nearPD"))){
        stop("in Radj: invalid adj_method")
    }

    if(adj_method=="nearPD" && !is.null(eig.min)){
        stop("in Radj: eig.min applies only with adj_method=='spectral'")
    }

    if(adj_method=="spectral"){
        if(is.null(eig.min)){
            eig.min <- 0.0
        }
        if(eig.tol < eig.min){
            stop("in Radj: eig.tol must be at least as large as eig.min")
        }
    }

    Madj <- M    
    e    <- eigen(M)
    lam  <- e$values
    U    <- e$vectors
    
    if(any(lam < (eig.tol * max(lam)))){
        done <- FALSE
        while(!done){
            if(adj_method=="spectral"){
                lam[which(lam < (eig.tol * max(lam)))] <- eig.min * max(lam)
                tmp <- U %*% diag(lam) %*% t(U)
            } else {
                .res <- nearPD2(Madj, fix0s=TRUE, do2eigen=FALSE, eig.tol=eig.tol, conv.tol=1e-11, maxit=5000)
                if(!.res$converged){
                    stop("in Radj: nearPD2 did not converge; try control$Radj_method = 'spectral'")
                }
                tmp <- as.matrix(.res$mat)
            }
            
            Madj <- sparseMatrix(i = nzlocs[,1], j = nzlocs[,2], x = tmp[nzlocs], dims=c(.B,.B), symmetric=TRUE)
            e    <- eigen(Madj)
            lam  <- e$values
            U    <- e$vectors
            
            done <- all(lam >= (ifelse(adj_method=="spectral", (eig.min * max(lam)), -.Machine$double.eps)))
        }
    }
    return(Madj)
}

Gadj <- function(M, adj_method, eig.tol, eig.min=NULL){
    ## function for G adjustment to PSD, accounting for sum-to-zero constraints
    ##
    ## if adj_method=="nearPD", calls nearPD2 with eig.tol
    ##
    ## if adj_method=="spectral", sets eigenvalues less than (eig.tol * max)
    ## equal to (eig.min * max), aside from the eigenvalue corresponding
    ## to the sum-to-zero constraints
    ##
    ## NOTE: fixed zeros are not allowed, unlike for R adjustment
    stopifnot( (class(M) == "dspMatrix") && is.character(adj_method) && (eig.tol >= 0.0) )
    
    if(!(adj_method %in% c("spectral","nearPD"))){
        stop("in Gadj: invalid adj_method")
    }

    if(adj_method=="nearPD" && !is.null(eig.min)){
        stop("in Gadj: eig.min applies only with adj_method=='spectral'")
    }

    if(adj_method=="spectral"){
        if(is.null(eig.min)){
            eig.min <- 0.0
        }
        if(eig.tol < eig.min){
            stop("in Gadj: eig.tol must be at least as large as eig.min")
        }
    }

    if(max(abs(c(apply(M, 1, sum), apply(M, 2, sum)))) > 1e-10){
        stop("in Gadj: M does not appear to satisfy sum-to-zero constraints")
    }

    Madj <- M
    e    <- eigen(M)
    lam  <- e$values
    U    <- e$vectors
    zloc <- which(apply(U, 2, var) < 1e-10)
    stopifnot(length(zloc)==1L)

    if(any(lam[-zloc] < (eig.tol * max(lam)))){

        if(adj_method=="spectral"){
            tofix      <- setdiff( which(lam < (eig.tol * max(lam))), zloc )
            lam[tofix] <- eig.min * max(lam)
            tmp        <- U %*% diag(lam) %*% t(U)
        } else {
            .res <- nearPD(M, do2eigen=FALSE, eig.tol=eig.tol, conv.tol=1e-11, maxit=5000)
            if(!.res$converged){
                stop("in Gadj: nearPD2 did not converge; try control$Gadj_method = 'spectral'")
            }
            tmp <- as.matrix(.res$mat)
        }

        Madj <- new("dspMatrix", Dim=rep(nrow(tmp),2), x=tmp[lower.tri(tmp, diag=TRUE)], uplo="L")
        if(max(abs(c(apply(Madj, 1, sum), apply(Madj, 2, sum)))) > 1e-10){
            stop("in Gadj: Madj does not satisfy sum-to-zero constraints")
        }
    }
    return(Madj)
}
