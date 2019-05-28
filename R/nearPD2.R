nearPD2	<- function (x,fix0s=FALSE,corr = FALSE, keepDiag = FALSE, do2eigen = TRUE, 
    doSym = FALSE, doDykstra = TRUE, only.values = FALSE, ensureSymmetry = !isSymmetric(x), 
    eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08, maxit = 100, 
    conv.norm.type = "I", trace = FALSE) 
{
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
    if (!converged) 
        warning(gettextf("'nearPD()' did not converge in %d iterations", 
            iter), domain = NA)
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
