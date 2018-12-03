## CHANGELOG:
## 11/27/2018:
##   -first commit to package
##
## 11/30/2018:
##   -fixed RAM-eating unique() operation
##   -changed "target" to be a named vector
##   -added some elements to return
##   -added ... argument that gets passed to nearPD()
##   -added arguments SigmaX, SigmaU, N that can be used to bypass variance component calculation
##
## 12/3/2018:
##   -relaxed requirement that "target" be final year; prints warning if target does not include
##    final year
##   -added some more checks while parsing "target"




## TODO:
## -option to parameterize SigmaX and estimate those parameters
##
## -maybe allow target to be a list of named vectors so that we get BLPs for a vector of outcomes simultaneously
##
## -fix global binding warnings in R CMD check

schoolgrowth <- function(d, target = NULL, control = list(), quietly=TRUE, SigmaX = NULL, SigmaU = NULL, N = NULL, ...){

    ## basic argument checks
    reqnames <- c("stuid","school","grade","year","subject","G")
    if(!all(reqnames %in% names(d))){
        stop("data 'd' does not contain all required variables")
    }
    if( !(nrow(na.omit(d[,reqnames])) == nrow(d)) ){
        stop("data 'd' has missing values in required variables")
    }

    if(!is.numeric(d$grade)){
        stop("'grade' variable in 'd' must be numeric")
    }

    if(!is.numeric(d$year)){
        stop("'year' variable in 'd' must be numeric")
    }

    if(!is.character(d$subject)){
        stop("'subject' variable in 'd' must be character")
    }

    if(!is.numeric(d$G)){
        stop("'G' variable in 'd' must be numeric")
    }

    ## SigmaX/SigmaU/N - more detailed checks later if these are not NULL
    tmp <- is.null(SigmaX) + is.null(SigmaU) + is.null(N)
    if(!(tmp %in% c(0,3))){
        stop("arguments 'SigmaX','SigmaU' and 'N' must be either all NULL, or all specified")
    }
    if(tmp == 3){
        est.vc <- TRUE
    } else {
        est.vc <- FALSE
    }
    
    ## control parameters
    if(!is.list(control)){
        stop("'control' must be a list")
    }
    
    if(is.null(control$school_nmin)){
        control$school_nmin <- 10
    }

    if(is.null(control$min_eigenvalue)){
        control$min_eigenvalue <- 1e-06
    }

    if(is.null(control$SigmaU_nmin)){
        control$SigmaU_nmin <- 0
    }

    if(is.null(control$nearPD.keepDiag)){
        control$nearPD.keepDiag <- FALSE
    }

    if(is.null(control$return.N)){
        control$return.N <- FALSE
    }
    
    ## restrict data to schools with at least "school_nmin" total observations
    cat("Restricting data to schools with sufficient numbers of students...\n")
    d$n <- ave(rep(1,nrow(d)), d$school, FUN=sum)
    if(!quietly){
        print(c(nrec = nrow(d), nsch = length(unique(d$school))))
    }
    d   <- subset(d, n >= control$school_nmin)
    if(!quietly){
        print(c(nrec = nrow(d), nsch = length(unique(d$school))))
    }
    d$n <- NULL

    ## create variable indicating unique combinations of year*grade*subject
    d$cell   <- paste("t",d$year,"_g",d$grade,"_b",d$subject,sep="")
    allcells <- sort(unique(d$cell))
    K        <- length(allcells)
    K2       <- K*(K+1)/2

    ## check and parse "target"
    ##
    dtab <- expand.grid(unique(d$year), unique(d$grade), unique(d$subject), stringsAsFactors=FALSE)
    names(dtab) <- c("year","grade","subject")
    dtab$cell <- paste("t",dtab$year,"_g",dtab$grade,"_b",dtab$subject,sep="")
    dtab <- subset(dtab, cell %in% d$cell)    
    ## NOTE: the simplest syntax for the previous step is:
    ## dtab <- unique(d[,c("year","grade","subject","cell")])
    ##
    ## but this leads to major RAM usage that doesn't appear to be freed with gc().
    ## so we construct the table manually, which doesn't eat so much RAM

    ## if target is left NULL, set defaults
    if(is.null(target)){
        target <- c(years="final", subjects="all", grades="all", weights="n")
    }

    ## now parse target
    if(!is.character(target) || is.null(names(target)) ){
        stop("'target' must be a named character vector")
    }
    
    if(!all(names(target) %in% c("years","subjects","grades","weights"))){
        stop("invalid names of 'target'; names must be among 'years','subjects','grades','weights'")
    }
    
    ## years:
    if(!any(names(target) == "years")){
        target["years"] <- "final"
    }

    if(target["years"] == "final"){
        w1 <- (dtab$year == max(dtab$year))
    } else {
        .years <- as.numeric(unlist(strsplit(target["years"],",")))
        if(any(is.na(.years)) || !all(.years %in% dtab$year)){
            stop("'target' specifies years that do not occur in data")
        }
        if(!(max(dtab$year) %in% .years)){
            warning("'target' does not include final year in data - double check that this was intended")
        }
        w1 <- dtab$year %in% .years
    }
    
    ## subjects:
    if(!any(names(target) == "subjects")){
        target["subjects"] <- "all"
    }
    
    if(target["subjects"] == "all"){
        w2 <- rep(TRUE, nrow(dtab))
    } else {
        .subjects <- gsub(" ","",unlist(strsplit(target["subjects"],",")))
        if(any(is.na(.subjects)) || !all(.subjects %in% dtab$subject)){
            stop("'target' specifies subjects that do not occur in data")
        }
        w2 <- (dtab$subject %in% .subjects)
    }
    
    ## grades:
    if(!any(names(target) == "grades")){
        target["grades"] <- "all"
    }

    if(target["grades"] == "all"){
        w3 <- rep(TRUE, nrow(dtab))
    } else {
        .grades <- as.numeric(unlist(strsplit(target["grades"],",")))
        if(any(is.na(.grades)) || !all(.grades %in% dtab$grade)){
            stop("'target' specifies grades that do not occur in data")
        }
        w3 <- dtab$grade %in% .grades
    }

    ## weights:
    if(!any(names(target) == "weights")){
        target["weights"] <- "n"
    }

    if(!(target["weights"] %in% c("n","equal"))){
        stop("'target' weights must be 'n' or 'equal'")
    }
        
    dtab$intarget <- w1 & w2 & w3
    if(sum(dtab$intarget) < 1){
        stop("'target' implies no contributing cells")
    }
    target_cells <- dtab$cell[dtab$intarget]
    
    ## #################################################################
    ## create "dcell" dataframe that tracks cell pairs and associated data
    ##
    ## NOTE: don't use combn() because we also want the diagonals and we
    ## need them in the correct positions
    ##
    ## NOTE: we order the elements of dcell to fill the lower triangle
    ## of a KxK matrix by column, so that we can later use lower.tri()
    ## #################################################################
    dcell <- data.frame(index = 1:K2, celli = rep("",K2),stringsAsFactors=F)
    dcell$cellj <- dcell$celli
    wh <- 0
    for(j in 1:K){
        for(i in j:K){
            wh <- wh + 1
            dcell$celli[wh] <- allcells[i]
            dcell$cellj[wh] <- allcells[j]
        }
    }
    ## check:
    ##
    ## tmp <- matrix(0,ncol=K,nrow=K)
    ## tmp[lower.tri(tmp, diag=TRUE)] <- dcell$index
    ## head(tmp)

    ## ##################################################################
    ## create big matrix "N" that will hold counts we will need for various
    ## computations.  It has one column for each school, and one row for each row of
    ## "dcell".  for rows of dcell that correspond to variances, N provides the
    ## number of students in each school who contribute to the performance measure.
    ## for rows of dcell that correspond to covariances, N provides the number of
    ## students in each school who contribute to both performance measures for the
    ## school.
    ## ##################################################################
    allschools <- sort(unique(d$school))
    
    if(est.vc){
        N <- matrix(0, ncol=length(allschools), nrow=nrow(dcell))
        colnames(N) <- allschools
        rownames(N) <- dcell$index
    } else {
        if( !is.numeric(N) || any(is.na(N)) || any(N < 0) ){
            stop("N must be a numeric matrix with no missing values and non-negative entries")
        }
        if(ncol(N) != length(allschools)){
            stop("N has incorrect number of columns")
        }
        if(!all(sort(colnames(N)) == allschools)){
            stop("columns of N must be named by school")
        }
        if(nrow(N) != nrow(dcell)){
            stop("N has incorrect number of rows")
        }
        ## NOTE: we don't have an easy way to check that the rows of N are properly
        ## ordered so need to have caveat in help file
    }


    ## ###################################################################
    ## shrinkage will be based on residuals controlling for school
    ## fixed effects and a block of fixed effects for grade*year*subject.
    ## NOTE: we ignore the resulting loss of degrees of freedom in later
    ## calculations.
    ## we get residuals using two-stage regression absorbing schools
    ## ###################################################################
    cat("Computing residuals...\n")
    
    ## drop one column of design matrix to account for fact that fixed effects sum
    ## to intercept, and "cell" is a partition that also generates an intercept
    X   <- model.matrix(~cell - 1, data=d, contrasts.arg=list(cell = contr.treatment))[,-1]
    for(j in 1:ncol(X)){
        X[,j] <- X[,j] - ave(X[,j], d$school)
    }
    ## "bhat" is cell means, "ahat" is school FE, "muhat" is the fixed effects part
    ## of the model, "R" is the residuals after accounting for fixed effects, "Y"
    ## is the aggregate performance measures
    d$bhat  <- as.vector((model.matrix(~cell - 1, data=d, contrasts.arg=list(cell = contr.treatment))[,-1]) %*% (coef(lm( I(d$G - ave(d$G, d$school)) ~ X - 1))))
    rm(X); gc()
    d$ahat  <- ave(d$G, d$school) - ave(d$bhat, d$school)
    d$muhat <- d$ahat + d$bhat
    d$R     <- d$G - d$muhat
    d$Y     <- ave(d$G, d$school, d$cell)
    stopifnot(max(abs( (d$Y - d$muhat) - ave(d$R, d$school, d$cell))) < 1e-10)
    ## check (only with small dataset)
    ##
    ## d$Rchk <- resid(lm(G ~ as.factor(school) + cell, data=d))
    ## stopifnot(max(abs(d$R - d$Rchk)) < 1e-6)
    ## d$Rchk <- NULL

    ## #######################################
    ## estimate "noise" and "signal" variance components using method-of-moments
    ## #######################################

    if(est.vc){
        ## "signal" and "noise" variance-covariance estimates    
        dcell$vX     <- dcell$vU     <- 0
        ## numbers of schools contributing to vX, vU    
        dcell$nsch_x <- dcell$nsch_u <- 0
        ## number of students contributing to vX, vU
        dcell$nstu_x <- dcell$nstu_u <- 0

        cat("Estimating variance components...\n")
        for(wh in 1:nrow(dcell)){
            if(!quietly){
                print(wh)
            }
            ci <- dcell$celli[wh]
            cj <- dcell$cellj[wh]
            
            if(ci == cj){
                ## ###########
                ## variances
                ## ###########
                tmp    <- subset(d, cell == ci)
                dcell$nsch_x[wh] <- length(unique(tmp$school))
                dcell$nstu_x[wh] <- nrow(tmp)
                N[wh,] <- sapply(allschools, function(s){ sum(tmp$school == s) })
                rvals  <- split(tmp$R, tmp$school)
                
                ## restrict to schools with at least two observations to compute vU[wh]
                rvals2  <- rvals[sapply(rvals, length) > 1]
                dcell$nsch_u[wh] <- length(rvals2)
                dcell$nstu_u[wh] <- length(unlist(rvals2))
                if(length(rvals2) > 0){
                    dcell$vU[wh]     <- weighted.mean(sapply(rvals2, var), w = (sapply(rvals2, length) - 1))
                }
                
                ## compute vX[wh]
                dcell$vX[wh] <- var(sapply(rvals, mean)) - (dcell$vU[wh] * mean(1/sapply(rvals, length)))
            } else {
                ## ############
                ## covariances
                ##
                ## NOTE: here I required the student to be in the same school for both
                ## performance measures in order to contribute to the estimate of the
                ## noise covariance, since we really are interested in the covariance
                ## for this subpopulation of students
                ## ############
                tmpi    <- subset(d, cell == ci, select = c("stuid","school","R"))
                tmpj    <- subset(d, cell == cj, select = c("stuid","school","R"))
                
                rvalsi <- split(tmpi$R, tmpi$school)
                rvalsj <- split(tmpj$R, tmpj$school)
                sboth  <- intersect(names(rvalsi), names(rvalsj))
                dcell$nsch_x[wh] <- length(sboth)
                dcell$nstu_x[wh] <- length(unique(c(tmpi$stuid[which(tmpi$school %in% sboth)], tmpj$stuid[which(tmpj$school %in% sboth)])))
                
                if(length(sboth) > 0){
                    ## pieces we will need later
                    mi <- t(sapply(rvalsi[sboth], function(x){ c(mean(x), length(x)) }))
                    mj <- t(sapply(rvalsj[sboth], function(x){ c(mean(x), length(x)) }))
                    
                    ## compute provisional value of dcell$vX[wh], which may be adjusted later for shared students
                    dcell$vX[wh] <- cov(mi[,1], mj[,1])
                    
                    ## now merge values for same student and restrict
                    tmp <- subset(na.omit(merge(tmpi, tmpj, by="stuid")), school.x == school.y)
                    if(nrow(tmp) >= 1){
                        tmp$school   <- tmp$school.x
                        tmp$school.x <- tmp$school.y <- NULL
                        tmp$n12      <- ave(rep(1,nrow(tmp)), tmp$school, FUN=sum)
                        counts       <- unique(tmp[,c("school","n12")])
                        
                        N[wh,] <- sapply(allschools, function(s){ sum(tmp$school == s) })
                        rvals  <- split(tmp, tmp$school)
                        
                        ## restrict to schools with at least two observations to compute vU[wh]
                        rvals <- rvals[sapply(rvals, nrow) >= 2]
                        if(length(rvals) >= 1){
                            dcell$nsch_u[wh] <- length(rvals)
                            dcell$nstu_u[wh] <- sum(sapply(rvals, nrow))
                            
                            ## compute estimates of covariance components
                            dcell$vU[wh]     <- weighted.mean(sapply(rvals, function(x){ cov(x$R.x, x$R.y)}), w = (sapply(rvals, nrow) - 1))
                            if(dcell$nstu_u[wh] < control$SigmaU_nmin){
                                dcell$vU[wh] <- 0.0
                            }
                            tmp <- as.data.frame(cbind(mi, mj))
                            names(tmp) <- c("r1","n1","r2","n2")
                            tmp$school <- rownames(tmp)
                            tmp <- merge(tmp, counts, by="school", all.x=TRUE)
                            tmp$n12[which(is.na(tmp$n12))] <- 0
                            dcell$vX[wh]     <- dcell$vX[wh] - (dcell$vU[wh] * mean( tmp$n12 / (tmp$n1 * tmp$n2) ))
                        }
                    }
                }
            }
        }

        rm(tmp, tmpi, tmpj, rvals, rvals2, rvalsi, rvalsj, ci, cj, mi, mj, sboth, counts, i, j); gc()
    
        ## ######################################
        ## build, check, adjust covariance matrices.
        ## try to keep diagonals; if that fails, relax diagonals
        ## ######################################
        vX <- matrix(0, ncol=K, nrow=K)
        rownames(vX) <- colnames(vX) <- allcells
        vU <- vX
        
        vX[lower.tri(vX, diag=TRUE)] <- dcell$vX
        vX <- vX + t(vX)
        diag(vX) <- diag(vX)/2
        SigmaX.raw <- vX
        
        vU[lower.tri(vU, diag=TRUE)] <- dcell$vU
        vU <- vU + t(vU)
        diag(vU) <- diag(vU)/2
        SigmaU.raw <- vU
        
        if(any(eigen(vX)$values < control$min_eigenvalue)){
            m1 <-  nearPD(vX, keepDiag=control$nearPD.keepDiag, maxit = 2000, ...)
            if(!m1$converged){
                stop("nearPD failed for vX; consider changing control$nearPD.keepDiag")
            }
            vX <- as.matrix(m1$mat)
        }
        
        if(any(eigen(vU)$values < control$min_eigenvalue)){
            m1 <-  nearPD(vU, keepDiag=control$nearPD.keepDiag, maxit = 2000, ...)
            if(!m1$converged){
                stop("nearPD failed for vU; consider changing control$nearPD.keepDiag")
            }
            vU <- as.matrix(m1$mat)
        }
        rm(m1); gc()
    } else {
        cat("Bypassing estimating variance components...\n")

        if( !is.numeric(SigmaX) || any(is.na(SigmaX)) ){
            stop("SigmaX must be a numeric matrix with no missing values")
        }

        if( !((nrow(SigmaX) == K) && all(rownames(SigmaX) == allcells)) ) {
            stop("SigmaX rows not specified properly")
        }

        if( !((ncol(SigmaX) == K) && all(colnames(SigmaX) == allcells)) ) {
            stop("SigmaX columns not specified properly")
        }

        if( !is.numeric(SigmaU) || any(is.na(SigmaU)) ){
            stop("SigmaU must be a numeric matrix with no missing values")
        }

        if( !((nrow(SigmaU) == K) && all(rownames(SigmaU) == allcells)) ) {
            stop("SigmaU rows not specified properly")
        }

        if( !((ncol(SigmaU) == K) && all(colnames(SigmaU) == allcells)) ) {
            stop("SigmaU columns not specified properly")
        }

        SigmaX.raw <- NULL
        SigmaU.raw <- NULL
        vX         <- SigmaX
        vU         <- SigmaU
    }

    ## #################
    ## BLP calculation
    ## #################
    cat("Computing smoothed school measures...\n")
    d$n    <- ave(rep(1,nrow(d)), d$school, d$cell, FUN=sum)
    dsch   <- unique(d[,c("school","grade","year","subject","cell","n","muhat","Y")])
    
    blpit <- function(x){
        ## get matrix of counts for this school
        Ns <- matrix(0,ncol=K,nrow=K)
        rownames(Ns) <- colnames(Ns) <- allcells
        Ns[lower.tri(Ns, diag=TRUE)] <- N[,x$school[1]]
        Ns <- Ns + t(Ns)
        diag(Ns) <- diag(Ns)/2

        ## proceed if the school has measures for the target cells
        if(any(x$cell %in% target_cells)){
            Ns  <- Ns[x$cell, x$cell, drop=F]
            stopifnot(all(diag(Ns) == x$n))
            Ntilde <- Ns / (diag(Ns) %*% t(diag(Ns)))
            
            vXs <- vX[x$cell, x$cell,drop=F]
            vUs <- vU[x$cell, x$cell,drop=F] * Ntilde
            
            lambda <- rep(0, nrow(x))
            wh     <- which(x$cell %in% target_cells)
            if(target["weights"]=="n"){
                lambda[wh] <- x$n[wh] / sum(x$n[wh])
            } else {
                lambda[wh] <- 1/length(wh)
            }
            
            tmp <- blp(x$Y, lambda, x$muhat, vXs, vUs, etol = control$min_eigenvalue)
            return(data.frame(school = x$school[1], ntot = sum(x$n[wh]), as.data.frame(as.list(tmp)), stringsAsFactors=FALSE))
        }
    }
    b <- do.call("rbind",lapply( split(dsch, dsch$school), blpit))

    ## ####################
    ## RETURN
    ## ####################
    .r <- list(control         = control,
               target          = target,
               dcell           = dcell,
               dsch            = dsch,
               SigmaX.raw      = SigmaX.raw,
               SigmaX          = vX,
               SigmaU.raw      = SigmaU.raw,
               SigmaU          = vU,
               N               = NULL,
               target_cells    = target_cells,
               adjusted_growth = b)    
    if(control$return.N){
        .r$N <- N
    }
    return(.r)
}

