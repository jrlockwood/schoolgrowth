## CHANGELOG:
## 11/27/2018: first commit to package

## TODO:
## 1) more data checks and checks on arguments
## 2) option to parameterize SigmaX and estimate those parameters
## 3) add ability to pass additional arguments to nearPD, possibly using ...
## 4) maybe allow target to be a vector so that we get BLPs for a vector of outcomes simultaneously
## 5) fix global binding warnings in R CMD check

schoolgrowth <- function(d, target, control = list(), quietly=TRUE){

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
    dtab <- unique(d[,c("year","grade","subject","cell")])

    if(!is.character(target) || (length(target) > 1)){
        stop("'target' not correctly specified; see help file")
    }
    target <- unlist(strsplit(target, split=";"))
    specs  <- vector(4, mode="list")
    names(specs)    <- c("years","subjects","grades","weights")
    specs$years     <- tolower(gsub(" ","",unlist(strsplit(target[grep("years:", target)], ":"))[2]))
    specs$subjects  <- tolower(gsub(" ","",unlist(strsplit(target[grep("subjects:", target)], ":"))[2]))
    specs$grades    <- tolower(gsub(" ","",unlist(strsplit(target[grep("grades:", target)], ":"))[2]))
    specs$weights   <- tolower(gsub(" ","",unlist(strsplit(target[grep("weights:", target)], ":"))[2]))

    if(any(sapply(specs, length) == 0)){
        stop("'target' not correctly specified; check for 'years:', 'subjects:', 'grades:' and 'weights:'")
    }

    if(specs$years != "final"){
        stop("current implementation supports only 'final' target year")
    }
    w1 <- (dtab$year == max(dtab$year))

    specs$subjects <- unlist(strsplit(specs$subjects,","))
    if(!all(specs$subjects %in% dtab$subject)){
        stop("'target' specifies subjects that do not occur in data")
    }
    w2 <- (dtab$subject %in% specs$subjects)

    if(specs$grades == "all"){
        w3 <- rep(TRUE, nrow(dtab))
    } else {
        specs$grades <- as.numeric(unlist(strsplit(specs$grades,",")))
        w3 <- dtab$grade %in% specs$grades
    }

    dtab$intarget <- w1 & w2 & w3
    if(sum(dtab$intarget) < 1){
        stop("'target' implies no contributing cells")
    }
    target_cells <- dtab$cell[dtab$intarget]

    if(!(specs$weights %in% c("n","equal"))){
        stop("'target' weights must be 'n' or 'equal'")
    }
    
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
    N <- matrix(0, ncol=length(allschools), nrow=nrow(dcell))
    colnames(N) <- allschools
    rownames(N) <- dcell$index

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

    vU[lower.tri(vU, diag=TRUE)] <- dcell$vU
    vU <- vU + t(vU)
    diag(vU) <- diag(vU)/2
    
    if(any(eigen(vX)$values < control$min_eigenvalue)){
        m1 <-  nearPD(vX, keepDiag=control$nearPD.keepDiag, maxit = 2000)
        if(!m1$converged){
            stop("nearPD failed for vX; consider changing control$nearPD.keepDiag")
        }
        vX <- as.matrix(m1$mat)
    }
    
    if(any(eigen(vU)$values < control$min_eigenvalue)){
        m1 <-  nearPD(vU, keepDiag=control$nearPD.keepDiag, maxit = 2000)
        if(!m1$converged){
            stop("nearPD failed for vU; consider changing control$nearPD.keepDiag")
        }
        vU <- as.matrix(m1$mat)
    }
    rm(m1); gc()

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
            if(specs$weights=="n"){
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
    return(list(dcell = dcell, dsch  = dsch, SigmaX = vX, SigmaU = vU, target_cells = target_cells, adjusted_growth = b))
}

