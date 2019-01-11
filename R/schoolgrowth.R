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
##
## 12/4/2018:
##   -added some variance components and R^2 for model of G
##   -added option "schoolFE" to control list; defaults to TRUE; if FALSE we do not use school FE
##   -required each year*grade*subject combination to have at least two schools
##
## 12/5/2018:
##   -more checks on input field formats plus coercion to character
##   -return "dtab"
##
## 12/12/2018:
##   -added option "target_contrast" that gets negative weights, to allow estimation of contrasts
##
## 01/04/2019:
##   -some fixes so that function will still work, and produce BLP, when the number of cells is
##    only 1 or 2.  also changed syntax of rm() after computing variance component estimates
##
## 01/08/2019:
##   -some changes to residual calculations to try to reduce RAM footprint, including using
##    Matrix:::sparse.model.matrix to get the fitted values from the two levels of fixed effects
##
## 01/11/2019:
##   -updated blp() to return weights, and modified schoolgrowth() to track weights and return them
##   -changed NAMESPACE and DESCRIPTION so that Matrix library is not loaded







## TODO:
## -more validity checks on control parameters and/or input data
##
## -option to parameterize SigmaX and estimate those parameters
##
## -maybe allow target to be a list of named vectors so that we get BLPs for a vector of outcomes simultaneously
##
## -fix global binding warnings in R CMD check

schoolgrowth <- function(d, target = NULL, target_contrast = NULL, control = list(), quietly=TRUE, SigmaX = NULL, SigmaU = NULL, N = NULL, ...){

    ## basic argument checks
    reqnames <- c("stuid","school","grade","year","subject","G")

    if(!is.data.frame(d)){
        stop("data 'd' must be a data frame")
    }
    
    if(!all(reqnames %in% names(d))){
        stop("data 'd' does not contain all required variables")
    }

    if(!is.numeric(d$G)){
        stop("'G' variable in 'd' must be numeric")
    }
    
    ## get "final" year:
    if(any(is.na(as.numeric(d$year)))){
        stop("variable 'year' cannot be coerced to numeric without introducing missing values")
    }
    .finalyear <- as.character(d$year[which.max(as.numeric(d$year))])

    ## coercion to character, exiting if this creates missing values, or if there were
    ## missing values to begin with
    for(v in c("stuid","school","grade","year","subject")){
        if(!is.character(d[,v])){
            d[,v] <- as.character(d[,v])
            cat(paste("NOTE: variable '",v,"' in 'd' coerced to character.\n",sep=""))
        }
    }
    
    if( !(nrow(na.omit(d[,reqnames])) == nrow(d)) ){
        stop("data 'd' has missing values in required variables")
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
        control$school_nmin <- 1
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

    if(is.null(control$schoolFE)){
        control$schoolFE <- TRUE
    }
    
    ## restrict data to schools with at least "control$school_nmin" total observations
    cat(paste("Restricting data to schools with at least",control$school_nmin,"student(s)...\n"))
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

    ## make sure there are at least two schools per cell
    if(min(tapply(d$school, d$cell, function(x){ length(unique(x)) })) < 1.5){
        stop("Each year*grade*subject combination must have at least two schools")
    }

    ## if there is only one cell and we use schoolFE, then the machinery breaks down
    ## because the signal variance for the cells is 0 by definition.  we could continue
    ## the code and try to make exceptions but I think it makes more sense to stop and
    ## throw an error here
    if( (K==1) && control$schoolFE){
        stop("Function cannot use school fixed effects with only a single cell in the data; re-run with control$schoolFE=FALSE")
    }
    
    #################
    ## parse "target"
    #################    
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
        w1 <- (dtab$year == .finalyear)
    } else {
        .years <- gsub(" ","", unlist(strsplit(target["years"],",")))
        if(any(is.na(.years)) || !all(.years %in% dtab$year)){
            stop("'target' specifies years that do not occur in data")
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
        .subjects <- gsub(" ","", unlist(strsplit(target["subjects"],",")))
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
        .grades <- gsub(" ","", unlist(strsplit(target["grades"],",")))
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

    #################
    ## parse "target_contrast"
    #################
    target_contrast_cells <- NULL
    
    if(!is.null(target_contrast)){
        
        if(!is.character(target_contrast) || is.null(names(target_contrast)) || (length(target_contrast) != 4)){
            stop("'target_contrast' must be a named character vector of length 4")
        }
        
        if(!all(sort(names(target_contrast)) == c("grades","subjects","weights","years"))){
            stop("invalid names of 'target_contrast'; names must be 'years','subjects','grades','weights'")
        }
    
        ## years:
        .years <- gsub(" ","", unlist(strsplit(target_contrast["years"],",")))
        if(any(is.na(.years)) || !all(.years %in% dtab$year)){
            stop("'target_contrast' specifies years that do not occur in data")
        }
        w1 <- dtab$year %in% .years
    
        ## subjects:
        if(target_contrast["subjects"] == "all"){
            w2 <- rep(TRUE, nrow(dtab))
        } else {
            .subjects <- gsub(" ","", unlist(strsplit(target_contrast["subjects"],",")))
            if(any(is.na(.subjects)) || !all(.subjects %in% dtab$subject)){
                stop("'target_contrast' specifies subjects that do not occur in data")
            }
            w2 <- (dtab$subject %in% .subjects)
        }
    
        ## grades:
        if(target_contrast["grades"] == "all"){
            w3 <- rep(TRUE, nrow(dtab))
        } else {
            .grades <- gsub(" ","", unlist(strsplit(target_contrast["grades"],",")))
            if(any(is.na(.grades)) || !all(.grades %in% dtab$grade)){
                stop("'target_contrast' specifies grades that do not occur in data")
            }
            w3 <- dtab$grade %in% .grades
        }

        ## weights:
        if(!(target_contrast["weights"] %in% c("n","equal"))){
            stop("'target_contrast' weights must be 'n' or 'equal'")
        }

        dtab$intarget_contrast <- w1 & w2 & w3
        if(sum(dtab$intarget_contrast) < 1){
            stop("'target_contrast' implies no contributing cells")
        }
        target_contrast_cells <- dtab$cell[dtab$intarget_contrast]
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
    ## shrinkage will be based on residuals controlling for a block of fixed effects
    ## for grade*year*subject, and also school fixed effect if control$schoolFE==TRUE
    ##
    ## NOTE: we ignore the resulting loss of degrees of freedom in later
    ## calculations.
    ## we get residuals using two-stage regression absorbing schools if control$schoolFE==TRUE
    ## ###################################################################
    d$Y     <- ave(d$G, d$school, d$cell)
    Gmodel  <- c(varG = var(d$G))    
    
    cat("Computing residuals...\n")

    if(K > 1){ ## there is more than one cell:
        if(control$schoolFE){
            ## drop one column of design matrix to account for fact that fixed effects sum
            ## to intercept, and "cell" is a partition that also generates an intercept
            ##
            ## X   <- model.matrix(~cell - 1, data=d, contrasts.arg=list(cell = contr.treatment))[,-1,drop=FALSE]
            ## for(j in 1:ncol(X)){
            ##    X[,j] <- X[,j] - ave(X[,j], d$school)
            ## }
            ## "bhat" is cell means, "muhat" is the fixed effects part
            ## of the model, "R" are the residuals after accounting for fixed effects,
            ## "Y" are the aggregate performance measures (computed above)
            ## d$bhat  <- as.vector((model.matrix(~cell - 1, data=d, contrasts.arg=list(cell = contr.treatment))[,-1,drop=FALSE]) %*% (coef(lm( I(d$G - ave(d$G, d$school)) ~ X - 1))))
            ## d$muhat <- (ave(d$G, d$school) - ave(d$bhat, d$school)) + d$bhat
            ## d$bhat <- NULL
            .mdf    <- length(unique(d$school)) + length(unique(d$cell)) - 1
            .X      <- sparse.model.matrix(~school + cell -1, data=d)
            .xpx    <- crossprod(.X)
            .xpy    <- crossprod(.X, d$G)
            .bhat   <- solve(.xpx, .xpy)
            d$muhat <- as.vector(.X %*% .bhat)
            rm(.X,.xpx,.xpy,.bhat); gc()
            d$R     <- d$G - d$muhat
            if(max(abs(c(tapply(d$R, d$school, mean), tapply(d$R, d$cell, mean)))) > 1e-8){
                stop("Error in computing residuals; not orthogonal to design matrix")
            }
            Gmodel["varR"] <- sum(d$R^2) / (nrow(d) - .mdf)
            ## check (only with small dataset)
            ##
            ## cat("CHECKING MATRIX CALCS")
            ## d$muhat_chk <- fitted(lm(G ~ as.factor(school) + cell, data=d))
            ## print(max(abs(d$muhat - d$muhat_chk)))
            ## d$muhat_chk <- NULL
        } else {
            d$muhat <- ave(d$G, d$cell)
            d$R     <- d$G - d$muhat
            Gmodel["varR"] <- sum(d$R^2) / (nrow(d) - (length(unique(d$cell))))
        }
    } else{ ## there is only one cell, in which case control$schoolFE must be FALSE.
        d$muhat <- mean(d$G)
        d$R     <- d$G - d$muhat
        Gmodel["varR"] <- sum(d$R^2) / (nrow(d) - 1)
    }
    
    Gmodel["Rsq"]  <- 1.0 - (Gmodel["varR"]/Gmodel["varG"])
    
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
                    ## NOTE: if there is only one school, set covariance to 0
                    if(nrow(mi) > 1){
                        dcell$vX[wh] <- cov(mi[,1], mj[,1])
                    } else {
                        dcell$vX[wh] <- 0.0
                    }
                    
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

        rm(list = intersect(ls(), c("tmp","tmpi","tmpj","rvals","rvals2","rvalsi","rvalsj","ci","cj","mi","mj","sboth","counts","i","j")))
        gc()
        
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
            rm(m1); gc()
        }
        
        if(any(eigen(vU)$values < control$min_eigenvalue)){
            m1 <-  nearPD(vU, keepDiag=control$nearPD.keepDiag, maxit = 2000, ...)
            if(!m1$converged){
                stop("nearPD failed for vU; consider changing control$nearPD.keepDiag")
            }
            vU <- as.matrix(m1$mat)
            rm(m1); gc()
        }
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
    if(nrow(dsch) != nrow(unique(d[,c("school","grade","year","subject","cell")]))){
        stop("problem computing dsch")
    }
    
    blpit <- function(x){
        ## get matrix of counts for this school
        Ns <- matrix(0,ncol=K,nrow=K)
        rownames(Ns) <- colnames(Ns) <- allcells
        Ns[lower.tri(Ns, diag=TRUE)] <- N[,x$school[1]]
        Ns <- Ns + t(Ns)
        diag(Ns) <- diag(Ns)/2

        ## create matrix to store weights for direct and BLP estimators.
        ## also stored indicator "obs" which indicates whether there are observed data for the
        ## corresponding cell.
        weights           <- matrix(0, ncol=3, nrow=K)
        rownames(weights) <- allcells
        colnames(weights) <- c("obs","direct","blp")

        ## proceed if the school has measures for the target cells
        ## (and if needed, target_contrast cells)
        validschool <- any(x$cell %in% target_cells)
        if(!is.null(target_contrast)){
            validschool <- validschool && any(x$cell %in% target_contrast_cells)
        }
        
        if(validschool){
            Ns  <- Ns[x$cell, x$cell, drop=F]
            stopifnot(all(diag(Ns) == x$n))
            Ntilde <- Ns / (diag(Ns) %*% t(diag(Ns)))
            
            vXs <- vX[x$cell, x$cell,drop=F]
            vUs <- vU[x$cell, x$cell,drop=F] * Ntilde
            
            lambda <- rep(0, nrow(x))
            wh     <- which(x$cell %in% target_cells)
            ntot   <- sum(x$n[wh])
            if(target["weights"]=="n"){
                lambda[wh] <- x$n[wh] / ntot
            } else {
                lambda[wh] <- 1/length(wh)
            }

            ncontrast <- 0
            if(!is.null(target_contrast)){
                whc       <- which(x$cell %in% target_contrast_cells)
                ncontrast <- sum(x$n[whc])
                if(target_contrast["weights"]=="n"){
                    lambda[whc] <- -1.0 * (x$n[whc] / ncontrast)
                } else {
                    lambda[whc] <- -1.0 / length(whc)
                }
            }
            
            tmp <- blp(x$Y, lambda, x$muhat, vXs, vUs, etol = control$min_eigenvalue)
            weights[x$cell,"obs"]    <- 1
            weights[x$cell,"direct"] <- tmp$wgt[,"direct"]
            weights[x$cell,"blp"]    <- tmp$wgt[,"blp"]
            
            return(list(est = data.frame(school = x$school[1], ntot = ntot, ncontrast = ncontrast, as.data.frame(as.list(tmp$est)), stringsAsFactors=FALSE),
                        wgt = weights))
        }
    }
    .res    <- lapply(split(dsch, dsch$school), blpit)
    b       <- do.call("rbind", lapply(.res, function(x){ x$est }))
    weights <- lapply(.res, function(x){ x$wgt })
        
    ## ####################
    ## RETURN
    ## ####################
    .r <- list(control               = control,
               target                = target,
               target_contrast       = target_contrast,
               dtab                  = dtab,
               dcell                 = dcell,
               dsch                  = dsch,
               Gmodel                = Gmodel,
               SigmaX.raw            = SigmaX.raw,
               SigmaX                = vX,
               SigmaU.raw            = SigmaU.raw,
               SigmaU                = vU,
               N                     = NULL,
               target_cells          = target_cells,
               target_contrast_cells = target_contrast_cells,
               adjusted_growth       = b,
               weights               = weights)
    if(control$return.N){
        .r$N <- N
    }
    return(.r)
}
