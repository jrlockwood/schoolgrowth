schoolgrowth <- function(d, target = NULL, target_contrast = NULL, control = list(), vE = NULL, vX = NULL){

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

    ## ###########################
    ## parse/set control parameters
    ## ###########################
    if(!is.list(control)){
        stop("'control' must be a list")
    }

    if(is.null(control$quietly)){
        control$quietly <- TRUE
    }
    
    if(is.null(control$school_nmin)){
        control$school_nmin <- 1
    }
    
    if(is.null(control$pattern_nmin)){
        control$pattern_nmin <- 30
    }
    
    if(is.null(control$blockpair_student_nmin)){
        control$blockpair_student_nmin <- 50
    }
    
    if(is.null(control$blockpair_school_nmin)){
        control$blockpair_school_nmin <- 10
    }
    
    if(is.null(control$eig.tol)){
        control$eig.tol <- 1e-06
    }
    
    mean_supplied <- !is.null(control$mean_varname)
    vE_supplied   <- !is.null(vE)
    vX_supplied   <- !is.null(vX)
    
    if(mean_supplied){
        if(!is.character(control$mean_varname) || (length(control$mean_varname) > 1) ){
            stop("'mean_varname' must be a character vector of length 1")
        }
        if( !(control$mean_varname %in% names(d)) ){
            stop("'mean_varname' not present in data 'd'")
        }
        if( !is.numeric(d[,control$mean_varname]) ){
            stop("'mean_varname' in data 'd' must be numeric")
        }
        if( any(is.na(d[,control$mean_varname])) ){
            stop("'mean_varname' in data 'd' has missing values")
        }
        if( !(vE_supplied & vX_supplied) ){
            stop("current version can use mean_varname only when vE and vX are also supplied")
        }
        d$muhat <- d[,control$mean_varname]
    }
    
    ## restrict data to schools with at least "control$school_nmin" total observations
    cat(paste("Restricting data to schools with at least",control$school_nmin,"student(s)...\n"))
    d$n <- ave(rep(1,nrow(d)), d$school, FUN=sum)
    if(!control$quietly){
        print(c(nrec = nrow(d), nsch = length(unique(d$school))))
    }
    d   <- subset(d, n >= control$school_nmin)
    if(!control$quietly){
        print(c(nrec = nrow(d), nsch = length(unique(d$school))))
    }
    d$n <- NULL
    
    ## make sure that we don't have fewer observations than control$pattern_nmin
    if(nrow(d) < control$pattern_nmin){
        stop("control$pattern_nmin is fewer than the total number of observations")
    }
    
    ## create block-level dataset "dblock" providing unique combinations of year*grade*subject
    ## and their block labels and IDs
    ## NOTE: don't use unique() due to excessive RAM usage.
    d$block    <- factor(paste("t",d$year,"_g",d$grade,"_b",d$subject,sep=""))
    d$blockid  <- as.integer(d$block)
    dblock     <- d[which(!duplicated(d$block)),c("year","grade","subject","block","blockid")]
    dblock     <- dblock[order(dblock$blockid),]
    allblocks  <- dblock$block
    B          <- length(allblocks)
    B2         <- B*(B+1)/2
    rownames(dblock) <- 1:nrow(dblock)
    
    if(B <= 2){
        stop("Current function requires at least three blocks")
    }

    ## if applicable, checks on user-supplied vE
    if(vE_supplied){
        cat("Checking user-supplied value of vE...\n")
        
        if( !is(vE,"symmetricMatrix") || !is(vE,"sparseMatrix") ){
            stop("vE must be stored as a sparse, symmetric matrix from the Matrix library")
        }
        
        if( any(is.na(vE)) ){
            stop("vE cannot contain missing values")
        }
        
        if( !((nrow(vE) == B) && all(rownames(vE) == allblocks)) ) {
            stop("vE rows not specified properly")
        }
        
        if( !((ncol(vE) == B) && all(colnames(vE) == allblocks)) ) {
            stop("vE columns not specified properly")
        }

        if( min(eigen(vE)$values) < -sqrt(.Machine$double.eps) ){
            stop("vE appears to have negative eigenvalues")
        }
    }

    ## if applicable, checks on user-supplied vX
    if(vX_supplied){
        cat("Checking user-supplied value of vX...\n")

        if( !is(vX,"symmetricMatrix") || !is(vX,"sparseMatrix") ){
            stop("vX must be stored as a sparse, symmetric matrix from the Matrix library")
        }
        
        if( any(is.na(vX)) ){
            stop("vX cannot contain missing values")
        }
        
        if( !((nrow(vX) == B) && all(rownames(vX) == allblocks)) ) {
            stop("vX rows not specified properly")
        }
        
        if( !((ncol(vX) == B) && all(colnames(vX) == allblocks)) ) {
            stop("vX columns not specified properly")
        }

        if( min(eigen(vX)$values) < -sqrt(.Machine$double.eps) ){
            stop("vX appears to have negative eigenvalues")
        }

        if( any(abs(apply(vX, 1, sum)) > 1e-8) ){
            stop("vX does not appear to satisfy required sum-to-zero constraints")
        }
    }

    ## ########################################################
    ## parse "target"
    ## ########################################################
    
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
        w1 <- (dblock$year == .finalyear)
    } else {
        .years <- gsub(" ","", unlist(strsplit(target["years"],",")))
        if(any(is.na(.years)) || !all(.years %in% dblock$year)){
            stop("'target' specifies years that do not occur in data")
        }
        w1 <- dblock$year %in% .years
    }
    
    ## subjects:
    if(!any(names(target) == "subjects")){
        target["subjects"] <- "all"
    }
    
    if(target["subjects"] == "all"){
        w2 <- rep(TRUE, nrow(dblock))
    } else {
        .subjects <- gsub(" ","", unlist(strsplit(target["subjects"],",")))
        if(any(is.na(.subjects)) || !all(.subjects %in% dblock$subject)){
            stop("'target' specifies subjects that do not occur in data")
        }
        w2 <- (dblock$subject %in% .subjects)
    }
    
    ## grades:
    if(!any(names(target) == "grades")){
        target["grades"] <- "all"
    }
    
    if(target["grades"] == "all"){
        w3 <- rep(TRUE, nrow(dblock))
    } else {
        .grades <- gsub(" ","", unlist(strsplit(target["grades"],",")))
        if(any(is.na(.grades)) || !all(.grades %in% dblock$grade)){
            stop("'target' specifies grades that do not occur in data")
        }
        w3 <- dblock$grade %in% .grades
    }
    
    ## weights:
    if(!any(names(target) == "weights")){
        target["weights"] <- "n"
    }
    
    if(!(target["weights"] %in% c("n","equal"))){
        stop("'target' weights must be 'n' or 'equal'")
    }
    
    dblock$intarget <- w1 & w2 & w3
    if(sum(dblock$intarget) < 1){
        stop("'target' implies no contributing blocks")
    }
    target_blocks   <- dblock$block[dblock$intarget]
    target_blockids <- dblock$blockid[dblock$intarget]
    
    ## ########################################################
    ## parse "target_contrast"
    ## ########################################################    
    target_contrast_blocks <- NULL
    
    if(!is.null(target_contrast)){
        
        if(!is.character(target_contrast) || is.null(names(target_contrast)) || (length(target_contrast) != 4)){
            stop("'target_contrast' must be a named character vector of length 4")
        }
        
        if(!all(sort(names(target_contrast)) == c("grades","subjects","weights","years"))){
            stop("invalid names of 'target_contrast'; names must be 'years','subjects','grades','weights'")
        }
        
        ## years:
        .years <- gsub(" ","", unlist(strsplit(target_contrast["years"],",")))
        if(any(is.na(.years)) || !all(.years %in% dblock$year)){
            stop("'target_contrast' specifies years that do not occur in data")
        }
        w1 <- dblock$year %in% .years
        
        ## subjects:
        if(target_contrast["subjects"] == "all"){
            w2 <- rep(TRUE, nrow(dblock))
        } else {
            .subjects <- gsub(" ","", unlist(strsplit(target_contrast["subjects"],",")))
            if(any(is.na(.subjects)) || !all(.subjects %in% dblock$subject)){
                stop("'target_contrast' specifies subjects that do not occur in data")
            }
            w2 <- (dblock$subject %in% .subjects)
        }
        
        ## grades:
        if(target_contrast["grades"] == "all"){
            w3 <- rep(TRUE, nrow(dblock))
        } else {
            .grades <- gsub(" ","", unlist(strsplit(target_contrast["grades"],",")))
            if(any(is.na(.grades)) || !all(.grades %in% dblock$grade)){
                stop("'target_contrast' specifies grades that do not occur in data")
            }
            w3 <- dblock$grade %in% .grades
        }
        
        ## weights:
        if(!(target_contrast["weights"] %in% c("n","equal"))){
            stop("'target_contrast' weights must be 'n' or 'equal'")
        }
        
        dblock$intarget_contrast <- w1 & w2 & w3
        if(sum(dblock$intarget_contrast) < 1){
            stop("'target_contrast' implies no contributing blocks")
        }
        target_contrast_blocks   <- dblock$block[dblock$intarget_contrast]
        target_contrast_blockids <- dblock$blockid[dblock$intarget_contrast]
    }
    
    ## #################################################################
    ## create "dblockpairs" dataframe that tracks block pairs and associated data
    ##
    ## NOTE: don't use combn() because we also want the diagonals and we
    ## need them in the correct positions
    ##
    ## NOTE: we order the elements of dblockpairs to fill the lower triangle
    ## of a BxB symmetric matrix by column
    ## #################################################################
    dblockpairs <- data.frame(index = 1:B2, blockidi = 0L, blockidj = 0L)
    wh <- 0
    for(j in 1:B){
        for(i in j:B){
            wh <- wh + 1
            dblockpairs$blockidi[wh] <- dblock$blockid[i]
            dblockpairs$blockidj[wh] <- dblock$blockid[j]
        }
    }
    dblockpairs$blocki <- dblock$block[match(dblockpairs$blockidi, dblock$blockid)]
    dblockpairs$blockj <- dblock$block[match(dblockpairs$blockidj, dblock$blockid)]
    ## CHECK (on order):
    ## print(sparseMatrix(i=dblockpairs$blockidi, j=dblockpairs$blockidj, x=dblockpairs$index, symmetric=TRUE))
    
    ## ##################################################################
    ## create "dsch" master list that will hold all school-level information.
    ## ##################################################################
    allschools  <- sort(unique(d$school))
    dsch        <- vector(length(allschools), mode="list")
    names(dsch) <- allschools
    for(s in 1:length(dsch)){
        dsch[[s]]$school <- allschools[s]
    }
    
    ## ###################################################################
    ## put number of schools in each block pair onto dblockpairs, and set
    ## vX_est for each blockpair depending on whether there are a sufficient
    ## number of schools to estimate that covariance element.
    ## ###################################################################
    cat("Computing counts of schools in each block pair...\n")
    .tab <- table(d$school, d$blockid)
    stopifnot(all(colnames(.tab) == 1:B) && all(rownames(.tab) == allschools))
    dblockpairs$nsch   <- 0L    
    dblockpairs$vX_est <- TRUE
    
    for(wh in 1:B2){
        dblockpairs$nsch[wh]   <- sum( (.tab[,dblockpairs$blockidi[wh]] * .tab[,dblockpairs$blockidj[wh]]) >= 0.99 )
        dblockpairs$vX_est[wh] <- ifelse(vX_supplied, FALSE, dblockpairs$nsch[wh] >= control$blockpair_school_nmin)
    }
    
    ## #####################################################################
    ## compute "N" for each school which is a sparse symmetric matrix that gives
    ## the number of students in each block pair.  I.e., diagonals are the
    ## number of students in a given school contributing to the growth measure
    ## for each block, and off-diagonals are the number of students in a given
    ## school contributing to the growth measures in a pair of blocks
    ## #####################################################################
    cat("Computing counts of students in each block pair by school...\n")
    tmp <- split(d[,c("school","stuid","blockid")], d$school)
    stopifnot(all(names(tmp) == names(dsch)))
    
    for(s in 1:length(dsch)){
        .tab <- table(tmp[[s]]$stuid, tmp[[s]]$blockid)
        oblocks <- as.integer(colnames(.tab))
        stopifnot(all(oblocks %in% dblock$blockid) && all(diff(oblocks) > 0))
        dsch[[s]]$oblocks <- oblocks
        dsch[[s]]$nblock  <- length(oblocks)
        
        trips <- list()
        wh       <- 0
        for(j in 1:ncol(.tab)){
            for(i in j:ncol(.tab)){
                wh <- wh + 1
                trips[[wh]] <- c(i = oblocks[i], j = oblocks[j], n = as.integer(sum(.tab[,i]*.tab[,j])))
            }
        }
        trips <- subset(as.data.frame(do.call("rbind", trips)), n > 0)
        dsch[[s]]$N <- sparseMatrix(i=trips$i, j=trips$j, x=trips$n, dims=c(B,B), symmetric=TRUE)
    }
    stopifnot(sum(sapply(dsch, function(x){ sum(diag(x$N))})) == nrow(d))
    rm(.tab,tmp); gc()
    
    ## ##################################################################
    ## define student pattern variable, where "pattern" refers to whether
    ## a gain is observed, or not observed, in each block.  gets counts of
    ## patterns and then combine small patterns until all patterns have
    ## at least control$pattern_nmin students.
    ##
    ## Also add field to "dblockpairs" that indicates the total number
    ## of students that are common to a given pair of blocks, and then define
    ## vE_est as indicator of whether there are a sufficient number of
    ## students in that block pair to estimate the corresponding error covariance.
    ## ##################################################################
    stupat           <- table(d$stuid, d$blockid)
    stopifnot( (all(colnames(stupat) == as.character(1:B))) && (nrow(stupat) == length(unique(d$stuid))) )
    dblockpairs$nstu   <- 0L
    dblockpairs$vE_est <- TRUE
    
    for(wh in 1:B2){
        dblockpairs$nstu[wh]   <- sum( stupat[,dblockpairs$blockidi[wh]] * stupat[,dblockpairs$blockidj[wh]] )
        dblockpairs$vE_est[wh] <- dblockpairs$nstu[wh] >= control$blockpair_student_nmin
    }

    ## stop if any block has too few students
    tmp <- subset(dblockpairs, blockidi == blockidj)
    if(!all(tmp$vE_est)){
        stop("At least one block has fewer than control$blockpair_student_nmin students")
    }
    
    stupat           <- data.frame(stuid = rownames(stupat), pattern = apply(as.matrix(stupat), 1, paste, collapse=""), stringsAsFactors=FALSE)
    stupat$stuid     <- as.character(stupat$stuid)
    stupat$pcount    <- ave(rep(1,nrow(stupat)), stupat$pattern, FUN=sum)    
    cat(paste("Number of patterns before collapsing:",length(unique(stupat$pattern)),"\n"))
    
    ## apply collapsing rule:
    ## 1) collapse all patterns with counts less than control$pattern_nmin into one pattern
    ## 2) if this combined pattern has count still below threshold, keep pulling in the
    ##    least-populated patterns until it crosses
    tmp <- stupat[which(!duplicated(stupat$pattern)),c("pattern","pcount")]
    stopifnot(nrow(tmp) == length(unique(stupat$pattern)))
    tmp <- tmp[order(-tmp$pcount),]
    tmp$cpattern <- tmp$pattern
    tmp$cpattern[which(tmp$pcount < control$pattern_nmin)] <- "collapsed"
    if(any(tmp$cpattern=="collapsed")){
        .w <- which(tmp$cpattern == "collapsed")
        .n <- sum(tmp$pcount[.w])
        while(.n < control$pattern_nmin){
            nextsmallest <- .w[1]-1
            tmp$cpattern[nextsmallest] <- "collapsed"
            .n <- .n + tmp$pcount[nextsmallest]
            .w  <- c(nextsmallest,.w)
        }
        tmp$pcount[.w] <- .n
    }
    cat(paste("Number of patterns after collapsing (control$pattern_nmin=",control$pattern_nmin,"): ", length(unique(tmp$cpattern)),"\n",sep=""))
    tab_patterns <- tmp
    
    ## assign the collapsed patterns to "collapsed"
    if(any(tmp$cpattern=="collapsed")){
        .cpats <- tmp$pattern[which(tmp$cpattern=="collapsed")]
        stupat$pattern[which(stupat$pattern %in% .cpats)] <- "collapsed"
    }
    
    ## create pattern ID and merge onto d
    stupat$patternid <- as.integer(as.factor(stupat$pattern))
    stopifnot(length(unique(stupat$patternid)) == length(unique(tmp$cpattern)))
    d$patternid <- stupat$patternid[match(d$stuid, stupat$stuid)]
    rm(stupat); gc()
    
    ## #############################################################
    ## checking for stratification of schools * (block*pattern indicators)
    ## #############################################################
    cat("Checking for stratification...\n")
    d$schoolid <- as.integer(as.factor(d$school))
    d$bpid     <- as.integer(as.factor(paste(d$blockid, d$patternid)))
    
    grp1.set <- d$schoolid[1]
    grp2.set <- unique(d$bpid[d$schoolid %in% grp1.set ])
    l.old    <- c(length(grp1.set), length(grp2.set))
    
    grp1.set <- unique(d$schoolid[ d$bpid %in% grp2.set ])
    grp2.set <- unique(d$bpid[ d$schoolid %in% grp1.set ])
    
    while( max( c(length(grp1.set), length(grp2.set)) - l.old ) > 0 ){
        l.old <- c(length(grp1.set), length(grp2.set))
        grp1.set <- unique(d$schoolid[ d$bpid %in% grp2.set ])
        grp2.set <- unique(d$bpid[ d$schoolid %in% grp1.set ])
    }
    connected <- (length(unique(d$schoolid)) == length(grp1.set)) && (length(unique(d$bpid)) == length(grp2.set))
    if(!connected){
        stop("schools and (block*pattern groups) are not connected; see help file for more information")
    }
    
    ## #######################################################
    ## compute estimate of residual covariance matrix using within-school,
    ## within-block, within-pattern deviations
    ## ########################################################
    if(!vE_supplied){
        d$sbp  <- as.integer(as.factor(paste(d$school, d$blockid, d$patternid)))
        d$nsbp <- ave(rep(1,nrow(d)), d$sbp, FUN=sum)
        d$e    <- d$G - ave(d$G, d$sbp)
        
        ## residual variance-covariance estimates, and numbers of students contributing to each
        dblockpairs$vE      <- 0.0
        dblockpairs$vE_nstu <- 0L

        cat("Estimating residual variances...\n")
        .ss  <- tapply(d$e^2, d$blockid, sum)
        .n   <- tapply(rep(1,nrow(d)), d$blockid, sum)
        .rdf <- .n - tapply(d$sbp, d$blockid, function(x){ length(unique(x))})
        stopifnot(all(names(.ss) == 1:B) && all(names(.n) == 1:B) && all(names(.rdf) == 1:B) )
        wh   <- which(dblockpairs$blockidi == dblockpairs$blockidj)
        stopifnot(all(dblockpairs$blockidi[wh] == 1:B))
        dblockpairs$vE_nstu[wh] <- as.vector(.n)
        dblockpairs$vE[wh]      <- as.vector(.ss / .rdf)

        cat("Estimating residual covariances...\n")
        for(wh in subset(dblockpairs, (blockidi != blockidj) & vE_est)$index){
            bi <- dblockpairs$blockidi[wh]
            bj <- dblockpairs$blockidj[wh]
            if(!control$quietly){
                cat(paste(bi, bj, "\n"))
            }
            tmp <- subset(d, blockid %in% c(bi, bj), select = c("stuid","school","blockid","sbp","nsbp","e"))
            tmp$blockid[which(tmp$blockid==bi)] <- 0
            tmp$blockid[which(tmp$blockid==bj)] <- 1
            
            tmp <- subset(tmp, stuid %in% unique(tmp$stuid[duplicated(tmp$stuid)]))
            if(nrow(tmp) > 0){
                stopifnot(all(as.data.frame(table(tmp$stuid))$Freq == 2))
                tmp <- reshape(tmp, timevar="blockid", idvar="stuid", direction="wide")
                ## restrict to sbp with at least two observations
                tmp <- subset(tmp, (nsbp.0 >= 2) & (nsbp.1 >= 2))
                if(nrow(tmp) > 0){
                    dblockpairs$vE_nstu[wh] <- nrow(tmp)
                    
                    ## get sbp pairs and compute unbiased estimator of covariance (from notes)
                    tmp$pair  <- paste(tmp$sbp.0, tmp$sbp.1,sep="-")
                    tmp$e0e1  <- ave(tmp$e.0 * tmp$e.1, tmp$pair)
                    tmp$n01   <- ave(rep(1, nrow(tmp)), tmp$pair, FUN=sum)
                    tmp       <- tmp[!duplicated(tmp$pair), c("pair","nsbp.0","nsbp.1","e0e1","n01")]
                    tmp$n0n1  <- tmp$nsbp.0 * tmp$nsbp.1
                    tmp$chat  <- (tmp$n0n1 * tmp$e0e1) / (tmp$n0n1 - tmp$nsbp.0 - tmp$nsbp.1 + tmp$n01)
                    dblockpairs$vE[wh] <- weighted.mean(tmp$chat, tmp$n01)
                }
            }
        }
        d$sbp <- d$nsbp <- d$e <- NULL
    } else {        
        dblockpairs$vE             <- vE[lower.tri(vE, diag=TRUE)]
        dblockpairs$vE_nstu        <- 0L
        is.na(dblockpairs$vE_nstu) <- TRUE
        dblockpairs$vE_est         <- FALSE
    }
    
    ## ###################################################################
    ## compute OLS estimates of regression coefficients
    ## ###################################################################
    if(!mean_supplied){
        cat("Computing OLS estimates of regression coefficients...\n")    
        Gmodel     <- c(varG = var(d$G))
        d$schoolid <- factor(d$schoolid)
        d$bpid     <- factor(d$bpid)

        ## first fit school FE only, no blocks, for R^2 calculation letting schools
        ## have as much as possible
        d$muhat <- ave(d$G, d$school)
        .e      <- sum( (d$G - d$muhat)^2 ) / (nrow(d) - length(unique(d$school)))
        Gmodel["Rsq_sfe"] <- 1.0 - (.e / Gmodel["varG"])

        ## now fit the actual model with school FE and block/pattern means
        .mdf    <- length(unique(d$schoolid)) + length(unique(d$bpid)) - 1
        .X      <- sparse.model.matrix(~schoolid - 1 + bpid, data=d, contrasts.arg=list(schoolid="contr.treatment",bpid="contr.sum"))
        stopifnot(all(.X@x %in% c(-1,0,1)) && (ncol(.X) == .mdf))
        stopifnot(length(grep("schoolid",colnames(.X))) == length(unique(d$school)))
        stopifnot(length(grep("bpid",    colnames(.X))) == length(unique(d$bpid))-1)
        .xpx    <- crossprod(.X)
        .xpy    <- crossprod(.X, d$G)
        .bhat   <- solve(.xpx, .xpy)
        d$muhat <- as.vector(.X %*% .bhat)
        stopifnot(max(abs(tapply(d$muhat, d$school, mean) - tapply(d$G, d$school, mean))) < 1e-6)
        stopifnot(max(abs(tapply(d$muhat, d$bpid,   mean) - tapply(d$G, d$bpid,   mean))) < 1e-6)

        ## add estimated means based on block/pattern FE but not including the school fixed effects,
        ## which will be needed later for estimating
        wh         <- grep("bpid", colnames(.X))
        d$muhat_bp <- as.vector(.X[,wh] %*% .bhat[wh])
        .bhat      <- as.vector(.bhat)
        names(.bhat) <- colnames(.X)
        
        ## put schoolFE on dsch
        tmp <- d[!duplicated(d$school),c("school","schoolid")]
        tmp$schoolFE <- .bhat[paste0("schoolid",tmp$schoolid)]
        stopifnot( !any(is.na(tmp$schoolFE)) )
        for(s in 1:length(dsch)){
            dsch[[s]]$schoolFE <- tmp$schoolFE[which(tmp$school == dsch[[s]]$school)]
        }

        rm(.X,.xpx,.xpy); gc()
        Gmodel["varE"] <- sum( (d$G - d$muhat)^2 ) / (nrow(d) - .mdf)
        Gmodel["Rsq_tot"]  <- 1.0 - (Gmodel["varE"]/Gmodel["varG"])
    
        ## cat("CHECKING MATRIX CALCS")
        ## d$muhat_chk <- fitted(lm(G ~ as.factor(school) + bpid, data=d))
        ## print(max(abs(d$muhat - d$muhat_chk)))
        ## d$muhat_chk <- NULL
    } else {
        cat("Bypassing mean estimation...\n")
        Gmodel <- NULL
        .bhat  <- NULL
    }
    
    ## ############################################################
    ## calculate and check various school*block aggregates, which will be used
    ## for estimating vX and ultimately for doing the BLP calculation
    ## ############################################################
    cat("Computing and checking school*block aggregate measures...\n")
    d$Y        <- ave(d$G,              d$school, d$block)
    ## this is just a placeholder to minimize later branching
    if(!mean_supplied){
        d$Ytilde   <- ave(d$G - d$muhat_bp, d$school, d$block)
    } else {
        d$Ytilde   <- 0.0
    }
    d$nsb      <- ave(rep(1,nrow(d)),   d$school, d$block, FUN=sum)
    
    tmp        <- d[!duplicated(paste(d$school, d$block)), c("school","grade","year","subject","block","blockid","nsb","muhat","Y","Ytilde")]
    tmp        <- tmp[order(tmp$school, tmp$blockid),]
    stopifnot(nrow(tmp) == length(unique(paste(d$school, d$block))))
    tmp        <- split(tmp, tmp$school)
    stopifnot(all(sapply(tmp, function(x){ all(diff(x$blockid) > 0) })) && all(names(tmp) == names(dsch)))
    for(s in 1:length(dsch)){
        dsch[[s]]$tab <- tmp[[s]]
    }

    stopifnot(all(unlist(lapply(dsch, function(x){ diag(x$N)[x$tab$blockid] - x$tab$nsb })) == 0))
    if(!mean_supplied){
        stopifnot(max(abs(sapply(dsch, function(x){ x$schoolFE - weighted.mean(x$tab$Ytilde, w = x$tab$nsb) }))) < Gmodel["varG"]*1e-10)
    }

    ## ########################################################
    ## compute provisional variance/covariance matrix "vU" of errors in school-level
    ## measures, as well as pieces needed for the contribution of the school to
    ## the WLS estimator for the elements of vX.
    ##
    ## NOTE: do this with a loop and accumulate key results so that we don't need
    ## to store large pieces for each school
    ## ########################################################
    cat("Computing required school-level quantities (may be slow when estimating vX)...\n")

    ## create vE
    ## NOTE: this is provisional, used for moment estimation, and may not be PSD.
    if(!vE_supplied){
        tmp <- subset(dblockpairs, vE_est)
        vE  <- sparseMatrix(i=tmp$blockidi, j=tmp$blockidj, x=tmp$vE, dims=c(B,B), symmetric=TRUE)
        rownames(vE) <- colnames(vE) <- allblocks
    }

    ## create matrix to implement sum-to-zero constraints
    A  <- sparseMatrix(i=c(1:(B-1), rep(B,B-1)), j=c(1:(B-1),1:(B-1)), x=c(rep(1,B-1),rep(-1,B-1)), dims=c(B,B-1))

    ## create matrix "adj_observed_moments" to hold the sum across schools of each
    ## school's observed cross-product matrix, adjusted for measurement error. We store
    ## as a (BxB) symmetric matrix
    adj_observed_moments <- sparseMatrix(i=1,j=1,x=0, dims=c(B,B), symmetric=TRUE)

    ## create matrix to hold the design matrix for the WLS estimator of vX.
    ## there are B(B+1)/2 rows, corresponding to the column-major lower triangle
    ## of the observed moments, and there are (B-1)B/2 columns, corresponding to the
    ## column-major lower triangle of vX*, the (B-1)x(B-1) full-rank covariance matrix.
    vX_Z <- sparseMatrix(i=1,j=1,x=0,dims=c(B*(B+1)/2, (B-1)*B/2))

    ## create vectors that will be used to mash matrix elements together as needed.
    ## posd is vector of diagonal positions for (B-1)x(B-1) matrix lower triangle.
    pos1 <- unlist(lapply(1:(B-1), function(i){i:(B-1)}))
    pos2 <- unlist(lapply(1:(B-1), function(i){rep(i,B-i)}))
    posd <- lapply((B-1):1, function(i){ c(1L,rep(0L,i-1)) })
    posd <- which(unlist(posd) > 0L)
        
    ## loop over schools
    for(s in 1:length(dsch)){
        x        <- dsch[[s]]
        stopifnot(all(x$oblocks == x$tab$blockid))
        b        <- x$oblocks

        ## compute error variance/covariance matrix for the school
        Ds       <- sparseMatrix(i=b, j=b, x=1.0/x$tab$nsb, dims=c(B,B), symmetric=TRUE)
        x$Ntilde <- Ds %*% x$N %*% Ds
        x$vU     <- x$Ntilde * vE
        x$Ntilde <- as(x$Ntilde, "symmetricMatrix")
        x$vU     <- as(x$vU,     "symmetricMatrix")
        x$pis    <- sparseVector(x$tab$nsb/sum(x$tab$nsb), i = x$tab$blockid, length=B)

        ## compute school contribution to adj_observed_moments, restricting to
        ## schools with 2+ blocks
        if(!vX_supplied && (x$nblock > 1)){
            Pis <- sparseMatrix(i=1,j=1,x=0, dims=c(B,B))
            for(.b in b){
                Pis[.b,] <- x$pis
            }
            Is      <- sparseMatrix(i=b,j=b,x=rep(1,x$nblock), dims=c(B,B), symmetric=TRUE)
            me_adj  <- (Is - Pis) %*% x$vU %*% t(Is - Pis)
            .y      <- sparseMatrix(i = rep(1,x$nblock), j = b, x = x$tab$Ytilde - x$schoolFE, dims=c(1,B))
            adj_observed_moments <- adj_observed_moments + (crossprod(.y) - me_adj)
            ## NOTE: equivalent calculation:
            ##
            ## .R  <- (Is - Pis) %*% sparseVector(x$tab$Ytilde, i = b, length=B)
            ## print(max(abs(crossprod(.y) - ( (.R %*% t(.R))))))

            ## ##################################################
            ## CHECK ME adjustment (method 1)
            ## 
            ## v    <- x$tab$nsb/sum(x$tab$nsb)
            ## Pis  <- matrix(v, nrow=x$nblock, ncol=x$nblock, byrow=T)
            ## S    <- x$vU[b,b]
            ## me_adj1 <- S - (Pis %*% S) - (S %*% t(Pis)) + (Pis %*% S %*% t(Pis))
            ## trips <- list()
            ## wh       <- 0
            ## for(j in 1:ncol(me_adj1)){
            ##     for(i in j:ncol(me_adj1)){
            ##         wh <- wh + 1
            ##         trips[[wh]] <- c(i = b[i], j = b[j], s = me_adj1[i,j])
            ##     }
            ## }
            ## trips <- as.data.frame(do.call("rbind", trips))
            ## me_adj1 <- sparseMatrix(i=trips$i, j=trips$j, x=trips$s, dims=c(B,B), symmetric=TRUE)
            ## print(max(abs(me_adj - me_adj1)))
            
            ## CHECK ME adjustment (method 2)
            ##
            ## v    <- x$tab$nsb/sum(x$tab$nsb)
            ## Pis <- matrix(0.0, nrow=B, ncol=B)
            ## for(i in b){
            ##     Pis[i,b] <- v
            ## }
            ## me_adj2 <- x$vU - (Pis %*% x$vU) - (x$vU %*% t(Pis)) + (Pis %*% x$vU %*% t(Pis))
            ## print(max(abs(me_adj - me_adj2)))
            
            ## CHECK ME adjustment (method 3)
            ##
            ## wsums      <- x$vU %*% x$pis
            ## wsum_wsums <- sum(x$pis * wsums)
            ## ij         <- rbind(t(combn(b,2)),cbind(b,b))
            ## piece      <- sparseMatrix(i=ij[,1], j=ij[,2], x = wsum_wsums - wsums[ij[,1]] - wsums[ij[,2]], dims=c(B,B), symmetric=TRUE)
            ## me_adj3        <- x$vU + piece
            ## print(max(abs(me_adj - me_adj3)))
            ## ##################################################
            
            ## compute this school's contribution to the design matrix for the WLS estimator of vX.
            ## The target piece is Cs %*% vX* %*% t(Cs) which is linear in the elements
            ## of vX*.  Here we compute the coefficients, and check that they do what we want.
            ##
            ## NOTE: for this step, the use of sparse matrices ended up being notably slower than
            ## standard matrices, so we just use standard matrices

            Cs     <- (Is - Pis) %*% A
            ## Zs <- sparseMatrix(i=1,j=1,x=0,dims=c(B*(B+1)/2, (B-1)*B/2))
            Zs <- matrix(0.0, nrow=B*(B+1)/2, ncol=(B-1)*B/2)
                
            ## t1 <- proc.time()
            wh <- 0
            for(j in 1:B){
                for(i in j:B){
                    wh      <- wh + 1
                    ## SLOW
                    ## Zs[wh,] <- (as(Cs[i,pos1],"sparseVector") * as(Cs[j,pos2],"sparseVector")) + (as(Cs[i,pos2],"sparseVector") * as(Cs[j,pos1],"sparseVector"))
                    ##
                    ## proceed only if we can get a nonzero result for Zs[wh,]
                    if( (i %in% b) && (j %in% b) ){ 
                        Zs[wh,] <- (Cs[i,pos1]*Cs[j,pos2]) + (Cs[i,pos2]*Cs[j,pos1])
                    }
                }
            }
            Zs[,posd] <- Zs[,posd] / 2.0
            ## t2 <- proc.time()
            ## print(t2["elapsed"] -t1["elapsed"])
            vX_Z <- vX_Z + as(Zs,"sparseMatrix")
            
            ## CHECK: does Zs do what is intended?
            ##
            ## .V      <- matrix(rnorm(B-1,sd=2),nrow=(B-1),ncol=1)
            ## .V      <- .V %*% t(.V)
            ## .Vl     <- .V[lower.tri(.V,diag=TRUE)]
            ## .target <- (Cs %*% .V %*% t(Cs))
            ## .target <- .target[lower.tri(.target,diag=TRUE)]
            ## print(max(abs((Zs %*% .Vl) - .target)))
        }
        
        dsch[[s]] <- x
        if(!control$quietly){
            cat(paste("school:",s,"\n"))
        }
    }

    ## #########################################################
    ## compute WLS estimates of vX* elements, force to PSD, and
    ## translate to estimates of vX.
    ## #########################################################
    if(!vX_supplied){
        cat("Estimating school*block variance components...\n")
        .Y     <- as.matrix(adj_observed_moments)
        .Y     <- .Y[lower.tri(.Y, diag=TRUE)]
        .W     <- sparseMatrix(i=1:B2, j=1:B2, x=dblockpairs$nsch, dims=c(B2,B2), symmetric=TRUE)
        .xpx   <- t(vX_Z) %*% .W %*% vX_Z
        .xpy   <- t(vX_Z) %*% .W %*% .Y
        vXstar <- matrix(0.0, ncol=(B-1), nrow=(B-1))
        vXstar[lower.tri(vXstar,diag=TRUE)] <- as.vector(solve(.xpx, .xpy))
        vXstar <- vXstar + t(vXstar)
        diag(vXstar) <- 0.5 * diag(vXstar)
        
        ## create "vXraw" which is based on raw vXstar, even if not PSD.
        ## this is just to save and return
        vXraw             <- A %*% vXstar %*% t(A)
        dblockpairs$vXraw <- vXraw[lower.tri(vXraw, diag=TRUE)]
        
        ## force vXstar to be PSD and create "vX" from that, which will be used for
        ## later calculations.
        ## NOTE: we don't use nearPD2 here because there are currently no forced zeros
        e      <- eigen(vXstar)
        e$values[which(e$values < max(e$values)*control$eig.tol)] <- 0.0
        vXstar.adj <- e$vectors %*% diag(e$values) %*% t(e$vectors)
        vX     <- A %*% vXstar.adj %*% t(A)
        dblockpairs$vX <- vX[lower.tri(vX, diag=TRUE)]
        vX     <- sparseMatrix(i=dblockpairs$blockidi, j=dblockpairs$blockidj, x=dblockpairs$vX, dims=c(B,B), symmetric=TRUE)
        rownames(vX) <- colnames(vX) <- allblocks
        rm(.Y,.W,.xpx,.xpy,vXstar,vXstar.adj)        
    } else {
        dblockpairs$vX <- vX[lower.tri(vX, diag=TRUE)]
    }

    ## ###################################################
    ## now that vX estimation is done, go back and force vE to be PSD if needed.
    ## NOTE: now using nearPD2 function to maintain fixed zeros
    ## ###################################################
    if(!vE_supplied){    
        names(dblockpairs)[which(names(dblockpairs)=="vE")] <- "vEraw"
        dblockpairs$vE <- dblockpairs$vEraw

        e <- eigen(vE)$values
        if(any(e < max(e)*control$eig.tol)){
            cat("Adjusting vE to make PSD...\n")
            .vE <- vE
            tmp <- nearPD2(vE, fix0s=TRUE, do2eigen=FALSE, eig.tol=control$eig.tol, conv.tol=1e-11, maxit=1000)$mat
            vE  <- as(tmp,"sparseMatrix")
            rownames(vE) <- colnames(vE) <- allblocks
            cat(paste0("Smallest eigenvalue of adjusted vE: ",min(eigen(vE)$values),"\n"))
            cat("Summary of differences between original and adjusted vE:\n")
            print(summary(c(as.matrix(vE - .vE))))
            dblockpairs$vE <- vE[lower.tri(vE, diag=TRUE)]

            ## adjust vU elements of dsch as well, since these are used in BLP
            for(s in 1:length(dsch)){
                x$vU     <- x$Ntilde * vE
                x$vU     <- as(x$vU, "symmetricMatrix")
            }
        }
    }

    ## #############################################
    ## BLP calculation
    ## #############################################
    cat("Computing smoothed school measures...\n")
    for(s in 1:length(dsch)){
        x <- dsch[[s]]
        b <- x$oblocks

        ## proceed if the school has measures for the target blocks
        ## (and if needed, target_contrast blocks)
        validschool <- any(b %in% target_blockids)
        if(!is.null(target_contrast)){
            validschool <- validschool && any(b %in% target_contrast_blockids)
        }
        
        if(validschool){            
            ## make lambda
            n      <- x$tab$nsb
            lambda <- rep(0.0, x$nblock)
            wh     <- which(b %in% target_blockids)
            ntot   <- sum(n[wh])
            if(target["weights"]=="n"){
                lambda[wh] <- n[wh] / ntot
            } else {
                lambda[wh] <- 1/length(wh)
            }
             
            ncontrast <- 0
            if(!is.null(target_contrast)){
                whc       <- which(b %in% target_contrast_blockids)
                ncontrast <- sum(n[whc])
                if(target_contrast["weights"]=="n"){
                    lambda[whc] <- -1.0 * (n[whc] / ncontrast)
                } else {
                    lambda[whc] <- -1.0 / length(whc)
                }
            }

            ## get BLP
            tmp <- blp(x$tab$Y, lambda, x$tab$muhat, as.matrix(vX[b,b,drop=F]), as.matrix(x$vU[b,b,drop=F]), eig.tol=control$eig.tol)

            ## create matrix to store weights for direct and BLP estimators.
            ## also stored indicator "obs" which indicates whether there
            ## are observed data for the corresponding block.
            weights           <- matrix(0, ncol=3, nrow=B)
            rownames(weights) <- allblocks
            colnames(weights) <- c("obs","direct","blp")

            weights[b,"obs"]    <- 1
            weights[b,"direct"] <- lambda
            ## for BLP, we want to account for the fact that, aside from the adjustment
            ## for pattern means, the school mean is sum(.pi * x$Y).  the weights will
            ## sum to 1
            .pi <- as.vector(x$pis[b])
            .Pi <- matrix(.pi, nrow=x$nblock, ncol=x$nblock, byrow=T)
            weights[b,"blp"] <- as.vector( t(lambda) %*% ( (tmp$IminusQ %*% .Pi) + tmp$Q ) )
            
            ## package up and save
            x$est     <- data.frame(school = x$school, ntot = ntot, ncontrast = ncontrast, as.data.frame(as.list(tmp$est)), stringsAsFactors=FALSE)
            x$weights <- weights
            dsch[[s]] <- x
        }
    }

    ## ####################
    ## RETURN
    ## ####################
    return(list(control                = control,
                target                 = target,
                target_contrast        = target_contrast,
                dblock                 = dblock,
                dblockpairs            = dblockpairs,
                dsch                   = dsch,
                bhat                   = .bhat,
                Gmodel                 = Gmodel,
                tab_patterns           = tab_patterns,
                vX                     = vX,
                vE                     = vE,
                target_blocks          = target_blocks,
                target_contrast_blocks = target_contrast_blocks,
                adjusted_growth        = do.call("rbind",lapply(dsch, function(x){ x$est }))))
}
