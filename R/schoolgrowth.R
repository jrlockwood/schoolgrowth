schoolgrowth <- function(d, target = NULL, target_contrast = NULL, control = list()){

    ## basic argument checks
    reqnames <- c("stuid","school","grade","year","subject","Y")
    
    if(!is.data.frame(d)){
        stop("data 'd' must be a data frame")
    }
    
    if(!all(reqnames %in% names(d))){
        stop("data 'd' does not contain all required variables")
    }
    
    if(!is.numeric(d$Y)){
        stop("'Y' variable in 'd' must be numeric")
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
    R_supplied    <- !is.null(control$R)
    G_supplied    <- !is.null(control$G)
    
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
        if( !(R_supplied & G_supplied) ){
            stop("current version can use mean_varname only when R and G are also supplied")
        }
        d$muhat <- d[,control$mean_varname]
    }

    ## stop if there are any schools with fewer than control$school_nmin records
    if(any(as.vector(table(d$school)) < control$school_nmin)){
        stop("there are schools with fewer than control$school_nmin records")
    }
    
    ## stop if there are fewer tota records than control$pattern_nmin
    if(nrow(d) < control$pattern_nmin){
        stop("control$pattern_nmin is larger than the total number of observations")
    }
    
    ## create block-level dataset "dblock" providing unique combinations of year*grade*subject
    ## and their block labels and IDs
    ## NOTE: don't use unique() due to excessive RAM usage.
    d$block     <- factor(paste("t",d$year,"_g",d$grade,"_b",d$subject,sep=""))
    d$blockid   <- as.integer(d$block)
    d$block_n   <- ave(rep(1,nrow(d)), d$blockid, FUN=sum)
    dblock      <- d[which(!duplicated(d$block)),c("year","grade","subject","block","blockid","block_n")]
    dblock      <- dblock[order(dblock$blockid),]
    B           <- nrow(dblock)
    B2          <- B*(B+1)/2
    .blocknames <- as.character(dblock$block)
    rownames(dblock) <- 1:B
    d$block_n   <- NULL

    ## sort data by school, blockid and student, which will faciliate later
    ## bookkeeping during the second-stage mixed-model estimation
    d <- d[order(d$school, d$blockid, d$stuid),]

    ## check that data have at most one record per student per block
    if(any(duplicated(paste(d$blockid, d$stuid)))){
        stop("Data must have at most one record per student per block")
    }

    ## current code needs at least 3 blocks
    if(B <= 2){
        stop("Current function requires at least three blocks")
    }

    ## if applicable, checks on user-supplied R
    if(R_supplied){
        R <- control$R
        if(!control$quietly){
            cat("Checking user-supplied value of R...\n")
        }
        
        if( !is(R,"symmetricMatrix") || !is(R,"sparseMatrix") ){
            stop("R must be a sparse, symmetric matrix from the Matrix library")
        }
        
        if( any(is.na(R)) ){
            stop("R cannot contain missing values")
        }
        
        if( !((nrow(R) == B) && all(rownames(R) == .blocknames)) ) {
            stop("R rows not specified properly")
        }
        
        if( !((ncol(R) == B) && all(colnames(R) == .blocknames)) ) {
            stop("R columns not specified properly")
        }

        .e <- eigen(R)$values
        if(any(.e < -sqrt(.Machine$double.eps))){
            stop("R appears to have negative eigenvalues")
        }
    }

    ## if applicable, checks on user-supplied G
    if(G_supplied){
        G <- control$G
        if(!control$quietly){
            cat("Checking user-supplied value of G...\n")
        }

        if( !is(G,"symmetricMatrix") || !is(G,"sparseMatrix") ){
            stop("G must a sparse, symmetric matrix from the Matrix library")
        }
        
        if( any(is.na(G)) ){
            stop("G cannot contain missing values")
        }
        
        if( !((nrow(G) == B) && all(rownames(G) == .blocknames)) ) {
            stop("G rows not specified properly")
        }
        
        if( !((ncol(G) == B) && all(colnames(G) == .blocknames)) ) {
            stop("G columns not specified properly")
        }

        .e <- eigen(G)$values
        if(any(.e < -sqrt(.Machine$double.eps))){
            stop("G appears to have negative eigenvalues")
        }

        if( any(abs(apply(G, 1, sum)) > 1e-8) ){
            stop("G does not appear to satisfy required sum-to-zero constraints")
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
    ## of a BxB symmetric matrix using column-major order
    ## #################################################################
    dblockpairs <- data.frame(iB = 1:B2, blockidi = 0L, blockidj = 0L)
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
    ## CHECK on order:
    ## print(sparseMatrix(i=dblockpairs$blockidi, j=dblockpairs$blockidj, x=dblockpairs$iB, symmetric=TRUE))

    ## ###################################################################
    ## add another index to dblockpairs.  The existing index is "iB" which
    ## indexes elements of lower triangle of BxB symmetric matrix using
    ## column-major order.  Add "iBminus1" which will index elements
    ## of (B-1)x(B-1) upper left submatrix of symmetric BxB matrix, corresponding
    ## to parameters that are directly estimated, rather than implicitly
    ## defined by sum-to-zero constraints. This is needed to implement
    ## fixing certain elements of G* to 0 based on G_est (defined later)
    ## ###################################################################
    dblockpairs$iBminus1 <- 0L
    is.na(dblockpairs$iBminus1) <- TRUE
    wh <- 0
    for(j in 1:(B-1)){
        for(i in j:(B-1)){
            wh <- wh + 1
            loc <- which( (dblockpairs$blockidi == i) & (dblockpairs$blockidj == j) )
            dblockpairs$iBminus1[loc] <- wh
        }
    }
    dblockpairs <- dblockpairs[,c("blocki","blockj","blockidi","blockidj","iB","iBminus1")]
    stopifnot(all(is.na(subset(dblockpairs, (blockidi == B) | (blockidj == B))$iBminus1)))
    stopifnot(all(na.omit(dblockpairs$iBminus1) == 1:(B*(B-1)/2)))
    
    ## ##################################################################
    ## create "dsch" master list that will hold all school-level information.
    ## ##################################################################
    tmp  <- sort(unique(d$school))
    dsch <- vector(length(tmp), mode="list")
    names(dsch) <- tmp
    for(s in 1:length(dsch)){
        dsch[[s]]$school <- tmp[s]
    }
    
    ## ###################################################################
    ## put number of schools in each block pair onto dblockpairs, and set
    ## G_est for each blockpair depending on whether there are a sufficient
    ## number of schools to estimate that covariance element.
    ## ###################################################################
    if(!control$quietly){
        cat("Computing counts of schools in each block pair...\n")
    }
    .tab <- table(d$school, d$blockid)
    stopifnot(all(colnames(.tab) == 1:B) && all(rownames(.tab) == names(dsch)))
    dblockpairs$nsch  <- 0L
    dblockpairs$G_est <- TRUE
    
    for(wh in 1:B2){
        dblockpairs$nsch[wh]  <- sum( (.tab[,dblockpairs$blockidi[wh]] * .tab[,dblockpairs$blockidj[wh]]) >= 0.99 )
        dblockpairs$G_est[wh] <- ifelse(G_supplied, FALSE, dblockpairs$nsch[wh] >= control$blockpair_school_nmin)
    }

    ## stop if any block has too few schools
    if(!G_supplied){
        tmp <- subset(dblockpairs, blockidi == blockidj)
        if(!all(tmp$G_est)){
            stop("At least one block has fewer than control$blockpair_school_nmin schools")
        }
    }
    
    ## #####################################################################
    ## compute "N" for each school which is a sparse symmetric matrix that gives
    ## the number of students in each block pair.  I.e., diagonals are the
    ## number of students in a given school contributing to the growth measure
    ## for each block, and off-diagonals are the number of students in a given
    ## school contributing to the growth measures in a pair of blocks
    ## #####################################################################
    if(!control$quietly){
        cat("Computing counts of students in each block pair by school...\n")
    }
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
    ## R_est as indicator of whether there are a sufficient number of
    ## students in that block pair to estimate the corresponding error covariance.
    ## ##################################################################
    stupat            <- table(d$stuid, d$blockid)
    stopifnot( (all(colnames(stupat) == as.character(1:B))) && (nrow(stupat) == length(unique(d$stuid))) )    
    dblockpairs$nstu  <- 0L
    dblockpairs$R_est <- TRUE
    
    for(wh in 1:B2){
        dblockpairs$nstu[wh]  <- sum( stupat[,dblockpairs$blockidi[wh]] * stupat[,dblockpairs$blockidj[wh]] )
        dblockpairs$R_est[wh] <- dblockpairs$nstu[wh] >= control$blockpair_student_nmin
    }

    ## stop if any block has too few students
    if(!R_supplied){
        tmp <- subset(dblockpairs, blockidi == blockidj)
        if(!all(tmp$R_est)){
            stop("At least one block has fewer than control$blockpair_student_nmin students")
        }
    }
    
    stupat           <- data.frame(stuid = rownames(stupat), pattern = apply(as.matrix(stupat), 1, paste, collapse=""), stringsAsFactors=FALSE)
    stupat$stuid     <- as.character(stupat$stuid)
    stupat$pcount    <- ave(rep(1,nrow(stupat)), stupat$pattern, FUN=sum)
    if(!control$quietly){
        cat(paste("Number of patterns before collapsing:",length(unique(stupat$pattern)),"\n"))
    }
    
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
    if(!control$quietly){
        cat(paste("Number of patterns after collapsing (control$pattern_nmin=",control$pattern_nmin,"): ", length(unique(tmp$cpattern)),"\n",sep=""))
    }
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
    if(!control$quietly){
        cat("Checking for stratification...\n")
    }
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
    if(!R_supplied){
        d$sbp  <- as.integer(as.factor(paste(d$school, d$blockid, d$patternid)))
        d$nsbp <- ave(rep(1,nrow(d)), d$sbp, FUN=sum)
        d$e    <- d$Y - ave(d$Y, d$sbp)
        
        ## residual variance-covariance estimates, and numbers of students contributing to each
        dblockpairs$R      <- 0.0
        dblockpairs$R_nstu <- 0L

        if(!control$quietly){
            cat("Estimating residual variances...\n")
        }
        .ss  <- tapply(d$e^2, d$blockid, sum)
        .n   <- tapply(rep(1,nrow(d)), d$blockid, sum)
        .rdf <- .n - tapply(d$sbp, d$blockid, function(x){ length(unique(x))})
        stopifnot(all(names(.ss) == 1:B) && all(names(.n) == 1:B) && all(names(.rdf) == 1:B) )
        wh   <- which(dblockpairs$blockidi == dblockpairs$blockidj)
        stopifnot(all(dblockpairs$blockidi[wh] == 1:B))
        dblockpairs$R_nstu[wh] <- as.vector(.n)
        dblockpairs$R[wh]      <- as.vector(.ss / .rdf)

        if(!control$quietly){
            cat("Estimating residual covariances...\n")
        }
        for(wh in subset(dblockpairs, (blockidi != blockidj) & R_est)$iB){
            bi <- dblockpairs$blockidi[wh]
            bj <- dblockpairs$blockidj[wh]
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
                    dblockpairs$R_nstu[wh] <- nrow(tmp)
                    
                    ## get sbp pairs and compute unbiased estimator of covariance (from notes)
                    tmp$pair  <- paste(tmp$sbp.0, tmp$sbp.1,sep="-")
                    tmp$e0e1  <- ave(tmp$e.0 * tmp$e.1, tmp$pair)
                    tmp$n01   <- ave(rep(1, nrow(tmp)), tmp$pair, FUN=sum)
                    tmp       <- tmp[!duplicated(tmp$pair), c("pair","nsbp.0","nsbp.1","e0e1","n01")]
                    tmp$n0n1  <- tmp$nsbp.0 * tmp$nsbp.1
                    tmp$chat  <- (tmp$n0n1 * tmp$e0e1) / (tmp$n0n1 - tmp$nsbp.0 - tmp$nsbp.1 + tmp$n01)
                    dblockpairs$R[wh] <- weighted.mean(tmp$chat, tmp$n01)
                }
            }
        }
        d$sbp <- d$nsbp <- d$e <- NULL
    } else {        
        dblockpairs$R             <- R[lower.tri(R, diag=TRUE)]
        dblockpairs$R_nstu        <- 0L
        is.na(dblockpairs$R_nstu) <- TRUE
        dblockpairs$R_est         <- FALSE
    }
    
    ## ###################################################################
    ## compute OLS estimates of regression coefficients
    ## ###################################################################
    if(!mean_supplied){
        if(!control$quietly){
            cat("Computing OLS estimates of regression coefficients...\n")
        }
        modstats   <- c(ntot = nrow(d), nstu = length(unique(d$stuid)), nsch = length(unique(d$school)), varY = var(d$Y))
        d$schoolid <- factor(d$school)
        d$bpid     <- factor(d$bpid)

        ## first fit school FE only, no blocks, for R^2 calculation letting schools
        ## have as much as possible
        d$muhat <- ave(d$Y, d$school)
        .e      <- sum( (d$Y - d$muhat)^2 ) / (nrow(d) - length(unique(d$school)))
        modstats["Rsq_sfe"] <- 1.0 - (.e / modstats["varY"])

        ## now fit the actual model with school FE and block/pattern FE
        .mdf    <- length(unique(d$schoolid)) + length(unique(d$bpid)) - 1
        .X      <- sparse.model.matrix(~schoolid - 1 + bpid, data=d, contrasts.arg=list(schoolid="contr.treatment",bpid="contr.sum"))
        stopifnot(all(.X@x %in% c(-1,0,1)) && (ncol(.X) == .mdf))
        stopifnot(length(grep("schoolid",colnames(.X))) == length(unique(d$school)))
        stopifnot(length(grep("bpid",    colnames(.X))) == length(unique(d$bpid))-1)
        .xpx    <- crossprod(.X)
        .xpy    <- crossprod(.X, d$Y)
        .bhat   <- solve(.xpx, .xpy)
        d$muhat <- as.vector(.X %*% .bhat)
        stopifnot(max(abs(tapply(d$muhat, d$school, mean) - tapply(d$Y, d$school, mean))) < 1e-6)
        stopifnot(max(abs(tapply(d$muhat, d$bpid,   mean) - tapply(d$Y, d$bpid,   mean))) < 1e-6)

        ## add estimated means based on block/pattern FE but not including the school fixed effects,
        ## which are treated as fixed and known for the remainder of the estimation steps
        wh           <- grep("bpid", colnames(.X))
        d$muhat_bp   <- as.vector(.X[,wh] %*% .bhat[wh])
        .bhat        <- as.vector(.bhat)
        names(.bhat) <- colnames(.X)
        
        ## put schoolFE on dsch (provisional for now; replaced with GLS estimator later)
        tmp <- d[!duplicated(d$school),c("school","schoolid")]
        tmp$schoolFE <- .bhat[paste0("schoolid",tmp$schoolid)]
        stopifnot( !any(is.na(tmp$schoolFE)) )
        for(s in 1:length(dsch)){
            dsch[[s]]$schoolFE <- tmp$schoolFE[which(tmp$school == dsch[[s]]$school)]
        }

        rm(.X,.xpx,.xpy); gc()
        modstats["varE"]     <- sum( (d$Y - d$muhat)^2 ) / (nrow(d) - .mdf)
        modstats["Rsq_tot"]  <- 1.0 - (modstats["varE"]/modstats["varY"])
    
        ## cat("CHECKING MATRIX CALCS")
        ## d$muhat_chk <- fitted(lm(Y ~ as.factor(school) + bpid, data=d))
        ## print(max(abs(d$muhat - d$muhat_chk)))
        ## d$muhat_chk <- NULL
    } else {
        if(!control$quietly){
            cat("Bypassing mean estimation...\n")
        }
        modstats <- NULL
        .bhat    <- NULL
    }
    
    ## ############################################################
    ## calculate and check various school*block aggregates, which will be used
    ## for estimating G and ultimately for doing the BLP calculation
    ## ############################################################
    if(!control$quietly){
        cat("Computing and checking school*block aggregate measures...\n")
    }
    d$Y_sb <- ave(d$Y, d$school, d$block)

    if(!mean_supplied){
        d$Y_sb_tilde   <- ave(d$Y - d$muhat_bp, d$school, d$block)
    } else {
        ## this is just a placeholder to minimize later branching
        d$Y_sb_tilde   <- 0.0
        ## key step: if mean_supplied, it generally will vary across students within school*blocks
        ## due to patterns, and we need to average it to the school*block level to get the implied
        ## school*block mean
        d$muhat <- ave(d$muhat, d$school, d$block)
    }
    d$nsb      <- ave(rep(1,nrow(d)), d$school, d$block, FUN=sum)
    
    tmp        <- d[!duplicated(paste(d$school, d$block)), c("school","grade","year","subject","block","blockid","nsb","muhat","Y_sb","Y_sb_tilde")]
    tmp        <- tmp[order(tmp$school, tmp$blockid),]
    stopifnot(nrow(tmp) == length(unique(paste(d$school, d$block))))
    tmp        <- split(tmp, tmp$school)
    stopifnot(all(sapply(tmp, function(x){ all(diff(x$blockid) > 0) })) && all(names(tmp) == names(dsch)))
    for(s in 1:length(dsch)){
        dsch[[s]]$tab <- tmp[[s]]
    }

    stopifnot(all(unlist(lapply(dsch, function(x){ diag(x$N)[x$tab$blockid] - x$tab$nsb })) == 0))
    if(!mean_supplied){
        stopifnot(max(abs(sapply(dsch, function(x){ x$schoolFE - weighted.mean(x$tab$Y_sb_tilde, w = x$tab$nsb) }))) < modstats["varY"]*1e-10)
    }

    ## ########################################################
    ## compute provisional variance/covariance matrix "R_sb" of errors in aggregate
    ## measures, as well as pieces needed for the contribution of the school to
    ## the WLS estimator for the elements of G.
    ##
    ## NOTE: do this with a loop and accumulate key results so that we don't need
    ## to store large pieces for each school
    ## ########################################################
    if(!control$quietly){
        cat("Computing required school-level quantities (may be slow when estimating G)...\n")
    }

    ## create R
    ## NOTE: this is provisional, used for moment estimation, and may not be PSD.
    if(!R_supplied){
        tmp <- subset(dblockpairs, R_est)
        R   <- sparseMatrix(i=tmp$blockidi, j=tmp$blockidj, x=tmp$R, dims=c(B,B), symmetric=TRUE)
        rownames(R) <- colnames(R) <- .blocknames
    }

    ## create matrix to implement sum-to-zero constraints
    A  <- sparseMatrix(i=c(1:(B-1), rep(B,B-1)), j=c(1:(B-1),1:(B-1)), x=c(rep(1,B-1),rep(-1,B-1)), dims=c(B,B-1))

    ## create matrix "adj_observed_moments" to hold the sum across schools of each
    ## school's observed cross-product matrix, adjusted for measurement error. We store
    ## as a (BxB) symmetric matrix
    adj_observed_moments <- sparseMatrix(i=1,j=1,x=0, dims=c(B,B), symmetric=TRUE)

    ## create matrix to hold the design matrix for the WLS estimator of G.
    ## there are B(B+1)/2 rows, corresponding to the column-major lower triangle
    ## of the observed moments, and there are (B-1)B/2 columns, corresponding to the
    ## column-major lower triangle of G*, the (B-1)x(B-1) full-rank covariance matrix.
    ##
    ## NOTE: sparse tended to be slower here so we make it dense, which is fine
    ## because the matrix actually ends up dense.  G_Z will accumulate across
    ## schools and each school's contribution will be initialized to Zs0.
    G_Z <- Zs0  <- matrix(0.0, nrow=B*(B+1)/2, ncol=(B-1)*B/2)
    
    ## create vectors that will be used to mash matrix elements together as needed.
    ## posd is vector of diagonal positions for (B-1)x(B-1) matrix lower triangle.
    pos1 <- unlist(lapply(1:(B-1), function(i){i:(B-1)}))
    pos2 <- unlist(lapply(1:(B-1), function(i){rep(i,B-i)}))
    posd <- lapply((B-1):1, function(i){ c(1L,rep(0L,i-1)) })
    posd <- which(unlist(posd) > 0L)

    ## create a matrix to hold positions of (i,j) element of (BxB) symmetric
    ## matrix stored as column-major lower triangle.
    posB <- matrix(0,ncol=B,nrow=B)
    posB[lower.tri(posB, diag=TRUE)] <- 1:B2
    
    ## loop over schools
    ## t1 <- proc.time()
    for(s in 1:length(dsch)){
        x        <- dsch[[s]]
        stopifnot(all(x$oblocks == x$tab$blockid))
        b        <- x$oblocks

        ## compute error variance/covariance matrix for the school
        Ds       <- sparseMatrix(i=b, j=b, x=1.0/x$tab$nsb, dims=c(B,B), symmetric=TRUE)
        x$Ntilde <- Ds %*% x$N %*% Ds
        x$R_sb   <- x$Ntilde * R
        x$Ntilde <- as(x$Ntilde, "symmetricMatrix")
        x$R_sb   <- as(x$R_sb,   "symmetricMatrix")
        x$pis    <- sparseVector(x$tab$nsb/sum(x$tab$nsb), i = x$tab$blockid, length=B)

        ## compute school contribution to adj_observed_moments, restricting to
        ## schools with 2+ blocks
        if(!G_supplied && (x$nblock > 1)){
            Pis <- sparseMatrix(i=1,j=1,x=0, dims=c(B,B))
            for(.b in b){
                Pis[.b,] <- x$pis
            }
            Is      <- sparseMatrix(i=b,j=b,x=rep(1,x$nblock), dims=c(B,B), symmetric=TRUE)
            me_adj  <- (Is - Pis) %*% x$R_sb %*% t(Is - Pis)
            .y      <- sparseMatrix(i = rep(1,x$nblock), j = b, x = x$tab$Y_sb_tilde - x$schoolFE, dims=c(1,B))
            adj_observed_moments <- adj_observed_moments + (crossprod(.y) - me_adj)
            ## NOTE: equivalent calculation:
            ##
            ## .R  <- (Is - Pis) %*% sparseVector(x$tab$Y_sb_tilde, i = b, length=B)
            ## print(max(abs(crossprod(.y) - ( (.R %*% t(.R))))))

            ## ##################################################
            ## CHECK ME adjustment (method 1)
            ## 
            ## v    <- x$tab$nsb/sum(x$tab$nsb)
            ## Pis  <- matrix(v, nrow=x$nblock, ncol=x$nblock, byrow=T)
            ## S    <- x$R_sb[b,b]
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
            ## me_adj2 <- x$R_sb - (Pis %*% x$R_sb) - (x$R_sb %*% t(Pis)) + (Pis %*% x$R_sb %*% t(Pis))
            ## print(max(abs(me_adj - me_adj2)))
            
            ## CHECK ME adjustment (method 3)
            ##
            ## wsums      <- x$R_sb %*% x$pis
            ## wsum_wsums <- sum(x$pis * wsums)
            ## ij         <- rbind(t(combn(b,2)),cbind(b,b))
            ## piece      <- sparseMatrix(i=ij[,1], j=ij[,2], x = wsum_wsums - wsums[ij[,1]] - wsums[ij[,2]], dims=c(B,B), symmetric=TRUE)
            ## me_adj3        <- x$R_sb + piece
            ## print(max(abs(me_adj - me_adj3)))
            ## ##################################################
            
            ## compute this school's contribution to the design matrix for the WLS estimator of G.
            ## The target piece is Cs %*% G* %*% t(Cs) which is linear in the elements
            ## of G*.  Here we compute the coefficients, and check that they do what we want.
            ##
            ## NOTE: for this step, the use of sparse matrices ended up being notably slower than
            ## standard matrices, so we just use standard matrices
            ##
            ## NOTE: need to loop over only observed block pairs for this school because other
            ## contributions are zero
            Cs  <- as.matrix((Is - Pis) %*% A)
            Zs  <- Zs0
            
            for(j in 1:length(b)){
                for(i in j:length(b)){
                    Zs[posB[b[i],b[j]],] <- (Cs[b[i],pos1]*Cs[b[j],pos2]) + (Cs[b[i],pos2]*Cs[b[j],pos1])
                }
            }
            Zs[,posd] <- Zs[,posd] / 2.0
            G_Z <- G_Z + Zs
            
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
    }
    ## t2 <- proc.time()
    ## print(t2["elapsed"] -t1["elapsed"])
    rm(Zs0); gc()
    
    ## #########################################################
    ## compute WLS estimates of G* elements, force to PSD,
    ## and translate to estimates of G.  Note that fixed zeros
    ## for G* are introduced depending on G_est
    ## #########################################################
    if(!G_supplied){
        if(!control$quietly){
            cat("Estimating school*block variance components...\n")
        }

        ## define pieces needed to get WLS estimator, restricting to observations (rows)
        ## and parameters (columns) depending on G_est

        ## observation (row) restrictions
        rkeep  <- dblockpairs$G_est
        .Y     <- as.matrix(adj_observed_moments)
        .Y     <- .Y[lower.tri(.Y, diag=TRUE)]
        stopifnot(length(rkeep) == length(.Y))
        .Y     <- .Y[rkeep]
        G_Z    <- G_Z[rkeep,]

        tmp    <- subset(dblockpairs, G_est)
        stopifnot(all(diff(tmp$iB) > 0))
        .W     <- sparseMatrix(i=1:nrow(tmp), j=1:nrow(tmp), x=tmp$nsch, dims=c(nrow(tmp),nrow(tmp)), symmetric=TRUE)

        ## parameter (column) restrictions
        tmp    <- subset(dblockpairs, (blockidi < B) & (blockidj < B) )
        stopifnot(all(diff(tmp$iBminus1) == 1))
        ckeep  <- tmp$G_est
        G_Z    <- G_Z[,ckeep]

        ## solve for estimable parameters and construct Gstar
        tmp    <- subset(dblockpairs, (blockidi < B) & (blockidj < B) & G_est)
        stopifnot(all(diff(tmp$iBminus1) > 0))        
        .xpx   <- t(G_Z) %*% .W %*% G_Z
        .xpy   <- t(G_Z) %*% .W %*% .Y
        tmp$gstar <- as.vector(solve(.xpx, .xpy))
        Gstar  <- sparseMatrix(i=tmp$blockidi, j=tmp$blockidj, x=tmp$gstar, dims=c(B-1,B-1), symmetric=TRUE)
        
        ## create "Graw" which is based on raw Gstar, even if not PSD.
        ## this is just to save and return
        Graw             <- as.matrix(A %*% Gstar %*% t(A))
        dblockpairs$Graw <- Graw[lower.tri(Graw, diag=TRUE)]
        
        ## force Gstar to be PSD and create "G" from that, which will be used for
        ## later calculations.
        ##
        ## NOTE: sometimes nearPD2() does not converge.  if it does, we use
        ## it.  if it does not, we fall back to the spectral decomposition
        ## adjustment (which generally will not preserve fixed zeros).
        ## note also that we start with keepDiag=TRUE and revert to keepDiag=FALSE
        ## if necessary
        e <- eigen(Gstar)
        if(any(e$values < -sqrt(.Machine$double.eps))){
            if(!control$quietly){
                cat("Adjusting G* to make PSD...\n")
            }
            .Gstar <- Gstar
            tmp <- nearPD2(Gstar, fix0s=TRUE, do2eigen=FALSE, eig.tol=control$eig.tol, conv.tol=1e-11, maxit=1000, keepDiag=TRUE)
            if(tmp$converged){
                Gstar <- tmp$mat
            } else {
                tmp <- nearPD2(Gstar, fix0s=TRUE, do2eigen=FALSE, eig.tol=control$eig.tol, conv.tol=1e-11, maxit=1000, keepDiag=FALSE)
                if(tmp$converged){
                    Gstar <- tmp$mat
                } else {
                    e$values[which(e$values < max(e$values)*control$eig.tol)] <- 0.0
                    Gstar <- e$vectors %*% diag(e$values) %*% t(e$vectors)
                }
            }
            if(!control$quietly){
                cat(paste0("Smallest eigenvalue of adjusted G*: ",min(eigen(Gstar)$values),"\n"))
                cat("Summary of differences between original and adjusted G*:\n")
                print(summary(c(as.matrix(Gstar - .Gstar))))
            }
        }
        G     <- A %*% Gstar %*% t(A)
        dblockpairs$G <- G[lower.tri(G, diag=TRUE)]
        G     <- sparseMatrix(i=dblockpairs$blockidi, j=dblockpairs$blockidj, x=dblockpairs$G, dims=c(B,B), symmetric=TRUE)
        rownames(G) <- colnames(G) <- .blocknames
        rm(.Y,.W,.xpx,.xpy,Gstar,G_Z)
    } else {
        dblockpairs$G <- G[lower.tri(G, diag=TRUE)]
    }

    ## ###################################################
    ## now that G estimation is done, go back and force R to be PSD if needed.
    ## NOTE: now using nearPD2 function to maintain fixed zeros, and as with
    ## G*, there is a sequential decision about how to compute the adjusted matrix
    ## ###################################################
    if(!R_supplied){    
        names(dblockpairs)[which(names(dblockpairs)=="R")] <- "Rraw"
        dblockpairs$R <- dblockpairs$Rraw

        e <- eigen(R)
        if(any(e$values < -sqrt(.Machine$double.eps))){
            if(!control$quietly){
                cat("Adjusting R to make PSD...\n")
            }
            .R  <- R
            tmp <- nearPD2(R, fix0s=TRUE, do2eigen=FALSE, eig.tol=control$eig.tol, conv.tol=1e-11, maxit=1000, keepDiag=TRUE)
            if(tmp$converged){
                R <- tmp$mat
            } else {
                tmp <- nearPD2(R, fix0s=TRUE, do2eigen=FALSE, eig.tol=control$eig.tol, conv.tol=1e-11, maxit=1000, keepDiag=FALSE)
                if(tmp$converged){
                    R <- tmp$mat
                } else {
                    e$values[which(e$values < max(e$values)*control$eig.tol)] <- 0.0
                    R <- e$vectors %*% diag(e$values) %*% t(e$vectors)
                }
            }
            R   <- as(R,"sparseMatrix")
            rownames(R) <- colnames(R) <- .blocknames
            if(!control$quietly){
                cat(paste0("Smallest eigenvalue of adjusted R: ",min(eigen(R)$values),"\n"))
                cat("Summary of differences between original and adjusted R:\n")
                print(summary(c(as.matrix(R - .R))))
            }
            dblockpairs$R <- R[lower.tri(R, diag=TRUE)]

            ## adjust R_sb elements of dsch as well, since these are used in BLP
            for(s in 1:length(dsch)){
                x$R_sb     <- x$Ntilde * R
                x$R_sb     <- as(x$R_sb, "symmetricMatrix")
            }
        }
    }
    
    ## #############################################
    ## get GLS estimator of school FE, and BLPs of random effects,
    ## using brute-force estimation with sparse matrices
    ##
    ## NOTE: contribution of school*block FE treated as fixed and known
    ##
    ## NOTE: "Z" matrix in standard notation is I here because Y is
    ## directly additive in the school*block random effects
    ##
    ## NOTE: I tried the Henderson mixed model equations but found cases
    ## where this disagreed with the brute-force solution, possibly because
    ## the HMME require R^{-1} and G^{-1}, neither of which necessarily
    ## exist here, and it's possible the Schur results that the HMME rely
    ## on are not valid unless R and G are invertible.
    ##
    ## #############################################
    if(!mean_supplied){
        if(!control$quietly){
            cat("Computing GLS estimator and raw BLPs...\n")
        }
        
        stopifnot(all(sapply(dsch, function(x){ x$nblock == nrow(x$tab)})))
        N <- sum(sapply(dsch, function(x){ x$nblock }))
        d$sbid <- paste(d$school, gsub(" ","0",formatC(as.character(d$blockid))), sep="_")
        stopifnot(length(unique(d$sbid)) == N)

        ## block*pattern means which are treated as fixed offsets
        muhat_bp <- tapply(d$muhat_bp, d$sbid, mean)
        stopifnot(all(names(muhat_bp) == d$sbid[!duplicated(d$sbid)]))
        muhat_bp <- as.vector(muhat_bp)       

        ## "Y"
        Y <- matrix(as.vector(unlist(lapply(dsch, function(x){ x$tab$Y_sb_tilde }))), ncol=1)
        stopifnot( (nrow(Y) == N) && all(d$Y_sb_tilde[!duplicated(d$sbid)] == Y) )

        ## "X" design matrix for school FE only
        tmp  <- d[!duplicated(d$sbid),"schoolid",drop=FALSE]
        .X   <- sparse.model.matrix(~schoolid - 1, data=tmp, contrasts.arg=list(schoolid="contr.treatment"))
        stopifnot( (nrow(.X) == N) && all(gsub("schoolid","",colnames(.X)) == names(dsch)))

        ## "bigG" = G blocked to schools
        bigG  <- bdiag(lapply(dsch, function(x){
            .b <- x$oblocks
            as(G[.b,.b,drop=FALSE], "symmetricMatrix")
        }))

        ## V^{-1}
        Vinv <- bdiag(lapply(dsch, function(x){
            .b   <- x$oblocks
            .tmp <- solve(G[.b,.b,drop=FALSE] + x$R_sb[.b,.b,drop=FALSE])
            if(max(abs(.tmp - t(.tmp))) > 1e-8){
                warning("matrix of questionable symmetry arose in computation of V^{-1}")
            }
            .tmp <- forceSymmetric(.tmp)
            as(.tmp, "symmetricMatrix")
        }))
        
        ## brute force using expression for GLS estimator and BLPs of random
        ## effects at the school*block level
        .Xp_Vinv    <- crossprod(.X, Vinv)
        .Xp_Vinv_X  <- as(.Xp_Vinv %*% .X, "symmetricMatrix")
        .bgls       <- solve(.Xp_Vinv_X, (.Xp_Vinv %*% Y))
        for(s in 1:length(dsch)){
            dsch[[s]]$schoolFE  <- .bgls[s,1]
            dsch[[s]]$tab$muhat <- NULL
        }
        .Xbgls      <- .X %*% .bgls
        .blp        <- muhat_bp + as.vector(.Xbgls + (bigG %*% Vinv %*% (Y - .Xbgls)))

        ## summaries of BLUPs
        if(!control$quietly){
            cat(paste("(MM) range of raw BLPs:",round(min(.blp),digits=4),",",round(max(.blp),digits=4),"\n"))
            tmp <- as.vector(unlist(lapply(dsch, function(x){ x$tab$Y_sb })))
            cat(paste("(MM) cor(Y,raw BLPs)  :",round(as.vector(cor(.blp, tmp)),digits=4),"\n"))
        }

        ## MSE estimators, treating block*pattern FE as known but school FE as unknown
        if(!control$quietly){
            cat("Computing MSE estimates of raw BLPs by school (may be slow)...\n")
        }
        .pos <- 0
        for(s in 1:length(dsch)){
            x    <- dsch[[s]]
            nb   <- x$nblock
            .i   <- (.pos+1):(.pos+nb)
            .pos <- max(.i)

            x$tab$blp <- .blp[.i]

            .vinv <- Vinv[.i,.i,drop=F]
            .G    <- bigG[.i,.i,drop=F]
            H     <- (matrix(1.0, ncol=nb, nrow=nb) %*% .vinv) / sum(.vinv)
            ImH   <- diag(nb) - H
            Q     <- H + ((.G %*% .vinv) %*% ImH) ## NOTE: rowsums==1, useful for weights
            ImQ   <- diag(nb) - Q

            x$mse_blp <- (ImQ %*% .G %*% t(ImQ)) + (Q %*% x$R_sb[x$oblocks,x$oblocks,drop=F] %*% t(Q))
            x$Q       <- Q
            x$var_schoolFE <- 1.0/sum(.vinv) ## variance of GLS estimator of school FE = (1'V^{-1}1)^{-1}
            dsch[[s]] <- x
        }
        stopifnot(.pos == N)
        rm(.blp,Y,.Xbgls,bigG,Vinv,.X,.Xp_Vinv_X,.Xp_Vinv); gc()

        ## use schoolFE and var_schoolFE to estimate variance in individual
        ## growth attributable to schoolFE, where schools are weighted according
        ## to the total number of attached growth scores, and where we assume
        ## that errors for estimated FE by school are uncorrelated, with
        ## is an approximation
        .l    <- sapply(dsch, function(x){ sum(x$tab$nsb) })
        .l    <- .l/sum(.l)
        .y    <- sapply(dsch, function(x){ x$schoolFE })
        .s    <- sapply(dsch, function(x){ x$var_schoolFE })
        modstats["estimated_variance_among_schools"] <- (sum(.l * .y^2) - sum(.l * .s)) - ( (sum(.l * .y))^2 - sum(.l^2 * .s) )
        modstats["estimated_percvar_among_schools"]     <- modstats["estimated_variance_among_schools"] / modstats["varY"]
    }
    
    ## #############################################
    ## BLP and MSE calculation for weighted composite
    ## #############################################
    if(!control$quietly){
        cat("Computing composite BLPs and MSE estimates...\n")
    }
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

            ## create matrix to store weights for direct and BLP estimators.
            ## also stored indicator "obs" which indicates whether there
            ## are observed data for the corresponding block.
            weights           <- matrix(0, ncol=3, nrow=B)
            rownames(weights) <- .blocknames
            colnames(weights) <- c("obs","direct","blp")
            weights[b,"obs"]    <- 1
            weights[b,"direct"] <- lambda
            
            ## get BLP, either by calling blp() function for known mean, or combining the pieces
            ## from the GLS procedure
            if(mean_supplied){
                tmp <- blp(x$tab$Y_sb, lambda, x$tab$muhat, as.matrix(G[b,b,drop=F]), as.matrix(x$R_sb[b,b,drop=F]), eig.tol=control$eig.tol)
                .pi <- as.vector(x$pis[b])
                .Pi <- matrix(.pi, nrow=x$nblock, ncol=x$nblock, byrow=T)
                weights[b,"blp"] <- as.vector( t(lambda) %*% ( (tmp$IminusQ %*% .Pi) + tmp$Q ) )
            } else{
                est.direct   <- sum(lambda * x$tab$Y_sb)
                mse.direct   <- as.vector(t(lambda) %*% x$R_sb[b,b,drop=F] %*% lambda)
                est.blp      <- sum(lambda * x$tab$blp)
                mse.blp      <- as.vector(t(lambda) %*% x$mse_blp %*% lambda)
                mse.null     <- as.vector(t(lambda) %*% G[b,b,drop=F] %*% lambda)
                prmse.null   <- 1 - mse.blp/mse.null
                prmse.direct <- 1 - mse.blp/mse.direct
                tmp <- list(est=c(est.direct = est.direct, mse.direct = mse.direct, est.blp = est.blp, mse.blp = mse.blp, prmse.null = prmse.null, prmse.direct = prmse.direct))
                weights[b,"blp"] <- as.vector( t(lambda) %*% x$Q )
            }
            
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
                bhat_ols               = .bhat,
                modstats               = modstats,
                tab_patterns           = tab_patterns,
                G                      = G,
                R                      = R,
                target_blocks          = target_blocks,
                target_contrast_blocks = target_contrast_blocks,
                adjusted_growth        = do.call("rbind",lapply(dsch, function(x){ x$est }))))
}
