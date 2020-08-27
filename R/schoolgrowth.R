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
        
    if(is.null(control$pattern_nmin)){
        control$pattern_nmin <- 100
    }

    if(is.null(control$alpha_zero)){
        control$alpha_zero <- FALSE
    }
    
    if(is.null(control$blockpair_student_nmin)){
        control$blockpair_student_nmin <- 100
    }

    if(is.null(control$return_d)){
        control$return_d <- FALSE
    }

    if(is.null(control$mse_blp_chk)){
        control$mse_blp_chk <- FALSE
    }

    if(is.null(control[["jackknife"]])){
        control$jackknife <- TRUE
    }

    if(is.null(control$return_schjack)){
        control$return_schjack <- TRUE
    }
        
    if(is.null(control$patterns_only)){
        control$patterns_only <- FALSE
    }

    R_supplied    <- !is.null(control[["R"]])
    G_supplied    <- !is.null(control[["G"]])
    
    ## Radj control parameters
    if(!R_supplied){
        if(is.null(control$Radj_method)){
            control$Radj_method <- "nearPD"
        }
        if(!(control$Radj_method %in% c("nearPD","spectral"))){
            stop("Invalid specification of control$Radj_method")
        }
        if(is.null(control$Radj_eig_tol)){
            control$Radj_eig_tol <- 1e-10
        }
        if( (control$Radj_method == "spectral") && is.null(control$Radj_eig_min) ){
            control$Radj_eig_min <- 0.0
        }
        if( (control$Radj_method == "spectral") && (control$Radj_eig_tol < control$Radj_eig_min) ){
            stop("control$Radj_eig_tol must be greater than or equal to control$Radj_eig_min")
        }
    }

    ## Gadj control parameters
    if(!G_supplied){
        if(is.null(control$Gadj_method)){
            control$Gadj_method <- "nearPD"
        }
        if(!(control$Gadj_method %in% c("nearPD","spectral"))){
            stop("Invalid specification of control$Gadj_method")
        }
        if(is.null(control$Gadj_eig_tol)){
            control$Gadj_eig_tol <- 1e-10
        }
        if( (control$Gadj_method == "spectral") && is.null(control$Gadj_eig_min) ){
            control$Gadj_eig_min <- 0.0
        }
        if( (control$Gadj_method == "spectral") && (control$Gadj_eig_tol < control$Gadj_eig_min) ){
            stop("control$Gadj_eig_tol must be greater than or equal to control$Gadj_eig_min")
        }
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

    ## checks on user-supplied R
    if(R_supplied){
        R <- control$R
        if(!control$quietly){
            cat("Checking user-supplied value of R...\n")
        }

        if(class(R) != "dsCMatrix"){
            stop("R must be a member of class 'dsCMatrix' from the Matrix library")
        }

        if(R@uplo != "L"){
            stop("R@uplo must be 'L'")
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
    
    ## checks on user-supplied G
    if(G_supplied){
        G <- control$G
        if(!control$quietly){
            cat("Checking user-supplied value of G...\n")
        }

        if(class(G) != "dspMatrix"){
            stop("G must be a member of class 'dspMatrix' from the Matrix library")
        }

        if(G@uplo != "L"){
            stop("G@uplo must be 'L'")
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

        if( any(abs(apply(G, 1, sum)) > 1e-10) ){
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
    } else if(target["years"] == "all"){
        w1 <- rep(TRUE, B)
    } else{
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
        w2 <- rep(TRUE, B)
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
        w3 <- rep(TRUE, B)
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
    target_blockids <- dblock$blockid[dblock$intarget]
    
    ## ########################################################
    ## parse "target_contrast"
    ## ########################################################    
    if(!is.null(target_contrast)){
        
        if(!is.character(target_contrast) || is.null(names(target_contrast)) || (length(target_contrast) != 4)){
            stop("'target_contrast' must be a named character vector of length 4")
        }
        
        if(!all(sort(names(target_contrast)) == c("grades","subjects","weights","years"))){
            stop("invalid names of 'target_contrast'; names must be 'years','subjects','grades','weights'")
        }
        
        ## years:
        if(target_contrast["years"] == "final"){
            w1 <- (dblock$year == .finalyear)
        } else if(target_contrast["years"] == "all"){
            w1 <- rep(TRUE, B)
        } else{
            .years <- gsub(" ","", unlist(strsplit(target_contrast["years"],",")))
            if(any(is.na(.years)) || !all(.years %in% dblock$year)){
                stop("'target_contrast' specifies years that do not occur in data")
            }
            w1 <- dblock$year %in% .years
        }
        
        ## subjects:
        if(target_contrast["subjects"] == "all"){
            w2 <- rep(TRUE, B)
        } else {
            .subjects <- gsub(" ","", unlist(strsplit(target_contrast["subjects"],",")))
            if(any(is.na(.subjects)) || !all(.subjects %in% dblock$subject)){
                stop("'target_contrast' specifies subjects that do not occur in data")
            }
            w2 <- (dblock$subject %in% .subjects)
        }
        
        ## grades:
        if(target_contrast["grades"] == "all"){
            w3 <- rep(TRUE, B)
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
        target_contrast_blockids <- dblock$blockid[dblock$intarget_contrast]
    }
    
    ## #################################################################
    ## create "dblockpairs" dataframe that tracks block pairs and associated data
    ##
    ## NOTE: don't use combn() because we also want the diagonals and we
    ## need them in the correct positions
    ##
    ## NOTE: elements of dblockpairs ordered to fill the lower triangle
    ## of a BxB symmetric matrix in column-major order, which is why
    ## the outer loop is over columns and inner loop over rows
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
    dblockpairs        <- dblockpairs[,c("blocki","blockj","blockidi","blockidj","iB")]
    ## CHECK on order:
    ## print(sparseMatrix(i=dblockpairs$blockidi, j=dblockpairs$blockidj, x=dblockpairs$iB, symmetric=TRUE))
    
    ## ##################################################################
    ## create "dsch" master list that will hold all school-level information.
    ## ##################################################################
    tmp  <- sort(unique(d$school))
    S    <- length(tmp)
    dsch <- vector(S, mode="list")
    names(dsch) <- tmp
    for(s in 1:S){
        dsch[[s]]$school <- tmp[s]
    }
    
    ## ###################################################################
    ## put number of schools in each block pair onto dblockpairs
    ## ###################################################################
    if(!control$quietly){
        cat("Computing counts of schools in each block pair...\n")
    }
    .tab <- table(d$school, d$blockid)
    stopifnot(all(colnames(.tab) == 1:B) && all(rownames(.tab) == names(dsch)))
    dblockpairs$nsch  <- 0L
    
    for(wh in 1:B2){
        dblockpairs$nsch[wh]  <- as.integer(sum( (.tab[,dblockpairs$blockidi[wh]] * .tab[,dblockpairs$blockidj[wh]]) >= 0.99 ))
    }

    ## stop if any block has too few schools for G estimation
    if(!G_supplied){
        if(!all(dblockpairs$nsch >= 1L)){
            stop("For G estimation, all pairs of blocks must have at least 1 school with growth scores in both blocks")
        }
    }
    
    ## #####################################################################
    ## compute "N" for each school which is a sparse symmetric matrix that gives
    ## the number of students in each block pair.  The diagonals are the
    ## number of students in a given school contributing to the growth measure
    ## for each block, and off-diagonals are the number of students in a given
    ## school contributing to the growth measures in a pair of blocks
    ## #####################################################################
    if(!control$quietly){
        cat("Computing counts of students in each block pair by school...\n")
    }
    tmp <- split(d[,c("school","stuid","blockid")], d$school)
    stopifnot(all(names(tmp) == names(dsch)))
    
    for(s in 1:S){
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
        trips <- as.data.frame(do.call("rbind", trips))
        trips <- trips[which(trips$n > 0),]
        dsch[[s]]$N <- sparseMatrix(i=trips$i, j=trips$j, x=trips$n, dims=c(B,B), symmetric=TRUE)
        stopifnot(max(abs(as.matrix(dsch[[s]]$N[oblocks,oblocks]) - crossprod(.tab))) < 1e-10)
    }
    stopifnot(sum(sapply(dsch, function(x){ sum(diag(x$N))})) == nrow(d))
    rm(.tab,tmp); gc()
    
    ## ##################################################################
    ## define student pattern variable, where "pattern" refers to whether
    ## a gain is observed, or not observed, in each block.  gets counts of
    ## patterns and then combine small patterns until all patterns have
    ## at least control$pattern_nmin students.
    ##
    ## NOTE: if control$pattern_nmin="min", it will find and use the smallest
    ## possible value of control$pattern_nmin such that there is no
    ## stratification.
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
        dblockpairs$nstu[wh]  <- as.integer(sum( stupat[,dblockpairs$blockidi[wh]] * stupat[,dblockpairs$blockidj[wh]] ))
        dblockpairs$R_est[wh] <- dblockpairs$nstu[wh] >= control$blockpair_student_nmin
    }

    ## stop if any block has too few students for R estimation
    if(!R_supplied){
        .locs <- which( dblockpairs$blockidi == dblockpairs$blockidj )        
        if(!all(dblockpairs$R_est[.locs])){
            stop("At least one block has fewer than control$blockpair_student_nmin students")
        }
    }

    ## reformat stupat
    stupat           <- data.frame(stuid = rownames(stupat), pattern = apply(as.matrix(stupat), 1, paste, collapse=""), stringsAsFactors=FALSE)
    stupat$stuid     <- as.character(stupat$stuid)
    stupat$pcount    <- ave(rep(1,nrow(stupat)), stupat$pattern, FUN=sum)
    if(!control$quietly){
        cat(paste("Number of patterns before collapsing:",length(unique(stupat$pattern)),"\n"))
    }
    
    ## find the minimum pattern count such that maintaining all patterns
    ## with at least that count ensures no stratification of
    ## schools * (block*pattern indicators)
    .uc       <- as.integer(c(sort(unique(stupat$pcount)), nrow(d)+1))
    .pos      <- 0L
    connected <- FALSE
    while(!connected){
        .pos             <- .pos + 1L
        stupat$cpattern  <- stupat$pattern
        stupat$cpattern[which(stupat$pcount < .uc[.pos])] <- "collapsed"
        stupat$patternid <- as.integer(as.factor(stupat$cpattern))
        d$patternid      <- stupat$patternid[match(d$stuid, stupat$stuid)]    

        d$schoolid <- as.integer(as.factor(d$school))
        d$bpid     <- as.integer(as.factor(paste(d$blockid, d$patternid)))
    
        grp1.set <- d$schoolid[1]
        grp2.set <- unique(d$bpid[d$schoolid %in% grp1.set ])
        l.old    <- c(length(grp1.set), length(grp2.set))
    
        grp1.set <- unique(d$schoolid[ d$bpid %in% grp2.set ])
        grp2.set <- unique(d$bpid[ d$schoolid %in% grp1.set ])
        
        while( max( c(length(grp1.set), length(grp2.set)) - l.old ) > 0 ){
            l.old    <- c(length(grp1.set), length(grp2.set))
            grp1.set <- unique(d$schoolid[ d$bpid %in% grp2.set ])
            grp2.set <- unique(d$bpid[ d$schoolid %in% grp1.set ])
        }
        connected <- (length(unique(d$schoolid)) == length(grp1.set)) && (length(unique(d$bpid)) == length(grp2.set))
    }
    smallest_pattern_nmin <- .uc[.pos]
    if(!control$quietly){
        cat(paste("Minimum pattern count required for no stratification:",smallest_pattern_nmin,"\n"))
    }
    d$bpid <- NULL
    stupat <- stupat[,c("stuid","pattern","pcount")]
    
    ## create pattern groupings depending on control$pattern_nmin, using collapsing if needed
    if(control$pattern_nmin=="min"){
        control$pattern_nmin <- smallest_pattern_nmin
        if(!control$quietly){
            cat(paste("Setting control$pattern_nmin to",smallest_pattern_nmin,"\n"))
        }
    }

    if(control$pattern_nmin < smallest_pattern_nmin){
        stop(paste("control$pattern_nmin is too small; data indicate it must be at least",smallest_pattern_nmin))
    }

    ## apply collapsing rule if needed:
    ## 1) collapse all patterns with counts less than control$pattern_nmin into one pattern
    ## 2) if this combined pattern has count still below threshold, keep pulling in the
    ##    least-populated patterns until it crosses
    tmp <- stupat[which(!duplicated(stupat$pattern)),c("pattern","pcount")]
    stopifnot(nrow(tmp) == length(unique(stupat$pattern)))
    tmp <- tmp[order(-tmp$pcount),]
    tmp$cpattern <- tmp$pattern
    tmp$cpattern[which(tmp$pcount < control$pattern_nmin)] <- "collapsed"
    if(any(tmp$cpattern=="collapsed")){
        if(!all(tmp$cpattern=="collapsed")){
            .w <- which(tmp$cpattern == "collapsed")
            .n <- sum(tmp$pcount[.w])
            while(.n < control$pattern_nmin){
                nextsmallest <- .w[1]-1
                tmp$cpattern[nextsmallest] <- "collapsed"
                .n <- .n + tmp$pcount[nextsmallest]
                .w  <- c(nextsmallest,.w)
            }
        }
    }
    if(!control$quietly){
        cat(paste("Number of patterns after collapsing (control$pattern_nmin=",control$pattern_nmin,"): ", length(unique(tmp$cpattern)),"\n",sep=""))
    }
    tab_patterns <- tmp
    rownames(tab_patterns) <- 1:nrow(tab_patterns)

    if(control$patterns_only){
        return(tab_patterns)
    }
    
    ## assign the collapsed patterns to "collapsed"
    if(any(tmp$cpattern=="collapsed")){
        .w <- which(tmp$cpattern=="collapsed")
        stupat$pattern[which(stupat$pattern %in% tmp$pattern[.w])] <- "collapsed"
        stopifnot(sum(stupat$pattern=="collapsed") == sum(tmp$pcount[.w]))
    }
    
    ## create pattern ID and merge onto d; create "bpid" = block*pattern indicator
    stupat$patternid <- as.integer(as.factor(stupat$pattern))
    stopifnot(length(unique(stupat$patternid)) == length(unique(tmp$cpattern)))
    d$patternid <- stupat$patternid[match(d$stuid, stupat$stuid)]
    rm(stupat); gc()
    d$bpid      <- as.integer(as.factor(paste(d$blockid, d$patternid)))
    
    ## #######################################################
    ## Estimate R using within-school, within-block, within-pattern deviations
    ## ########################################################
    if(!R_supplied){
        d$sbp  <- as.integer(as.factor(paste(d$school, d$blockid, d$patternid)))
        d$nsbp <- ave(rep(1,nrow(d)), d$sbp, FUN=sum)
        d$e    <- d$Y - ave(d$Y, d$sbp)
        
        dblockpairs$R      <- 0.0
        dblockpairs$R_nstu <- 0L

        if(!control$quietly){
            cat("Estimating R (variances)...\n")
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
            cat("Estimating R (covariances)...\n")
        }

        .locs <- which( (dblockpairs$blockidi != dblockpairs$blockidj) & (dblockpairs$R_est) )
        for(wh in dblockpairs$iB[.locs]){
            bi  <- dblockpairs$blockidi[wh]
            bj  <- dblockpairs$blockidj[wh]
            tmp <- d[which(d$blockid %in% c(bi, bj)), c("stuid","school","blockid","sbp","nsbp","e")]
            tmp$blockid[which(tmp$blockid==bi)] <- 0
            tmp$blockid[which(tmp$blockid==bj)] <- 1

            .dups <- unique(tmp$stuid[duplicated(tmp$stuid)])
            tmp   <- tmp[which(tmp$stuid %in% .dups),]

            if(nrow(tmp) > 0){
                stopifnot(all(as.data.frame(table(tmp$stuid))$Freq == 2))
                tmp <- reshape(tmp, timevar="blockid", idvar="stuid", direction="wide")
                ## restrict to sbp with at least two observations
                tmp <- tmp[which( (tmp$nsbp.0 >= 2) & (tmp$nsbp.1 >= 2) ),]
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
        dblockpairs$R             <- as.matrix(R)[lower.tri(R, diag=TRUE)]
        dblockpairs$R_nstu        <- 0L
        is.na(dblockpairs$R_nstu) <- TRUE
        dblockpairs$R_est         <- FALSE
    }
    
    ## ###################################################################
    ## compute OLS estimates of (alpha, mu)
    ## ###################################################################
    if(!control$quietly){
        cat("Computing OLS estimates of regression coefficients...\n")
    }
    modstats   <- list(ntot = nrow(d), nstu = length(unique(d$stuid)), nsch = length(unique(d$school)), varY = var(d$Y))
    d$schoolid <- factor(d$school)
    d$bpid     <- factor(d$bpid)

    if(!control$alpha_zero){
        ## fit model with school FE and block/pattern FE
        .mdf    <- length(unique(d$schoolid)) + length(unique(d$bpid)) - 1
        .X      <- sparse.model.matrix(~schoolid - 1 + bpid, data=d, contrasts.arg=list(schoolid="contr.treatment",bpid="contr.sum"))
        stopifnot(all(.X@x %in% c(-1,0,1)) && (ncol(.X) == .mdf))
        stopifnot(length(grep("schoolid",colnames(.X))) == length(unique(d$school)))
        stopifnot(length(grep("bpid",    colnames(.X))) == length(unique(d$bpid))-1)
        .xpx    <- crossprod(.X)
        .xpy    <- crossprod(.X, d$Y)
        .bhat   <- solve(.xpx, .xpy)
        tmp     <- as.vector(.X %*% .bhat)
        stopifnot(max(abs(tapply(tmp, d$school, mean) - tapply(d$Y, d$school, mean))) < 1e-6)
        stopifnot(max(abs(tapply(tmp, d$bpid,   mean) - tapply(d$Y, d$bpid,   mean))) < 1e-6)
        
        ## add estimated means based on block/pattern FE (but not including the school fixed effects),
        ## which are treated as fixed and known for the remainder of the estimation steps
        wh           <- grep("bpid", colnames(.X))
        d$alpha      <- as.vector(.X[,wh] %*% .bhat[wh])
        .bhat        <- as.vector(.bhat)
        names(.bhat) <- colnames(.X)
    
        ## put estimated school FE on dsch (provisional for now; replaced with GLS estimator later)
        tmp    <- d[!duplicated(d$school),c("school","schoolid")]
        tmp$mu <- .bhat[paste0("schoolid",tmp$schoolid)]
        stopifnot( all(!is.na(tmp$mu)) )
        for(s in 1:S){
            dsch[[s]]$mu <- tmp$mu[which(tmp$school == dsch[[s]]$school)]
        }
        rm(.X,.xpx,.xpy,tmp); gc()
    } else {
        .bhat   <- NULL
        d$alpha <- 0.0
        ## when alpha==0, the provisional estimator of mu is just sample mean of Y by school
        for(s in 1:S){
            dsch[[s]]$mu <- mean(d$Y[which(d$school == dsch[[s]]$school)])
        }
    }
        
    ## ############################################################
    ## calculate and check various school*block aggregates, which will be used
    ## for estimating G and ultimately for doing the EBLP calculation
    ## ############################################################
    if(!control$quietly){
        cat("Computing and checking school*block aggregate measures...\n")
    }
    d$alpha_sb   <- ave(d$alpha,          d$school, d$block)
    d$Y_sb       <- ave(d$Y,              d$school, d$block)
    d$Y_sb_tilde <- d$Y_sb - d$alpha_sb
    d$nsb        <- ave(rep(1,nrow(d)), d$school, d$block, FUN=sum)
    
    tmp        <- d[!duplicated(paste(d$school, d$block)), c("school","grade","year","subject","block","blockid","nsb","alpha_sb","Y_sb","Y_sb_tilde")]
    tmp        <- tmp[order(tmp$school, tmp$blockid),]
    stopifnot(nrow(tmp) == length(unique(paste(d$school, d$block))))
    tmp        <- split(tmp, tmp$school)
    stopifnot(all(sapply(tmp, function(x){ all(diff(x$blockid) > 0) })) && all(names(tmp) == names(dsch)))
    for(s in 1:S){
        dsch[[s]]$tab <- tmp[[s]]
    }

    stopifnot(all(unlist(lapply(dsch, function(x){ diag(x$N)[x$tab$blockid] - x$tab$nsb })) == 0))
    stopifnot(max(abs(sapply(dsch, function(x){ x$mu - weighted.mean(x$tab$Y_sb_tilde, w = x$tab$nsb) }))) < modstats[["varY"]]*1e-8)

    ## ########################################################
    ## if jackknife, get number of schools contributing to G estimation (those with
    ## nblock > 1), and among those schools, assign which jackknife batch each will
    ## be excluded from.  default number of batches is min(#contributing schools,
    ## 50), or use user-supplied J
    ## #########################################################
    if(control$jackknife){
        Sj <- sum(sapply(dsch, function(x){ x$nblock > 1 }))

        ## determine number of batches
        if(!is.null(control$jackknife_J)){
            if(control$jackknife_J == "max"){
                J <- Sj
            } else {
                J <- as.integer(control$jackknife_J)
                if(is.na(J)){
                    stop("invalid specification of control$jackknife_J")
                }
                if(J < 2){
                    stop("control$jackknife_J is too small")
                }
                if(J > Sj){
                    stop("control$jackknife_J is too large")
                }
            }
            if(J > 100){
                warning("control$jackknife_J is large, which may result in excessive RAM usage")
            }
        } else {
            J  <- min(50, Sj)
        }

        ## determine number of schools excluded from each batch
        .n <- rep(0, J)
        .i <- .cumi <- 1L
        while(.cumi <= Sj){
            .n[.i] <- .n[.i] + 1L
            .i     <- .i + 1L
            if(.i > J){
                .i <- 1L
            }
            .cumi <- .cumi + 1L
        }
        stopifnot(Sj == sum(.n))

        ## randomly assign each school to the batch that will exclude it
        ## (schools with only one block are nominally assigned batch 0)
        tmp <- sample(rep(as.integer(1:J),times=.n))
        .i      <- 0
        for(s in 1:S){
            if(dsch[[s]]$nblock == 1){
                dsch[[s]]$jackex <- 0L
            } else {
                .i <- .i + 1
                dsch[[s]]$jackex <- tmp[.i]
            }
        }
    } else {
        J <- 0L
    }
    
    ## ########################################################
    ## compute provisional variance/covariance matrix "R_sb" of errors in aggregate
    ## measures, as well as pieces needed for the contribution of the school to
    ## the WLS estimator for the elements of G.
    ##
    ## NOTE: do this with a loop and accumulate key results so that we don't need
    ## to store large pieces for each school, and can also accumulate pieces
    ## needed for jackknife
    ## ########################################################
    if(!control$quietly){
        cat("Computing required school-level quantities (may be slow when estimating G)...\n")
    }

    ## create R
    ## NOTE: this is provisional, used for moment estimation, and may not be PSD.
    if(!R_supplied){
        tmp <- dblockpairs[which(dblockpairs$R_est),]
        R   <- sparseMatrix(i=tmp$blockidi, j=tmp$blockidj, x=tmp$R, dims=c(B,B), symmetric=TRUE)
        rownames(R) <- colnames(R) <- .blocknames
    }

    ## create matrix to implement sum-to-zero constraints
    A  <- sparseMatrix(i=c(1:(B-1), rep(B,B-1)), j=c(1:(B-1),1:(B-1)), x=c(rep(1,B-1),rep(-1,B-1)), dims=c(B,B-1))

    ## create matrix to hold the sum across schools of each
    ## school's observed cross-product matrix, adjusted for measurement error.
    ## We store as a (BxB) symmetric matrix.
    ##
    ## Also create matrix to hold the design matrix for the WLS estimator of G.
    ## there are B(B+1)/2 rows, corresponding to the column-major lower triangle
    ## of the observed moments, and there are (B-1)B/2 columns, corresponding to the
    ## column-major lower triangle of G*, the (B-1)x(B-1) full-rank covariance matrix.
    ##
    ## NOTE: sparse tended to be slower here so we make it dense, which is fine
    ## because the matrix actually ends up dense.  Zsum will accumulate across
    ## schools and each school's contribution will be initialized to Zs0.
    ##
    ## NOTE: expanded to list to accommodate jackknife calculations.  the
    ## first element of the list is for the full data, and the other elements
    ## are for jackknife samples
    Zs0  <- matrix(0.0, nrow=B*(B+1)/2, ncol=(B-1)*B/2)
    .tmp <- sparseMatrix(i=1,j=1,x=0, dims=c(B,B), symmetric=TRUE)
    Zsum <- Ysum <- vector(J+1,mode="list")
    for(j in 1:(J+1)){
        Ysum[[j]] <- .tmp
        Zsum[[j]] <- Zs0
    }
    
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
    for(s in 1:S){
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

        ## compute school contribution to Ysum, for schools with 2+ blocks
        if(!G_supplied && (x$nblock > 1)){
            Pis <- sparseMatrix(i=1,j=1,x=0, dims=c(B,B))
            for(.b in b){
                Pis[.b,] <- x$pis
            }
            Is        <- sparseMatrix(i=b,j=b,x=rep(1,x$nblock), dims=c(B,B), symmetric=TRUE)
            me_adj    <- (Is - Pis) %*% x$R_sb %*% t(Is - Pis)
            .y        <- sparseMatrix(i = rep(1,x$nblock), j = b, x = x$tab$Y_sb_tilde - x$mu, dims=c(1,B))
            .piece    <- crossprod(.y) - me_adj
            Ysum[[1]] <- Ysum[[1]] + .piece
            if(control$jackknife){
                ## NOTE: need x$jackex + 1 because j=2...J+1 are positions of jackknife sums
                Ysum[[x$jackex + 1]] <- Ysum[[x$jackex + 1]] + .piece
            }
            
            ## NOTE: equivalent calculation:
            ##
            ## .R  <- (Is - Pis) %*% sparseVector(x$tab$Y_sb_tilde, i = b, length=B)
            ## print(max(abs(crossprod(.y) - ( (.R %*% t(.R))))))

            ## ##################################################
            ## CHECK ME adjustment (method 1)
            ## 
            ## v    <- x$tab$nsb/sum(x$tab$nsb)
            ## Pis  <- matrix(v, nrow=x$nblock, ncol=x$nblock, byrow=T)
            ## .S    <- x$R_sb[b,b]
            ## me_adj1 <- .S - (Pis %*% .S) - (.S %*% t(Pis)) + (Pis %*% .S %*% t(Pis))
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
            Zsum[[1]] <- Zsum[[1]] + Zs
            if(control$jackknife){
                Zsum[[x$jackex + 1]] <- Zsum[[x$jackex + 1]] + Zs
            }
            
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
    rm(Zs0); gc()

    ## adjust jackknife sums to be the total minus the relevant sum for excluded schools,
    ## so that the result equals the appropriate sum for the schools that are included
    ## in each jackknife batch
    if(control$jackknife){
        for(j in 2:(J+1)){
            Ysum[[j]] <- Ysum[[1]] - Ysum[[j]]
            Zsum[[j]] <- Zsum[[1]] - Zsum[[j]]
        }
    }
    
    ## #########################################################
    ## compute WLS estimates of G* elements, translate to G, and
    ## force to PSD. This is done for the base sample and (if applicable)
    ## all jackknife samples, storing everything in a list "G" of matrices
    ## #########################################################
    if(G_supplied){
        dblockpairs$G <- G@x
        G             <- list(G)
    } else {
        if(!control$quietly){
            cat("Estimating G...\n")
        }
        G <- Gstar <- vector(J+1, mode="list")

        ## weight matrix .W for GLS estimation
        .W  <- sparseMatrix(i=1:B2, j=1:B2, x=dblockpairs$nsch, dims=c(B2,B2), symmetric=TRUE)

        ## #########################################################
        ## compute Gstar for full sample and each jackknife sample.
        ## for each element of Gstar, compute G and force to PSD
        ## #########################################################        
        for(j in 1:(J+1)){
            .Y     <- as.matrix(Ysum[[j]])
            .Y     <- .Y[lower.tri(.Y, diag=TRUE)]

            ## solve for estimable parameters and construct Gstar and G
            .xpx       <- t(Zsum[[j]]) %*% .W %*% Zsum[[j]]
            .xpy       <- t(Zsum[[j]]) %*% .W %*% .Y
            Gstar[[j]] <- new("dspMatrix", Dim=c(B-1L,B-1L), uplo="L", x=as.vector(solve(.xpx, .xpy)))
            .G         <- as.matrix(A %*% Gstar[[j]] %*% t(A))
            .G         <- new("dspMatrix", Dim=c(B,B), uplo="L", x=.G[lower.tri(.G, diag=TRUE)])
            ## PSD adjustment
            G[[j]]     <- Gadj(.G, control$Gadj_method, control$Gadj_eig_tol, control$Gadj_eig_min)
            rownames(G[[j]]) <- colnames(G[[j]]) <- .blocknames

            ## add full-sample pieces to dblockpairs
            if(j==1){
                dblockpairs$Graw <- .G@x
                dblockpairs$G    <- G[[j]]@x
            }
        }
        rm(.Y, .W, .xpx, .xpy, Zsum, Ysum, .G)
        
        ## ############################################################################        
        ## if control$jackknife, use jackknife samples to estimate variance of G
        ## ############################################################################
        if(control$jackknife){
            .m <- as.matrix(Reduce("+",G[-1]) / J)
            varhat_G <- ((J-1)/J) * Reduce("+", lapply(G[-1], function(x){ (as.matrix(x) - .m)^2 }))
        }
    }
    
    ## ###################################################
    ## now that G estimation is done, go back and force R to be PSD if needed,
    ## and adjust R_sb elements of dsch as well
    ## ###################################################
    if(!R_supplied){    
        names(dblockpairs)[which(names(dblockpairs)=="R")] <- "Rraw"
        nzlocs        <- as.matrix(dblockpairs[dblockpairs$R_est, c("blockidi","blockidj")])
        stopifnot(nrow(nzlocs) == length(R@x))
        R             <- Radj(R, nzlocs, control$Radj_method, control$Radj_eig_tol, control$Radj_eig_min)
        rownames(R)   <- colnames(R) <- .blocknames
        dblockpairs$R <- as.matrix(R)[lower.tri(R, diag=TRUE)]

        for(s in 1:S){
            dsch[[s]]$R_sb <- dsch[[s]]$Ntilde * R
            dsch[[s]]$R_sb <- as(dsch[[s]]$R_sb, "symmetricMatrix")
        }
    }
    
    ## #############################################
    ## get GLS estimator of school FE, and EBLPs of block means,
    ## treating block*pattern contributions for each school as fixed
    ##
    ## NOTE: "Z" matrix in standard notation is I here because Y is
    ## directly additive in the school*block random effects, and
    ## "X" matrix in standard notation is just a vector of 1.
    ##
    ## NOTE: (4/28/2020) reduce cost of Matrix overhead using
    ## ".Gm" and ".Rm" to try to speed up calculations with large S*J
    ## #############################################
    if(!control$quietly){
        cat("Computing GLS estimators, raw EBLPs, and MSEs (may be slow, especially with jackknife)...\n")
    }

    for(s in 1:S){
        x     <- dsch[[s]]
        nb    <- x$nblock
        b     <- x$oblocks
        
        if(control$jackknife){
            x$jack_blp  <- matrix(nrow=nb, ncol=J)
            x$jack_mu   <- rep(-99.0, length=J)
            x$jack_vinv <- matrix(nrow=3, ncol=J)
            rownames(x$jack_vinv) <- c("smallest","largest","trace")
        }
        
        .Y    <- matrix(x$tab$Y_sb_tilde, ncol=1)
        .R    <- x$R_sb[b,b,drop=FALSE]
        .Rm   <- as.matrix(.R)

        for(j in 1:(J+1)){
            ## GLS estimator of mu, plus EBLP at school*block level, for this "j"
            .Gm   <- as.matrix(G[[j]][b,b,drop=FALSE])
            .vinv <- solve(.Gm + .Rm)
            .mu   <- sum(as.vector(.vinv %*% .Y)) / sum(.vinv)
            .blp  <- x$tab$alpha_sb + as.vector(.mu + (.Gm %*% .vinv %*% (.Y - .mu)))

            if(j==1){
                x$mu      <- .mu
                x$tab$blp <- .blp
                
                ## MSE estimators, treating alpha_sb as known, and not accounting for jackknife
                .G    <- as(G[[j]][b,b,drop=FALSE], "symmetricMatrix")
                H     <- (matrix(1.0, ncol=nb, nrow=nb) %*% .vinv) / sum(.vinv)
                ImH   <- diag(nb) - H
                Q     <- H + ((.G %*% .vinv) %*% ImH) ## NOTE: rowsums==1, useful for weights
                ImQ   <- diag(nb) - Q
                x$mse_blp <- (ImQ %*% .G %*% t(ImQ)) + (Q %*% .R %*% t(Q))
        
                if(control$mse_blp_chk){
                    ## check MSE using alt formula from Das, Jiang and Rao (2004)
                    g1  <- .G - (.G %*% .vinv %*% .G)
                    tmp <- diag(nb) - (.vinv %*% .G)
                    g2  <- t(tmp) %*% (matrix(1.0, ncol=nb, nrow=nb) / sum(.vinv)) %*% tmp
                    mse_blp_chk <- g1 + g2
                    if(max(abs(x$mse_blp - mse_blp_chk)) > (modstats[["varY"]]*(1e-8))){
                        if(!control$quietly){
                            print(max(abs(x$mse_blp - mse_blp_chk)))
                        }
                        stop(paste0("mse_blp_chk failed: school ",s))
                    }
                }
                
                x$Q         <- Q
                
                ## variance of GLS estimator of school FE = (1'V^{-1}1)^{-1}, conditional on known V        
                x$var_muhat <- 1.0/sum(.vinv)
            } else {
                x$jack_mu[j-1]    <- .mu
                x$jack_blp[,j-1]  <- .blp
                x$jack_vinv[,j-1] <- c(min(.vinv), max(.vinv), sum(diag(.vinv)))
            }
        }
        
        ## add jackknife variance/MSE pieces
        ## "plugin" refers to the first-order plug-in estimate, "so" refers to second-order jackknife
        ##
        ## NOTE: originally we subtracted x$mu and x$tab$blp, but this seemed to result in overcorrection,
        ## so we changed the code to center with respect to jackknife sample means
        if(control$jackknife){
            ## GLS estimator of school fixed effect
            x$var_muhat_plugin <- x$var_muhat
            .mean              <- mean(x$jack_mu)
            x$var_muhat_so     <- ((J-1)/J) * sum( (x$jack_mu - .mean)^2 )
            x$var_muhat        <- x$var_muhat_plugin + x$var_muhat_so

            ## vector of EBLP of school*block means
            x$mse_blp_plugin   <- x$mse_blp
            .mean              <- as.vector(apply(x$jack_blp, 1, mean))
            .tmp <- Reduce("+",lapply(1:J, function(j){ crossprod( matrix(x$jack_blp[,j] - .mean, nrow=1) )}))
            x$mse_blp_so       <- ((J-1)/J) * .tmp
            x$mse_blp          <- x$mse_blp_plugin + x$mse_blp_so
        }

        dsch[[s]]   <- x

        if( !control$quietly && ((s %% 100) == 0) ){
            cat(paste(s,"schools done\n"))
        }
    }

    ## summaries of EBLPS
    if(!control$quietly){
        .blp <- as.vector(unlist(lapply(dsch, function(x){ x$tab$blp })))
        tmp  <- as.vector(unlist(lapply(dsch, function(x){ x$tab$Y_sb })))
        cat(paste("range of block-level EBLPs:",round(min(.blp),digits=4),",",round(max(.blp),digits=4),"\n"))
        cat(paste("cor(Y,EBLP) at block level:",round(as.vector(cor(.blp, tmp)),digits=4),"\n"))
    }
    
    #####################################################################################
    ## use estimated school fixed effects and associated estimated error variances to
    ## estimate variance in individual growth attributable to school fixed effects,
    ## where schools are weighted according to the total number of attached growth
    ## scores, and where we assume that errors for estimated FE by school are
    ## uncorrelated, which is an approximation
    #####################################################################################
    .l    <- sapply(dsch, function(x){ sum(x$tab$nsb) })
    .l    <- .l/sum(.l)
    .y    <- sapply(dsch, function(x){ x$mu })
    .s    <- sapply(dsch, function(x){ x$var_muhat })
    modstats[["estimated_variance_among_schools"]] <- (sum(.l * .y^2) - sum(.l * .s)) - ( (sum(.l * .y))^2 - sum(.l^2 * .s) )
    modstats[["estimated_percvar_among_schools"]]  <- modstats[["estimated_variance_among_schools"]] / modstats[["varY"]]
    
    ## #############################################
    ## EBLP and MSE calculation for target
    ## #############################################
    if(!control$quietly){
        cat("Computing EBLPs and MSE estimates for targets...\n")
    }
    for(s in 1:S){
        x <- dsch[[s]]
        b <- x$oblocks
        n <- x$tab$nsb
        gconfig <- paste(sort(unique(x$tab$grade)), collapse="")
        
        ## proceed if the school has measures for the target blocks
        ## (and if needed, target_contrast blocks)
        validschool <- any(b %in% target_blockids)
        if(!is.null(target_contrast)){
            validschool <- validschool && any(b %in% target_contrast_blockids)
        }
        
        if(validschool){
            ## make lambda
            lambda  <- rep(0.0, x$nblock)
            wh      <- which(b %in% target_blockids)
            ntarget <- sum(n[wh])
            if(target["weights"]=="n"){
                lambda[wh] <- n[wh] / ntarget
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

            ## create matrix to store weights for direct and EBLP estimators.
            ## also stored indicator "obs" which indicates whether there
            ## are observed data for the corresponding block.
            weights             <- matrix(0, ncol=3, nrow=B)
            rownames(weights)   <- .blocknames
            colnames(weights)   <- c("obs","direct","blp")
            weights[b,"obs"]    <- 1
            weights[b,"direct"] <- lambda
            
            ## get direct, EBLP and MSE estimators for target
            est.direct   <- sum(lambda * x$tab$Y_sb)
            mse.direct   <- as.vector(t(lambda) %*% x$R_sb[b,b,drop=F] %*% lambda)
            est.blp      <- sum(lambda * x$tab$blp)
            mse.blp      <- as.vector(t(lambda) %*% x$mse_blp %*% lambda)
            est.hybrid   <- ifelse(mse.blp < mse.direct, est.blp, est.direct)
            mse.hybrid   <- ifelse(mse.blp < mse.direct, mse.blp, mse.direct)
            prmse.direct <- 1 - mse.blp/mse.direct
            tmp <- list(est=c(est.direct = est.direct, mse.direct = mse.direct, est.blp = est.blp, mse.blp = mse.blp, est.hybrid = est.hybrid, mse.hybrid = mse.hybrid, prmse.direct = prmse.direct))
            weights[b,"blp"] <- as.vector( t(lambda) %*% x$Q )
            
            ## package up and save
            x$est         <- data.frame(school = x$school, gconfig = gconfig, ntotal = sum(n), ntarget = ntarget, ncontrast = ncontrast, as.data.frame(as.list(tmp$est)), stringsAsFactors=FALSE)
            x$weights     <- weights
        } else { ## school has insufficient observed data for direct estimator of target
            x$est         <- data.frame(school = x$school, gconfig = gconfig, ntotal = sum(n), ntarget = 0, ncontrast = 0, est.direct = NA, mse.direct = NA, est.blp = NA, mse.blp = NA, est.hybrid = NA, mse.hybrid = NA, prmse.direct = NA, stringsAsFactors=FALSE)
        }
        dsch[[s]] <- x
    }
    
    ## ####################
    ## RETURN (after some clean-up)
    ## ####################
    dblockpairs$iB      <- NULL
    .agg <- do.call("rbind",lapply(dsch, function(x){ x$est }))
    rownames(.agg) <- 1:nrow(.agg)

    if(control$jackknife && !control$return_schjack){
        dsch <- lapply(dsch, function(x){
            x$jack_mu   <- summary(x$jack_mu)
            x$jack_blp  <- t(apply(x$jack_blp, 1, summary))
            x$jack_vinv <- t(apply(x$jack_vinv, 1, summary))
            x
        })
    }
    
    .r <- list(control                = control,
               target                 = target,
               target_contrast        = target_contrast,
               dblock                 = dblock,
               dblockpairs            = dblockpairs,
               dsch                   = dsch,
               bhat_ols               = .bhat,
               modstats               = modstats,
               tab_patterns           = tab_patterns,
               G                      = G[[1]],
               R                      = R,
               aggregated_growth      = .agg)
    
    if(control$return_d){
        .r$d <- d
    }

    if(control$jackknife){
        .r$G_jack     <- G[-1]
        .r$Gstar_jack <- Gstar[-1]
        .r$varhat_G   <- varhat_G
    }

    return(.r)
}
