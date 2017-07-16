# Extra parameters to HOPACH as alpha1, alpha2, a, pow, and ... (for user fn)
# Options for HOPACH distances include:
# USER, BINARY, S-FUN (has optional extra parameters for S-function),
# SEPARATELY - JACCARD, BINARY 3

#----------------------------------------------------------------------------#
#          DISTANCES: user, binary, jaccard (alpha1 = 1, alpha2 = 0),        #    
#                              binary3, S-function                           #
#----------------------------------------------------------------------------#

# user.distance <- function(x, y, ...) {}

# General binary
binary.distance <- function(x, y, alpha1, alpha2) {
    n <- length(x)
    ones <- which(x == 1)
    a <- sum(y[ones])
    d <- n - length(ones) - sum(y[-ones])
    bc <- n - a - d
    bc/(alpha1*a + alpha2*d + bc)
}

# Special case binary with alpha1 = 1, alpha2 = 0
jaccard.distance <- function(x, y) {
    non.zero <- sum(x + y != 0)
    diff <- sum(x + y == 1)
    return(diff/non.zero)
}

# d(1, 1) = 0; d(0, 0) = 1/2; d(1, 0) = d(0, 1) = 1,
# where 1 indicates significant and 0 - not significant
b3 <- function(x, y) {
    s <- rowSums(cbind(x, y))
    s[!s] <- 3/2
    sum(2 - s)
}

# S-function
# x, y are p-values
sfun.distance <- function(x, y, a, pow) { 
    sfun <- function(x) 1 - exp(-a * x^pow)          # shape of the curve
    signif <- which((x>0.2) + (y>0.2) < 2)
    newx <- x[signif]
    newy <- y[signif]
    sum(abs(sfun(newx) - sfun(newy)))/length(signif)
}

#----------------------------------------------------------------------------#
#                        distance matrix functions                           #
#----------------------------------------------------------------------------#

binary.distance.matrix <- function(X, alpha1 = 1, alpha2 = 0) {
    n <- nrow(X)
    dist.mat <- matrix(NA, n, n)
    for (i in 2:n) {
        for (j in 1:(i-1)) {
            dist.mat[i, j] <- dist.mat[j, i] <- binary.distance(X[i, ], X[j, ], 
                                                                alpha1, alpha2)
        }
    }
    diag(dist.mat) <- 0
    return(dist.mat)
}

jaccard.distance.matrix <- function(X) {
    n <- nrow(X)
    dist.mat <- matrix(NA, n, n)
    for (i in 2:n) {
        for (j in 1:(i-1)) {
            dist.mat[i, j] <- dist.mat[j, i] <- jaccard.distance(X[i, ], X[j, ])
        }
    }
    diag(dist.mat) <- 0
    return(dist.mat)
}

binary3.distance.matrix <- function(X) {
    n <- nrow(X)
    dist.mat <- matrix(NA, n, n)
    for (i in 2:n) {
        for (j in 1:(i-1)) {
            dist.mat[i, j] <- dist.mat[j, i] <- binary3.distance(X[i, ], X[j, ])
        }
    }
    diag(dist.mat) <- 0
    return(dist.mat)
}


sfun.distance.matrix <- function(X, a = 150, pow = 2) {
    n <- nrow(X)
    dist.mat <- matrix(NA, n, n)
    for (i in 2:n) {
        for (j in 1:(i-1)) {
            dist.mat[i, j] <- dist.mat[j, i] <- sfun.distance(X[i, ], X[j, ], 
                                                              a, pow)
        }
    }
    diag(dist.mat) <- 0
    return(dist.mat)
}

user.distance.matrix <- function(X, ...) {
    n <- nrow(X)
    dist.mat <- matrix(NA, n, n)
    for (i in 2:n) {
        for (j in 1:(i-1)) {
            dist.mat[i, j] <- dist.mat[j, i] <- user.distance(X[i, ], X[j, ], ...)
        }
    }
    diag(dist.mat) <- 0
    return(dist.mat)
}

##############################################################################
####################    updated distancematrix function    ###################
##############################################################################

distancematrix <- function(X, d, alpha1 = 1, alpha2 = 0, a = 150, pow = 2, 
                           na.rm = TRUE, ...) {
    X <- as.matrix(X)
    if (d == "euclid") 
        return(disseuclid(X, na.rm))
    if (d == "cor") 
        return(disscor(X, na.rm))
    if (d == "abscor") 
        return(dissabscor(X, na.rm))
    if (d == "cosangle") 
        return(disscosangle(X, na.rm))
    if (d == "abscosangle") 
        return(dissabscosangle(X, na.rm))
    if (d == "binary")
        return(dissbinary(X, alpha1, alpha2, na.rm))    
    if (d == "jaccard")
        return(dissjaccard(X, na.rm))    
    if (d == "binary3")
        return(dissbinary3(X, na.rm))    
    if (d == "sfun") # S shaped function
        return(disssfun(X, a, pow, na.rm))   
    if (d == "user")
        return(dissuser(X, na.rm, ...)) # user-specified distance
    stop("Distance metric ", d, " not available")
}

#----------------------------------------------------------------------------#
#  auxillary dissjaccard, dissbinary3, dissbinary, and disssfun functions    #
#----------------------------------------------------------------------------#

dissbinary <- function(X, alpha1, alpha2, na.rm = TRUE) {
    if (!is.matrix(X)) {
        stop(paste(sQuote("X"), "not a matrix"))
    }
    out <- as.dist(binary.distance.matrix(X, alpha1, alpha2))
    dmat <- new("hdist", Data = out[1:length(out)], Size = attr(out, 
                                                                "Size"), Labels = (1:(attr(out, "Size"))), Call = as.character(attr(out, 
                                                                                                                                    "call")[3]))
    return(dmat)
}

dissjaccard <- function(X, na.rm = TRUE) {
    if (!is.matrix(X)) {
        stop(paste(sQuote("X"), "not a matrix"))
    }
    out <- dist(X, method = "binary")
    dmat <- new("hdist", Data = out[1:length(out)], Size = attr(out, 
                                                                "Size"), Labels = (1:(attr(out, "Size"))), Call = as.character(attr(out, 
                                                                                                                                    "call")[3]))
    return(dmat)
}

dissbinary3 <- function(X, na.rm = TRUE) {
    if (!is.matrix(X)) {
        stop(paste(sQuote("X"), "not a matrix"))
    }
    out <- as.dist(binary3.distance.matrix(X))
    dmat <- new("hdist", Data = out[1:length(out)], Size = attr(out, 
                                                                "Size"), Labels = (1:(attr(out, "Size"))), Call = as.character(attr(out, 
                                                                                                                                    "call")[3]))
    return(dmat)
}

disssfun <- function(X, a, pow, na.rm = TRUE) {
    if (!is.matrix(X)) {
        stop(paste(sQuote("X"), "not a matrix"))
    }
    out <- as.dist(sfun.distance.matrix(X, a, pow))
    dmat <- new("hdist", Data = out[1:length(out)], Size = attr(out, 
                                                                "Size"), Labels = (1:(attr(out, "Size"))), Call = as.character(attr(out, 
                                                                                                                                    "call")[3]))
    return(dmat)
}

dissuser <- function(X, na.rm = TRUE, ...) {
    if (!is.matrix(X)) {
        stop(paste(sQuote("X"), "not a matrix"))
    }
    out <- as.dist(user.distance.matrix(X, ...))
    dmat <- new("hdist", Data = out[1:length(out)], Size = attr(out, 
                                                                "Size"), Labels = (1:(attr(out, "Size"))), Call = as.character(attr(out, 
                                                                                                                                    "call")[3]))
    return(dmat)
}

##########################################################################
#####              updated distancevector function                  ######
##########################################################################

distancevector <- function (X, y, d, alpha1 = 1, alpha2 = 0, a = 150, pow = 2, 
                            na.rm = TRUE, ...) {
    X <- as.matrix(X)
    y <- as.vector(y)
    if (d == "cosangle") 
        return(vdisscosangle(X, y, na.rm))
    if (d == "abscosangle") 
        return(vdissabscosangle(X, y, na.rm))
    if (d == "euclid") 
        return(vdisseuclid(X, y, na.rm))
    if (d == "abseuclid") 
        return(vdissabseuclid(X, y, na.rm))
    if (d == "cor") 
        return(vdisscor(X, y, na.rm))
    if (d == "abscor") 
        return(vdissabscor(X, y, na.rm))
    if (d == "binary") 
        return(vdissbinary(X, y, alpha1, alpha2)) 
    if (d == "jaccard")
        return(vdissjaccard(X, y))
    if (d == "binary3") 
        return(vdissbinary3(X, y)) 
    if (d == "sfun") 
        return(vdisssfun(X, y, a, pow))
    if (d == "user")
        return(vdissuser(X, y, ...))            
    stop("Distance metric ", d, " not available")
}
    
#----------------------------------------------------------------------------#
#      auxillary function for distancevector - vdissbinary and others        #
#----------------------------------------------------------------------------#

vdissbinary <- function(X, y, alpha1, alpha2) {
    if (!is.matrix(X)) 
        stop("First arg to vdissbinary() must be a matrix")
    if (!is.vector(y)) 
        stop("Second arg to vdissbinary() must be a vector")
    if (length(y) != length(X[1, ])) 
        stop("matrix and vector dimensions do not agree in vdissbinary()") 
    out <- apply(X, 1, binary.distance, y = y, alpha1 = alpha1, alpha2 = alpha2)
    return(out)
}

vdissjaccard <- function(X, y) {
    if (!is.matrix(X)) 
        stop("First arg to vdissjaccard() must be a matrix")
    if (!is.vector(y)) 
        stop("Second arg to vdissjaccard() must be a vector")
    if (length(y) != length(X[1, ])) 
        stop("matrix and vector dimensions do not agree in vdissjaccard()")
    out <- apply(X, 1, jaccard.distance, y = y)
    return(out)
}

vdissbinary3 <- function(X, y) {
    if (!is.matrix(X)) 
        stop("First arg to vdissbinary3() must be a matrix")
    if (!is.vector(y)) 
        stop("Second arg to vdissbinary3() must be a vector")
    if (length(y) != length(X[1, ])) 
        stop("matrix and vector dimensions do not agree in vdissbinary3()") 
    out <- apply(X, 1, binary3.distance, y = y)
    return(out)
}

vdisssfun <- function(X, y, a, pow) {
    if (!is.matrix(X)) 
        stop("First arg to vdisssfun() must be a matrix")
    if (!is.vector(y)) 
        stop("Second arg to vdisssfun() must be a vector")
    if (length(y) != length(X[1, ])) 
        stop("matrix and vector dimensions do not agree in vdisssfun()") 
    out <- apply(X, 1, sfun.distance, y = y, a = a, pow = pow)
    return(out)
}

vdissuser <- function(X, y, ...) {
    if (!is.matrix(X)) 
        stop("First arg to vdissuser() must be a matrix")
    if (!is.vector(y)) 
        stop("Second arg to vdissuser() must be a vector")
    if (length(y) != length(X[1, ])) 
        stop("matrix and vector dimensions do not agree in vdissuser()") 
    out <- apply(X, 1, user.distance, y = y, ...)
    return(out)
}

##############################################################################
##########            updated parameters for functions              ##########
##############################################################################

hopach <- function (data, dmat = NULL, d = "cosangle", clusters = "best", 
    K = 15, 
    kmax = 9, 
    khigh = 9, 
    coll = "seq", 
    newmed = "medsil", 
    mss = "med", 
    impr = 0, 
    initord = "co", 
    ord = "own", 
    verbose = FALSE, alpha1 = 1, alpha2 = 0, a = 150, pow = 2, ...) 
{
    if (inherits(data, "ExpressionSet")) 
        data <- exprs(data)
    data <- as.matrix(data)
    if (is.null(dmat)) {
        dmat <- distancematrix(data, d, alpha1, alpha2, a, pow, na.rm = TRUE, ...)
    }
    else if (is.matrix(dmat)) {
        dmat <- as(dmat, "hdist")
    }
    else if (class(dmat) == "dist") {
        dmat <- hdist(Data = as.numeric(dmat), Size = attr(dmat, 
            "Size"), Labels = (1:(attr(dmat, "Size"))), Call = as.character(attr(dmat, 
            "call"))[3])
    }
    else if (!is.hdist(dmat)) {
        stop("Distance matrix could not be created into hdist object.")
    }
    if (K > 15) {
        K <- 15
        warning("K set to 15 - can't do more than 15 splits")
    }
    if (K < 1) {
        K <- 1
        warning("K set to 1 - can't do less than 1 level")
    }
    if (clusters != "none") {
        cuttree <- mssrundown(data, K, kmax, khigh, d, dmat, 
            initord, coll, newmed, stop = (clusters == "greedy"), 
            finish = TRUE, within = mss, between = mss, impr, 
            verbose, alpha1, alpha2, a, pow, ...)
        if (cuttree[[1]] > 1) 
            cutord <- orderelements(cuttree, data, rel = ord, 
                d, dmat, alpha1, alpha2, a, pow, ...)[[2]]
        else cutord <- NULL
        out1 <- list(k = cuttree[[1]], medoids = cuttree[[2]], 
            sizes = cuttree[[3]], labels = cuttree[[4]], order = cutord)
        finaltree <- msscomplete(cuttree, data, K, khigh, d, 
            dmat, within = mss, between = mss, verbose, alpha1, alpha2, 
            a, pow, ...)
    }
    else {
        out1 <- NULL
        finaltree <- msscomplete(mssinitlevel(as.matrix(data), 
            kmax, khigh, d, dmat, within = mss, between = mss, 
            initord, verbose = FALSE, alpha1, alpha2, a, pow, ...), 
            data, K, khigh, d, dmat, within = mss, between = mss, 
            verbose, alpha1, alpha2, a, pow, ...)
    }
    dimnames(finaltree[[6]]) <- list(NULL, c("label", "medoid"))
    out2 <- list(labels = finaltree[[4]], order = orderelements(finaltree, 
        data, rel = ord, d, dmat, alpha1, alpha2, a, pow, ...)[[2]], 
        medoids = finaltree[[6]])
    return(list(clustering = out1, final = out2, call = match.call(), 
        metric = d))
}

mssrundown <- function (data, K = 16, kmax = 9, khigh = 9, d = "cosangle", 
    dmat = NULL, initord = "co", coll = "seq", newmed = "medsil", 
    stop = TRUE, finish = FALSE, within = "med", between = "med", 
    impr = 0, verbose = FALSE, alpha1, alpha2, a, pow, ...) 
{
    if (!is.matrix(data)) 
        stop("First arg to mssrundown() must be a matrix")
    bestlevel <- level <- mssinitlevel(data, kmax, khigh, d, 
        dmat, within, between, initord, verbose, alpha1, alpha2, a, pow, ...)
    bestmss <- mss <- labelstomss(level[[4]], dmat, khigh, within, 
        between)
    bestl <- l <- 1
    ind <- 0
    if (verbose) 
        cat("Searching for main clusters... \n")
    if (level[[5]] == 1) 
        return(level)
    while ((l <= K) && (ind == 0)) {
        if (verbose) 
            cat("Level ", l, "\n")
        if (coll == "seq") 
            levelc <- msscollap(data, level, khigh, d, dmat, 
                newmed, within, between, impr, alpha1, alpha2, a, pow, ...)
        if (coll == "all") 
            levelc <- mssmulticollap(data, level, khigh, d, dmat, 
                newmed, within, between, impr, alpha1, alpha2, a, pow, ...)
        mss <- labelstomss(levelc[[4]], dmat, khigh, within, 
            between)
        if (mss >= bestmss & stop == TRUE) 
            ind <- 1
        else {
            if (mss < bestmss) {
                bestlevel <- levelc
                bestmss <- mss
                bestl <- l
            }
        }
        l <- l + 1
        if (l <= K) {
            level <- mssnextlevel(data, levelc, dmat, kmax, khigh, 
                within, between)
            if (finish == TRUE) {
                if (sum(trunc(level[[4]]/10) * 10 == level[[4]]) == 
                  length(level[[4]]) & l <= K) {
                  ind <- 1
                  bestlevel <- levelc
                  bestmss <- mss
                  bestl <- (l - 1)
                }
            }
        }
    }
    if (verbose) 
        cat("Identified", bestlevel[[1]], " main clusters in level", 
            bestl, "with mSS =", bestmss, "\n")
    return(bestlevel)
}

orderelements <- function (level, data, rel = "own", d = "cosangle", 
    dmat = NULL, alpha1, alpha2, a, pow, ...) 
{
    idn <- 1:length(data[, 1])
    k <- level[[1]]
    labels <- level[[4]]
    medoids <- level[[2]]
    clussizes <- level[[3]]
    ord <- order(labels)
    idnord <- idn[ord]
    subdataord <- data[ord, ]
    if (is.null(dmat)) 
        dmat <- distancematrix(data, d, alpha1, alpha2, a, pow, na.rm = TRUE, ...)
    distord <- dmat[ord, ]
    labelsord <- labels[ord]
    count <- 1
    for (j in (1:k)) {
        start <- count
        end <- count + clussizes[j] - 1
        if (clussizes[j] > 2) {
            tempid <- idnord[start:end]
            if (rel == "co") {
                distj <- distord[, ord][start:end, start:end]
                idnord[start:end] <- tempid[improveordering(distj)]
            }
            else {
                if (rel == "neighbor") {
                  if (j < k) 
                    mednext <- medoids[j + 1]
                  else mednext <- medoids[j - 1]
                }
                else mednext <- medoids[j]
                dmednext <- distord[start:end, mednext]
                if (rel == "neighbor") {
                  if (j < k) 
                    ordtemp <- rev(order(dmednext))
                  else ordtemp <- order(dmednext)
                }
                else ordtemp <- order(dmednext)
                idnord[start:end] <- tempid[ordtemp]
            }
        }
        else idnord[start:end] <- idnord[start:end]
        count <- count + clussizes[j]
    }
    list(data[idnord, ], idnord)
}

msscomplete <- function (level, data, K = 16, khigh = 9, d = "cosangle", dmat = NULL, 
    within = "med", between = "med", verbose = FALSE, alpha1, alpha2, a, pow, ...) 
{                             # don't really need d and extra parameters for d here
    if (!is.matrix(data)) 
        stop("First arg to msscomplete() must be a matrix")
    count <- digits(level[[4]][1])
    if (verbose) 
        cat("Running down without collapsing from Level", count, 
            "\n")
    while ((max(level[[3]]) > 3) & (count < K)) {
        level <- newnextlevel(data, level, dmat, 2, khigh)
        count <- count + 1
        if (verbose) 
            cat("Level", count, "\n")
    }
    return(level)
}

mssinitlevel <- function (data, kmax = 9, khigh = 9, d = "cosangle", dmat = NULL, 
    within = "med", between = "med", ord = "co", verbose = FALSE, alpha1, alpha2, 
    a, pow, ...) 
{
    if (!is.matrix(data)) 
        stop("First arg to mssinitlevel() must be a matrix")
    p <- length(data[, 1])
    
#-------------------------- THE PROBLEM --------------------------------------#
# dmat has been converted to an "hdist" object already, so it is NOT a matrix #
# no need to calculate it again and use distancematrix() here                 #
#-----------------------------------------------------------------------------#

#   if (is.matrix(dmat)) 
#       dmat <- as.hdist(dmat)
#   else dmat <- distancematrix(data, d = d)
    if (dmat@Size != p) 
        stop("Data and distance matrix dimensions do not agree in mssinitlevel()")
    m <- msscheck(dmat, kmax, khigh, within, between)
    if (m[1] == 1) {
        if (verbose) 
            cat("No strong evidence for clusters in the first level - \n continuing to split root node anyway. \n")
        m <- msscheck(dmat, kmax, khigh, within, between, force = TRUE)
    }
    pamobj <- pam(dmat@Data, m[1], diss = TRUE)
    rowmedoids <- pamobj$medoids
    final <- ifelse(max(pamobj$clusinfo[, 1]) < 3, 1, 0)
    if (m[1] > 2) {
        medoidsdata <- as.matrix(data[rowmedoids, ])
        l <- length(rowmedoids)
        medoidsdist <- dmat[rowmedoids, rowmedoids]
        if (ord == "co") 
            medoidsord <- improveordering(medoidsdist)
        if (ord == "clust") {
            mpamobj <- pam(medoidsdist@Size, 2, diss = TRUE)
            labelsmed <- mpamobj$clustering
            medmed <- mpamobj$medoids
            clussizes <- mpamobj$clusinfo[, 1]
            prevlevel <- mssnextlevel(medoidsdata, list(2, medmed, 
                clussizes, labelsmed, 0, cbind(c(1, 2), medmed)), 
                dmat = medoidsdist, kmax, khigh, within, between)
            final <- prevlevel[[5]]
            if (final == 0) {
                depth <- (l - 1)
                for (j in (1:depth)) {
                  if (final == 0) {
                    clustnext <- mssnextlevel(medoidsdata, prevlevel, 
                      dmat = medoidsdist, kmax, khigh, within, 
                      between)
                    final <- clustnext[[5]]
                  }
                  if (final == 1) {
                    j <- depth
                    prevlevel <- clustnext
                  }
                }
            }
            medoidsord <- orderelements(prevlevel, medoidsdata, rel = "neighbor", 
                d = d, dmat = medoidsdist, alpha1, alpha2, a, pow, ...)[[2]]
        }
        k <- m[1]
        rowmedoids <- rowmedoids[medoidsord]
        labels <- lab2 <- pamobj$clustering
        for (j in (1:k)) lab2[labels == medoidsord[j]] <- j
        output <- list(k, rowmedoids, pamobj$clusinfo[, 1][medoidsord], 
            lab2, final, cbind(1:k, rowmedoids))
    }
    else output <- list(2, pamobj$medoids, pamobj$clusinfo[, 
        1], pamobj$clustering, final, cbind(1:2, pamobj$medoids))
    return(output)
}

msscollap <- function (data, level, khigh, d = "cosangle", dmat = NULL, newmed = "medsil", 
    within = "med", between = "med", impr = 0, alpha1, alpha2, a, pow, ...) 
{
    newk <- level[[1]]
    mss1 <- labelstomss(level[[4]], dmat, khigh, within, between)
    maxncoll <- max(0, newk - 2)
    ncoll <- 0
    coll <- 1
    if (newk <= 2) 
        coll <- 0
    while ((coll == 1) && (ncoll <= maxncoll)) {
        levelc <- collap(data, level, d, dmat, newmed, alpha1, alpha2, a, pow, ...)
        mss2 <- labelstomss(levelc[[4]], dmat, khigh, within, 
            between)
        if (mss1 == 0) 
            r <- 0
        else r <- (mss1 - mss2)/mss1
        if (r < impr) 
            coll <- 0
        else {
            mss1 <- mss2
            level <- levelc
            ncoll <- ncoll + 1
        }
    }
    return(level)
}

collap <- function (data, level, d = "cosangle", dmat = NULL, newmed = "medsil", 
    alpha1, alpha2, a, pow, ...) 
{
    k <- level[[1]]
    if (k < 3) {
        warning("Not enough medoids to use newmed='medsil' in collap() - \n using newmed='nn' instead \n")
        newmed <- "nn"
    }
    medoids <- level[[2]]
    clussizes <- level[[3]]
    if (sum(is.na(clussizes))) 
        warning("NA in clussizes")
    medoidsdata <- data[medoids, ]
    if (sum(is.na(medoidsdata)) > 0) 
        warning("missing value(s) in medoidsdata in collap()")
    distmed <- dmat[medoids, medoids]
    distv <- distmed@Data
    indexmin <- order(distv)[1]
    best <- vectmatrix(indexmin, k)
    clustfinal <- paircoll(best[1], best[2], data, level, d, 
        dmat, newmed, alpha1, alpha2, a, pow, ...)
    return(clustfinal)
}

paircoll <- function (i, j, data, level, d = "cosangle", dmat = NULL, 
    newmed = "medsil", alpha1, alpha2, a, pow, ...) 
{
    p <- length(data[, 1])
    k <- level[[1]]
    labels <- level[[4]]
    medoids <- level[[2]]
    clussizes <- level[[3]]
    N <- length(level[[6]][, 1])
    block <- level[[6]][(N - k + 1):N, ]
    if (N == k) 
        prevblock <- NULL
    else prevblock <- level[[6]][1:(N - k), ]
    if (max(i, j) > k) 
        stop("The clusters to collapse do not exist in paircoll()")
    labeli <- labels[medoids[i]]
    labelj <- labels[medoids[j]]
    oldlabels <- labels
    labels[labels == labelj] <- labeli
    trunclabels <- trunc(oldlabels/10)
    labelparents <- unique(trunclabels)
    parenti <- order(labelparents)[labelparents == trunc(labeli/10)]
    parentj <- order(labelparents)[labelparents == trunc(labelj/10)]
    if (newmed == "nn") 
        fakemed <- (data[medoids[i], ] * clussizes[i] + data[medoids[j], 
            ] * clussizes[j])/(clussizes[i] + clussizes[j])
    if (newmed == "uwnn") 
        fakemed <- (data[medoids[i], ] + data[medoids[j], ])/2
    if (newmed == "nn" || newmed == "uwnn") {
        rowsub <- (1:p)[labels == labeli]
        distsfm <- distancevector(data[rowsub, ], as.vector(fakemed), 
            d, alpha1, alpha2, a, pow, na.rm = TRUE, ...)
        medoids[i] <- rowsub[order(distsfm)[1]]
    }
    else {
        colldist <- as(dmat[labels == labeli, labels == labeli], 
            "matrix")
        rowsub <- (1:p)[labels == labeli]
        if (newmed == "center") {
            sumdist <- rowSums(colldist)
        }
        if (newmed == "medsil") {
            othermed <- medoids[-c(i, j)]
            collp <- length(labels[labels == labeli])
            othern <- length(othermed)
            if (othern == 0) 
                stop("Not enough medoids to use newmed='medsil' in paircoll()")
            if (is(dmat, "hdist")) {
                if (othern == 1) {
                  otherdist <- dmat[labels == labeli, othermed]
                }
                else {
                  otherdist <- dmat[labels == labeli, othermed]
                }
            }
            if (othern == 1) 
                b <- otherdist
            else b <- apply(otherdist, 1, min)
            b <- matrix(b, nrow = collp, ncol = collp)
            diag(b) <- 0
            b <- abs(b - colldist)/pmax(colldist, b)
            sumdist <- rowSums(b)
            medoids[i] <- rowsub[order(sumdist, decreasing = TRUE) == 
                1]
        }
    }
    k <- k - 1
    clussizes[i] <- clussizes[i] + clussizes[j]
    block[i, 2] <- medoids[i]
    if (j <= k) {
        for (l in (j:k)) {
            medoids[l] <- medoids[l + 1]
            clussizes[l] <- clussizes[l + 1]
            block[l, ] <- block[l + 1, ]
        }
    }
    medoids <- medoids[1:k]
    clussizes <- clussizes[1:k]
    block <- block[1:k, ]
    if (labels[labels == labeli][1]/10 == trunclabels[labels == 
        labeli][1]) {
        labels[labels == labeli] <- labels[labels == labeli] + 
            1
        labels[labels == labelj] <- labels[labels == labelj] + 
            1
        block[i, 1] <- block[i, 1] + 1
    }
    return(list(k, medoids, clussizes, labels, level[[5]], rbind(prevblock, 
        block)))
}

mssmulticollap <- function (data, level, khigh, d = "cosangle", dmat = NULL, newmed = "medsil", 
    within = "med", between = "med", impr = 0, alpha1, alpha2, a, pow, ...) 
{
    if (!is.matrix(data)) 
        stop("First arg to mssmulticollap() must be a matrix")
    medoids <- level[[2]]
    medoidsdata <- data[medoids, ]
    if (sum(is.na(medoidsdata)) > 0) 
        warning("missing value(s) in medoidsdata in mssmulticollap()")
    distmed <- dmat[medoids, medoids]
    k <- level[[1]]
    ord <- order(distmed@Data)
    mss1 <- labelstomss(level[[4]], dmat, khigh, within, between)
    maxncoll <- max(0, k * (k - 1)/2)
    ncoll <- 0
    i <- 1
    while (i <= maxncoll) {
        clusts <- vectmatrix(ord[i], k)
        if (k < 3) {
            if (newmed == "medsil") 
                warning("Can't use newmed=medsil with less than 3 clusters. \n Substituting newmed=nn")
            newmed <- "nn"
        }
        levelc <- paircoll(clusts[1], clusts[2], data, level, 
            d, dmat, newmed, alpha1, alpha2, a, pow, ...)
        mss2 <- labelstomss(levelc[[4]], dmat, khigh, within, 
            between)
        if (mss1 == 0) 
            r <- 0
        else r <- (mss1 - mss2)/mss1
        if (r >= impr) {
            mss1 <- mss2
            level <- levelc
            ncoll <- ncoll + 1
            k <- level[[1]]
            maxncoll <- max(0, k * (k - 1)/2)
            i <- 0
            medoids <- level[[2]]
            medoidsdata <- data[medoids, ]
            if (sum(is.na(medoidsdata)) > 0) 
                warning("missing value(s) in medoidsdata in mssmulticollap()")
            distmed <- dmat[medoids, medoids]
            ord <- order(distmed@Data)
        }
        i <- i + 1
    }
    return(level)
}

