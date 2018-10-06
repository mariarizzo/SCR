    #######################################################
    ###       Statistical Computing with R              ###
    ###       Maria L. Rizzo                            ###
    ###       Chapman & Hall / CRC                      ###
    ###       ISBN 9781584885450                        ###
    ###                                                 ###
    ###       R code for Chapter 8 Examples             ###
    #######################################################

### Example 8.1 (Permutation distribution of a statistic)

    attach(chickwts)
    x <- sort(as.vector(weight[feed == "soybean"]))
    y <- sort(as.vector(weight[feed == "linseed"]))
    detach(chickwts)

    R <- 999              #number of replicates
    z <- c(x, y)          #pooled sample
    K <- 1:26
    reps <- numeric(R)   #storage for replicates
    t0 <- t.test(x, y)$statistic

    for (i in 1:R) {
        #generate indices k for the first sample
        k <- sample(K, size = 14, replace = FALSE)
        x1 <- z[k]
        y1 <- z[-k]      #complement of x1
        reps[i] <- t.test(x1, y1)$statistic
        }
    p <- mean(c(t0, reps) >= t0)
    p

    hist(reps, main = "", freq = FALSE, xlab = "T (p = 0.202)",
        breaks = "scott")
    points(t0, 0, cex = 1, pch = 16)      #observed T


### Example 8.2 (Permutation distribution of the K-S statistic)    

    # continues Example 8.1
    R <- 999             #number of replicates
    z <- c(x, y)         #pooled sample
    K <- 1:26
    D <- numeric(R)      #storage for replicates
    options(warn = -1)
    D0 <- ks.test(x, y, exact = FALSE)$statistic
    for (i in 1:R) {
        #generate indices k for the first sample
        k <- sample(K, size = 14, replace = FALSE)
        x1 <- z[k]
        y1 <- z[-k]      #complement of x1
        D[i] <- ks.test(x1, y1, exact = FALSE)$statistic
        }
    p <- mean(c(D0, D) >= D0)
    options(warn = 0)
    p

    hist(D, main = "", freq = FALSE, xlab = "D (p = 0.46)",
        breaks = "scott")
    points(D0, 0, cex = 1, pch = 16)      #observed D


### Example 8.3 (Example 8.2, cont.)

    attach(chickwts)
    x <- sort(as.vector(weight[feed == "sunflower"]))
    y <- sort(as.vector(weight[feed == "linseed"]))
    detach(chickwts)
    summary(cbind(x, y))
    options(warn = -1)
    D0 <- ks.test(x, y, exact = FALSE)$statistic
    for (i in 1:R) {
        #generate indices k for the first sample
        k <- sample(K, size = 14, replace = FALSE)
        x1 <- z[k]
        y1 <- z[-k]      #complement of x1
        D[i] <- ks.test(x1, y1, exact = FALSE)$statistic
        }
    p <- mean(c(D0, D) >= D0)
    options(warn = 0)
    p


### Example 8.4 (Finding nearest neighbors)

    library(knnFinder)  #for nn function

    #generate a small multivariate data set
    x <- matrix(rnorm(12), 3, 4)
    y <- matrix(rnorm(12), 3, 4)

    z <- rbind(x, y)
    o <- rep(0, nrow(z))

    DATA <- data.frame(cbind(z, o))
    NN <- nn(DATA, p = nrow(z)-1)

    D <- dist(z)
    round(as.matrix(D), 2)
    
    NN$nn.idx
    round(NN$nn.dist, 2)
    detach(package:knnFinder)


### Example 8.5 (Nearest neighbor statistic)


    library(knnFinder)
    with(chickwts, {
    x <- as.vector(weight[feed == "sunflower"])
    y <- as.vector(weight[feed == "linseed"])})
    z <- c(x, y)
    o <- rep(0, length(z))
    z <- as.data.frame(cbind(z, o))

    NN <- nn(z, p=3)
    cbind(z, NN$nn.idx)

    block1 <- NN$nn.idx[1:12, ]
    block2 <- NN$nn.idx[13:24, ]
    i1 <- sum(block1 < 12.5)
    i2 <- sum(block2 > 12.5)

    c(i1, i2)


### Example 8.6 (Nearest neighbor test)

    library(boot)
    #continues the previous example
    #uses package knnFinder loaded in prev. example
    
    Tn3 <- function(z, ix, sizes) {
        n1 <- sizes[1]
        n2 <- sizes[2]
        n <- n1 + n2
        z <- z[ix, ]
        o <- rep(0, NROW(z))
        z <- as.data.frame(cbind(z, o))
        NN <- nn(z, p=3)
        block1 <- NN$nn.idx[1:n1, ]
        block2 <- NN$nn.idx[(n1+1):n, ]
        i1 <- sum(block1 < n1 + .5)
        i2 <- sum(block2 > n1 + .5)
        return((i1 + i2) / (3 * n))
    }
    N <- c(12, 12)
    
    boot.obj <- boot(data = z, statistic = Tn3,
        sim = "permutation", R = 999, sizes = N)
    boot.obj
    
    tb <- c(boot.obj$t, boot.obj$t0)
    mean(tb >= boot.obj$t0)

    hist(tb, freq=FALSE, main="",
         xlab="replicates of T(n,3) statistic")
    points(boot.obj$t0, 0, cex=1, pch=16)
    detach(package:knnFinder)


### Example 8.7 (Two-sample energy statistic)

    edist.2 <- function(x, ix, sizes) {
        # computes the e-statistic between 2 samples
        # x:          Euclidean distances of pooled sample
        # sizes:      vector of sample sizes
        # ix:         a permutation of row indices of x

        dst <- x
        n1 <- sizes[1]
        n2 <- sizes[2]
        ii <- ix[1:n1]
        jj <- ix[(n1+1):(n1+n2)]
        w <- n1 * n2 / (n1 + n2)

        # permutation applied to rows & cols of dist. matrix
        m11 <- sum(dst[ii, ii]) / (n1 * n1)
        m22 <- sum(dst[jj, jj]) / (n2 * n2)
        m12 <- sum(dst[ii, jj]) / (n1 * n2)
        e <- w * ((m12 + m12) - (m11 + m22))
        return (e)
    }

    d <- 3
    a <- 2 / sqrt(d)
    x <- matrix(rnorm(20 * d), nrow = 20, ncol = d)
    y <- matrix(rnorm(10 * d, a, 1), nrow = 10, ncol = d)
    z <- rbind(x, y)
    dst <- as.matrix(dist(z))

    edist.2(dst, 1:30, sizes = c(20, 10))


### Example 8.8 (Two-sample energy test)

    library(boot)  #for boot function
    dst <- as.matrix(dist(z))
    N <- c(20, 10)

    boot.obj <- boot(data = dst, statistic = edist.2,
        sim = "permutation", R = 999, sizes = N)
    boot.obj
    
    #calculate the ASL    
    e <- boot.obj$t0
    tb <- c(e, boot.obj$t)
    mean(tb >= e)

    hist(tb, main = "", breaks="scott", freq=FALSE,
        xlab="Replicates of e")
    points(e, 0, cex=1, pch=16)


    #energy test applied under F=G
    d <- 3
    a <- 0
    x <- matrix(rnorm(20 * d), nrow = 20, ncol = d)
    y <- matrix(rnorm(10 * d, a, 1), nrow = 10, ncol = d)
    z <- rbind(x, y)
    dst <- as.matrix(dist(z))

    N <- c(20, 10)
    dst <- as.matrix(dist(z))
    boot.obj <- boot(data = dst, statistic = edist.2,
        sim="permutation", R=999, sizes=N)
    boot.obj

    #calculate the ASL    
    e <- boot.obj$t0
    E <- c(boot.obj$t, e)
    mean(E >= e)

    hist(E, main = "", breaks="scott",
        xlab="Replicates of e", freq=FALSE)
    points(e, 0, cex=1, pch=16)

### Example 8.9 (k-sample energy distances)

    library(energy)  #for edist
    z <- iris[ , 1:4]
    dst <- dist(z)

    edist(dst, sizes = c(50, 50, 50), distance = TRUE)
    detach(package:energy)


### Example 8.11 (Distance covariance statistic)

    dCov <- function(x, y) {
        x <- as.matrix(x)
        y <- as.matrix(y)
        n <- nrow(x)
        m <- nrow(y)
        if (n != m || n < 2) stop("Sample sizes must agree")
        if (! (all(is.finite(c(x, y)))))
            stop("Data contains missing or infinite values")

        Akl <- function(x) {
            d <- as.matrix(dist(x))
            m <- rowMeans(d)
            M <- mean(d)
            a <- sweep(d, 1, m)
            b <- sweep(a, 2, m)
            return(b + M)
        }
        A <- Akl(x)
        B <- Akl(y)
        dCov <- sqrt(mean(A * B))
        dCov
    }

    z <- as.matrix(iris[1:50, 1:4])
    x <- z[ , 1:2]
    y <- z[ , 3:4]
    # compute the observed statistic
    dCov(x, y)

### Example 8.12 (Distance correlation statistic)

    DCOR <- function(x, y) {
        x <- as.matrix(x)
        y <- as.matrix(y)
        n <- nrow(x)
        m <- nrow(y)
        if (n != m || n < 2) stop("Sample sizes must agree")
        if (! (all(is.finite(c(x, y)))))
            stop("Data contains missing or infinite values")
        Akl <- function(x) {
            d <- as.matrix(dist(x))
            m <- rowMeans(d)
            M <- mean(d)
            a <- sweep(d, 1, m)
            b <- sweep(a, 2, m)
            return(b + M)
        }
        A <- Akl(x)
        B <- Akl(y)
        dCov <- sqrt(mean(A * B))
        dVarX <- sqrt(mean(A * A))
        dVarY <- sqrt(mean(B * B))
        dCor <- sqrt(dCov / sqrt(dVarX * dVarY))
        list(dCov=dCov, dCor=dCor, dVarX=dVarX, dVarY=dVarY)
    }

    z <- as.matrix(iris[1:50, 1:4])
    x <- z[ , 1:2]
    y <- z[ , 3:4]

    DCOR(x, y)
    

### Example 8.13 (Distance covariance test)

    ndCov2 <- function(z, ix, dims) {
        p <- dims[1]
        q1 <- p + 1
        d <- p + dims[2]
        x <- z[ , 1:p]     #leave x as is
        y <- z[ix, q1:d]   #permute rows of y
        return(nrow(z) * dCov(x, y)^2)
    }

    library(boot)
    z <- as.matrix(iris[1:50, 1:4])
    boot.obj <- boot(data = z, statistic = ndCov2, R = 999,
        sim = "permutation", dims = c(2, 2))

    tb <- c(boot.obj$t0, boot.obj$t)
    hist(tb, nclass="scott", xlab="", main="",
            freq=FALSE)
    points(boot.obj$t0, 0, cex=1, pch=16)

    mean(tb >= boot.obj$t0)
    boot.obj
