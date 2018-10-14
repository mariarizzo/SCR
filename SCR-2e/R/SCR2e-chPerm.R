#######################################################
###       Statistical Computing with R  2e          ###
###       Maria L. Rizzo                            ###
###       Chapman & Hall/CRC The R Series           ###
###       ISBN 9781466553323 - CAT# K15269          ###
###       January 2019                              ###
###                                                 ###
###       R code for Chapter 10                     ###
###       Permutation Tests                         ###
#######################################################

# packages to install: energy, yaImpute
    
### Example 10.1 (Permutation distribution of a statistic)

    attach(chickwts)
    x <- sort(weight[feed == "soybean"])
    y <- sort(weight[feed == "linseed"])
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


### Example 10.2 (Permutation distribution of the K-S statistic)    

    # continues Example 10.1
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


### Example 10.3 (Example 10.2, cont.)

    attach(chickwts)
    x <- sort(weight[feed == "sunflower"])
    y <- sort(weight[feed == "linseed"])
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


### Example 10.4 (Finding nearest neighbors)
### using yaImpute::ann
    
    set.seed(439)  
    library(yaImpute)  #for ann function
    
    #generate a small multivariate data set
    x <- matrix(rnorm(12), 3, 4)
    y <- matrix(rnorm(12), 3, 4)
    z <- rbind(x, y)
    k <- nrow(z)  #number of nearest neighbors desired
    
    ## Do an exact kd-tree search
    kd.exact <- ann(ref=z, target=z,
                    tree.type="kd", k=k, verbose=FALSE)
    kd.exact$knnIndexDist[,1:k]           #NN indices
    round(sqrt(kd.exact$knnIndexDist[,-(1:k)]),2)  #Euclidean distances
    
    ## Do an approximate kd-tree search
    kd.approx <- ann(ref=z, target=z,
                     tree.type="kd", k=k, eps=100, verbose=FALSE)
    kd.approx$knnIndexDist[,1:k]           #NN indices
    detach(package:yaImpute)
    
    
    ### Example 10.5 (Nearest neighbor statistic)
    ### using yaImpute::ann
    
    library(yaImpute)
    attach(chickwts)
    x <- weight[feed == "sunflower"]
    y <- weight[feed == "linseed"]
    z <- as.matrix(c(x, y))
    detach(chickwts)
    
    k <- 4  #want first 3 nearest neighbors
    NN <- ann(ref=z, target=z, tree.type="kd", k=k, verbose=FALSE)
    idx <- NN$knnIndexDist[,1:k]
    nn.idx <- idx[,-1]   #first NN is in column 2
    
    block1 <- nn.idx[1:12, ]
    block2 <- nn.idx[13:24, ]
    i1 <- sum(block1 < 12.5)
    i2 <- sum(block2 > 12.5)
    
    c(i1, i2)
    
    detach(package:yaImpute)
    
    
### Example 10.6 (Nearest neighbor test)
### using yaImpute::ann
    
    library(boot)
    #continues the previous example
    
    ## function to return the matrix of indices NN_j of nearest neighbors
    NN.idx <- function(x, tree.type="kd", k=NROW(x)) {
      x <- as.matrix(x)
      k <- min(c(k+1, NROW(x)))
      NN <- yaImpute::ann(ref=x, target=x, 
                          tree.type="kd", k=k, verbose=FALSE)
      idx <- NN$knnIndexDist[,1:k]
      nn.idx <- idx[,-1]   #first NN is in column 2
      row.names(nn.idx) <- idx[,1]
      nn.idx
    }
    
    ## function to compute the NN statistic T(n,3)
    Tn3 <- function(z, ix=1:NROW(z), sizes) {
      z <- as.matrix(z)
      n1 <- sizes[1]
      n2 <- sizes[2]
      n <- n1 + n2
      z <- as.matrix(z[ix, ])
      nn.idx <- NN.idx(z, k=3)
      block1 <- nn.idx[1:n1, ]
      block2 <- nn.idx[(n1+1):n, ]
      i1 <- sum(block1 < n1 + .5)
      i2 <- sum(block2 > n1 + .5)
      return((i1 + i2) / (3 * n))
    }
    
    attach(chickwts)
    x <- weight[feed == "sunflower"]
    y <- weight[feed == "linseed"]
    z <- c(x, y)
    detach(chickwts)
    
    N <- c(NROW(x), NROW(y))
    
    boot.obj <- boot(data = z, statistic = Tn3,
                     sim = "permutation", R = 999, sizes = N)
    boot.obj
    
    tb <- c(boot.obj$t, boot.obj$t0)
    mean(tb >= boot.obj$t0)
    
    hist(tb, freq=FALSE, main="",
         xlab="replicates of T(n,3) statistic")
    points(boot.obj$t0, 0, cex=1, pch=16)
    
    


### Example 10.7 (Two-sample energy statistic)

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


### Example 10.8 (Two-sample energy test)

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

### Example 10.9 (k-sample energy distances)

    z <- iris[ , 1:4]
    dst <- dist(z)
    energy::edist(dst, sizes = c(50, 50, 50), distance = TRUE)

### Example 10.10 (Distance Components (disco))

    set.seed(413)
    energy::disco(iris[ , 1:4], factors = iris$Species, R = 999)
    
### Example 10.11 (Power Comparison)
    
    # results of several simulations summarized in Table
    
### Example 10.12 (Distance covariance statistic)

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

### Example 10.13 (Distance correlation statistic)

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
    

### Example 10.14 (Distance covariance test)

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
