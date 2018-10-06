### Statistical Computing with R/CRC Press
###
### Revised Examples 8.4, 8.5, 8.6 
### package knnFinder was withdrawn from CRAN so these examples
### are revised to replace the nearest neighbor function nn (knnFinder)
### with ann() in package yaImpute
#
# Notes regarding ann() function:
# 1. Squared Euclidean distances are returned (docs say Euclidean distances)
# 2. Rows of $knnIndexDist correspond to ref
# 3. Columns 1:k of $knnIndexDist give the NN indices
# 4. For the first through r nearest neighbors, set k=r+1 in ann()
#    and ignore col. 1 in $knnIndexDist
# 5. exact vs approximate kd-tree search
#    exact kd-tree search in ann with tree.type="kd" and eps=0;
#    Set eps>0 for approximate search (nn was exact)
###########################################################################

### Example 8.4R (Finding nearest neighbors)
### Revised version using ann (yaImpute)

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


### Example 8.5R (Nearest neighbor statistic)
### Revised version using ann (yaImpute)

    library(yaImpute)
    attach(chickwts)
    x <- as.vector(weight[feed == "sunflower"])
    y <- as.vector(weight[feed == "linseed"])
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



### Example 8.6R (Nearest neighbor test)
### Revised version using ann (yaImpute)

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
    x <- as.vector(weight[feed == "sunflower"])
    y <- as.vector(weight[feed == "linseed"])
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


