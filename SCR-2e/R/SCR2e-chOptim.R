#######################################################
###       Statistical Computing with R  2e          ###
###       Maria L. Rizzo                            ###
###       Chapman & Hall/CRC The R Series           ###
###       ISBN 9781466553323 - CAT# K15269          ###
###       January 2019                              ###
###                                                 ###
###       R code for Chapter 14                     ###
###       Optimization                              ###
#######################################################


### Example 14.1 (One-dimensional optimization with optimize)

    f <- function(x)
        log(x + log(x))/log(1+x)

    curve(f(x), from = 2, to = 15, ylab = "f(x)")
    
    optimize(f, lower = 4, upper = 8, maximum = TRUE)


### Example 14.2 (MLE: Gamma distribution)

    m <- 20000
    est <- matrix(0, m, 2)
    n <- 200
    r <- 5
    lambda <- 2

    obj <- function(lambda, xbar, logx.bar) {
        digamma(lambda * xbar) - logx.bar - log(lambda)
        }

    for (i in 1:m) {
        x <- rgamma(n, shape=r, rate=lambda)
        xbar <- mean(x)
        u <- uniroot(obj, lower = .001, upper = 10e5,
                xbar = xbar, logx.bar = mean(log(x)))
        lambda.hat <- u$root
        r.hat <- xbar * lambda.hat
        est[i, ] <- c(r.hat, lambda.hat)
    }

    ML <- colMeans(est)

    hist(est[, 1], breaks="scott", freq=FALSE,
        xlab="r", main="")
    points(ML[1], 0, cex=1.5, pch=20)
    hist(est[, 2], breaks="scott", freq=FALSE,
         xlab=bquote(lambda), main="")
    points(ML[2], 0, cex=1.5, pch=20)


### Example 14.3 (Two-dimensional optimization with optim)

    LL <- function(theta, sx, slogx, n) {
        r <- theta[1]
        lambda <- theta[2]
        loglik <- n * r * log(lambda) + (r - 1) * slogx -
            lambda * sx - n * log(gamma(r))
        - loglik
        }

    n <- 200
    r <- 5;    lambda <- 2
    x <- rgamma(n, shape=r, rate=lambda)

    optim(c(1,1), LL, sx=sum(x), slogx=sum(log(x)), n=n)

    mlests <- replicate(20000, expr = {
      x <- rgamma(200, shape = 5, rate = 2)
      optim(c(1,1), LL, sx=sum(x), slogx=sum(log(x)), n=n)$par
      })
    colMeans(t(mlests))


### Example 14.4 (MLE for a quadratic form)

    LL <- function(lambda, y) {
        lambda3 <- 1 - sum(lambda)
        f1 <- dgamma(y, shape=1/2, rate=1/(2*lambda[1]))
        f2 <- dgamma(y, shape=1/2, rate=1/(2*lambda[2]))
        f3 <- dgamma(y, shape=1/2, rate=1/(2*lambda3))
        f <- f1/3 + f2/3 + f3/3   #density of mixture
        #returning -loglikelihood
        return( -sum(log(f)))
        }

    set.seed(543)
    m <- 2000
    lambda <- c(.6, .25, .15)  #rate is 1/(2lambda)
    lam <- sample(lambda, size = 2000, replace = TRUE)
    y <- rgamma(m, shape = .5, rate = 1/(2*lam))

    opt <- optim(c(.5,.3), LL, y=y)
    theta <- c(opt$par, 1 - sum(opt$par))

    as.data.frame(unlist(opt))

    theta


### Example 14.5 (EM algorithm for a mixture model)

    set.seed(543)
    lambda <- c(.6, .25, .15)  #rate is 1/(2lambda)
    lam <- sample(lambda, size = 2000, replace = TRUE)
    y <- rgamma(m, shape = .5, rate = 1/(2*lam))

    N <- 10000                #max. number of iterations
    L <- c(.5, .4, .1)        #initial est. for lambdas
    tol <- .Machine$double.eps^0.5
    L.old <- L + 1

    for (j in 1:N) {
        f1 <- dgamma(y, shape=1/2, rate=1/(2*L[1]))
        f2 <- dgamma(y, shape=1/2, rate=1/(2*L[2]))
        f3 <- dgamma(y, shape=1/2, rate=1/(2*L[3]))
        py <- f1 / (f1 + f2 + f3) #posterior prob y from 1
        qy <- f2 / (f1 + f2 + f3) #posterior prob y from 2
        ry <- f3 / (f1 + f2 + f3) #posterior prob y from 3

        mu1 <- sum(y * py) / sum(py) #update means
        mu2 <- sum(y * qy) / sum(qy)
        mu3 <- sum(y * ry) / sum(ry)
        L <- c(mu1, mu2, mu3)  #update lambdas
        L <- L / sum(L)

        if (sum(abs(L - L.old)/L.old) < tol) break
        L.old <- L
    }

    print(list(lambda = L/sum(L), iter = j, tol = tol))


### Example 14.6 (Simplex algorithm)

    library(boot)   #for simplex function
    A1 <- rbind(c(-2, 1, 1), c(4, -1, 3))
    b1 <- c(1, 3)
    a <- c(2, 2, 3)
    simplex(a = a, A1 = A1, b1 = b1, maxi = TRUE)
    detach(package:boot)


### Example 14.7 (Solving the Morra game)

    solve.game <- function(A) {
        #solve the two player zero-sum game by simplex method
        #optimize for player 1, then player 2
        #maximize v subject to ...
        #let x strategies 1:m, and put v as extra variable
        #A1, the <= constraints
        #
        min.A <- min(A)
        A <- A - min.A   #so that v >= 0
        max.A <- max(A)
        A <- A / max(A)
        m <- nrow(A)
        n <- ncol(A)
        it <- n^3
        a <- c(rep(0, m), 1) #objective function
        A1 <- -cbind(t(A), rep(-1, n)) #constraints <=
        b1 <- rep(0, n)
        A3 <- t(as.matrix(c(rep(1, m), 0))) #constraints sum(x)=1
        b3 <- 1
        sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
                    maxi=TRUE, n.iter=it)
        #the 'solution' is [x1,x2,...,xm | value of game]
        #
        #minimize v subject to ...
        #let y strategies 1:n, with v as extra variable
        a <- c(rep(0, n), 1) #objective function
        A1 <- cbind(A, rep(-1, m)) #constraints <=
        b1 <- rep(0, m)
        A3 <- t(as.matrix(c(rep(1, n), 0))) #constraints sum(y)=1
        b3 <- 1
        sy <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
                    maxi=FALSE, n.iter=it)

        soln <- list("A" = A * max.A + min.A,
                    "x" = sx$soln[1:m],
                    "y" = sy$soln[1:n],
                    "v" = sx$soln[m+1] * max.A + min.A)
        soln
        }


    #enter the payoff matrix
    A <- matrix(c(  0,-2,-2,3,0,0,4,0,0,
                    2,0,0,0,-3,-3,4,0,0,
                    2,0,0,3,0,0,0,-4,-4,
                    -3,0,-3,0,4,0,0,5,0,
                    0,3,0,-4,0,-4,0,5,0,
                    0,3,0,0,4,0,-5,0,-5,
                    -4,-4,0,0,0,5,0,0,6,
                    0,0,4,-5,-5,0,0,0,6,
                    0,0,4,0,0,5,-6,-6,0), 9, 9)

    library(boot)  #needed for simplex function

    s <- solve.game(A)
    round(cbind(s$x, s$y), 7)
