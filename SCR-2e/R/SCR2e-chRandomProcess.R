#######################################################
###       Statistical Computing with R  2e          ###
###       Maria L. Rizzo                            ###
###       Chapman & Hall/CRC The R Series           ###
###       ISBN 9781466553323 - CAT# K15269          ###
###       January 2019                              ###
###                                                 ###
###       R code for Chapter 4                      ###
###       Generating Random Processes               ###
#######################################################


### Example 4.1 (Poisson process)

    lambda <- 2
    t0 <- 3
    Tn <- rexp(100, lambda)       #interarrival times
    Sn <- cumsum(Tn)              #arrival times
    n <- min(which(Sn > t0))      #arrivals+1 in [0, t0]


### Example 4.2 (Poisson process, cont.)

    lambda <- 2
    t0 <- 3
    upper <- 100
    pp <- numeric(10000)
    for (i in 1:10000) {
        N <- rpois(1, lambda * upper)
        Un <- runif(N, 0, upper)      #unordered arrival times
        Sn <- sort(Un)                #arrival times
        n <- min(which(Sn > t0))      #arrivals+1 in [0, t0]
        pp[i] <- n - 1                #arrivals in [0, t0]
        }

    #alternately, the loop can be replaced by replicate function
    pp <- replicate(10000, expr = {
        N <- rpois(1, lambda * upper)
        Un <- runif(N, 0, upper)      #unordered arrival times
        Sn <- sort(Un)                #arrival times
        n <- min(which(Sn > t0))      #arrivals+1 in [0, t0]
        n - 1  })                     #arrivals in [0, t0]

    c(mean(pp), var(pp))


### Example 4.3 (Nonhomogeneous Poisson process)

    lambda <- 3
    upper <- 100
    N <- rpois(1, lambda * upper)
    Tn <- rexp(N, lambda)
    Sn <- cumsum(Tn)
    Un <- runif(N)
    keep <- (Un <= cos(Sn)^2)    #indicator, as logical vector
    Sn[keep]

    round(Sn[keep], 4)


### Example 4.4 (Renewal process)

    t0 <- 5
    Tn <- rgeom(100, prob = .2)   #interarrival times
    Sn <- cumsum(Tn)              #arrival times
    n <- min(which(Sn > t0))      #arrivals+1 in [0, t0]

    Nt0 <- replicate(1000, expr = {
        Sn <- cumsum(rgeom(100, prob = .2))
        min(which(Sn > t0)) - 1
        })
    table(Nt0)/1000
    Nt0

    t0 <- seq(0.1, 30, .1)
    mt <- numeric(length(t0))

    for (i in 1:length(t0)) {
        mt[i] <- mean(replicate(1000,
        {
        Sn <- cumsum(rgeom(100, prob = .2))
        min(which(Sn > t0[i])) - 1
        }))
    }
    plot(t0, mt, type = "l", xlab = "t", ylab = "mean")
    abline(0, .25)
    

### Example 4.5 (Symmetric random walk)

    n <- 400
    incr <- sample(c(-1, 1), size = n, replace = TRUE)
    S <- as.integer(c(0, cumsum(incr)))
    plot(0:n, S, type = "l", main = "", xlab = "i")

### Example 4.6 (Generator for the time until return to origin)
    
    set.seed(12345)
    
    #compute the probabilities directly
    n <- 1:10000
    p2n <- exp(lgamma(2*n-1)
              - log(n) - (2*n-1)*log(2) - 2*lgamma(n))
    #or compute using dbinom
    P2n <- (.5/n) * dbinom(n-1, size = 2*n-2, prob = 0.5)
    pP2n <- cumsum(P2n)

    #given n compute the time of the last return to 0 in (0,n]
    n <- 200
    sumT <- 0
    while (sumT <= n) {
        u <- runif(1)
        s <- sum(u > pP2n)
        if (s == length(pP2n))
            warning("T is truncated")
        Tj <- 2 * (1 + s)
        #print(c(Tj, sumT))
        sumT <- sumT + Tj
        }
    sumT - Tj

### Example 4.7 (Brownian motion)
    
    
    simBM <- function(n, T) {
      times <- seq(0, T, length = n+1)
      z <- rnorm(n)
      w <- rep(0, n)
      s <- sqrt(diff(times))
      for (k in 2:n) {
        w[k] <- w[k-1] + s[k] * z[k]
      }
      return (list(w=w, t=times))
    }
    
    set.seed(1)
    n <- 200
    x1 <- simBM(n, 1)
    x2 <- simBM(n, 1)
    x3 <- simBM(n, 1)
    r <- range(c(x1$w, x2$w, x3$w))
    plot(x1$w, type="l", main="", xlab="t", ylab="W", ylim=r)
    lines(x2$w, lty=2)
    lines(x3$w, lty=3)
    
    interpBM <- function(w, t0, times) {
      k1 <- sum(times < t0)
      k <- k1 + 1
      b <- (t0 - times[k1]) / (times[k] - times[k1])
      return (w[k1] + b * (w[k] - w[k1]))
    }
    
    plot(x1$t[1:10], x1$w[1:10], type="b", main="", xlab="t", ylab="W")
    tmids <- x1$t + 0.0025
    for (i in 1:10) {
      w <- interpBM(x1$w, tmids[i], x1$t)
      points(tmids[i], w, pch=2)
    }
    
    legend("topleft", c("Generated W", "Interpolated W"), pch=c(1,2), lty=1, bty="n")
    