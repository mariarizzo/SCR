#######################################################
###       Statistical Computing with R  2e          ###
###       Maria L. Rizzo                            ###
###       Chapman & Hall/CRC The R Series           ###
###       ISBN 9781466553323 - CAT# K15269          ###
###       January 2019                              ###
###                                                 ###
###       R code for Chapter 7                      ###
###       Monte Carlo Methods in Inference          ###
#######################################################

# packages to install: energy, ggplot2

### Example 7.1 (Basic Monte Carlo estimation)

    m <- 1000
    g <- numeric(m)
    for (i in 1:m) {
        x <- rnorm(2)
        g[i] <- abs(x[1] - x[2])
    }
    est <- mean(g)
    est


### Example 7.2 (Estimating the MSE of a trimmed mean)

    n <- 20
    m <- 1000
    tmean <- numeric(m)
    for (i in 1:m) {
        x <- sort(rnorm(n))
        tmean[i] <- sum(x[2:(n-1)]) / (n-2)
        }
    mse <- mean(tmean^2)
    mse
    sqrt(sum((tmean - mean(tmean))^2)) / m    #se

    n <- 20
    m <- 1000
    tmean <- numeric(m)
    for (i in 1:m) {
        x <- sort(rnorm(n))
        tmean[i] <- median(x)
        }
    mse <- mean(tmean^2)
    mse
    sqrt(sum((tmean - mean(tmean))^2)) / m    #se


### Example 7.3 (MSE of a trimmed mean, cont.)

     set.seed(522)
     n <- 20
     K <- n/2 - 1
     m <- 1000
     mse <- matrix(0, n/2, 6)

     trimmed.mse <- function(n, m, k, p) {
         #MC est of mse for k-level trimmed mean of
         #contaminated normal pN(0,1) + (1-p)N(0,100)
         tmean <- numeric(m)
         for (i in 1:m) {
             sigma <- sample(c(1, 10), size = n,
                 replace = TRUE, prob = c(p, 1-p))
             x <- sort(rnorm(n, 0, sigma))
             tmean[i] <- sum(x[(k+1):(n-k)]) / (n-2*k)
             }
         mse.est <- mean(tmean^2)
         se.mse <- sqrt(mean((tmean-mean(tmean))^2)) / sqrt(m)
         return(c(mse.est, se.mse))
     }

    for (k in 0:K) {
        mse[k+1, 1:2] <- trimmed.mse(n=n, m=m, k=k, p=1.0)
        mse[k+1, 3:4] <- trimmed.mse(n=n, m=m, k=k, p=.95)
        mse[k+1, 5:6] <- trimmed.mse(n=n, m=m, k=k, p=.9)
    }

### Example 7.4 (Confidence interval for variance)

    n <- 20
    alpha <- .05
    x <- rnorm(n, mean=0, sd=2)
    UCL <- (n-1) * var(x) / qchisq(alpha, df=n-1)


### Example 7.5 (MC estimate of confidence level)

    n <- 20
    alpha <- .05
    UCL <- replicate(1000, expr = {
        x <- rnorm(n, mean = 0, sd = 2)
        (n-1) * var(x) / qchisq(alpha, df = n-1)
        } )
    #count the number of intervals that contain sigma^2=4
    sum(UCL > 4)
    #or compute the mean to get the confidence level
    mean(UCL > 4)


### Example 7.6 (Empirical confidence level)

    n <- 20
    alpha <- .05
    UCL <- replicate(1000, expr = {
        x <- rchisq(n, df = 2)
        (n-1) * var(x) / qchisq(alpha, df = n-1)
        } )
    sum(UCL > 4)
    mean(UCL > 4)


### Example 7.7 (Empirical Type I error rate)

    n <- 20
    alpha <- .05
    mu0 <- 500
    sigma <- 100

    m <- 10000          #number of replicates
    p <- numeric(m)     #storage for p-values
    for (j in 1:m) {
        x <- rnorm(n, mu0, sigma)
        ttest <- t.test(x, alternative = "greater", mu = mu0)
        p[j] <- ttest$p.value
        }

    p.hat <- mean(p < alpha)
    se.hat <- sqrt(p.hat * (1 - p.hat) / m)
    print(c(p.hat, se.hat))


### Example 7.8 (Skewness test of normality)

    n <- c(10, 20, 30, 50, 100, 500) #sample sizes
    cv <- qnorm(.975, 0, sqrt(6/n))  #crit. values for each n

    sk <- function(x) {
        #computes the sample skewness coeff.
        xbar <- mean(x)
        m3 <- mean((x - xbar)^3)
        m2 <- mean((x - xbar)^2)
        return( m3 / m2^1.5 )
    }

    #n is a vector of sample sizes
    #we are doing length(n) different simulations

    p.reject <- numeric(length(n)) #to store sim. results
    m <- 10000                     #num. repl. each sim.

    for (i in 1:length(n)) {
        sktests <- numeric(m)       #test decisions
        for (j in 1:m) {
            x <- rnorm(n[i])
            #test decision is 1 (reject) or 0
            sktests[j] <- as.integer(abs(sk(x)) >= cv[i] )
            }
        p.reject[i] <- mean(sktests) #proportion rejected
    }

    p.reject


### Example 7.9 (Empirical power)

    #set.seed(521)
    n <- 20
    m <- 1000
    mu0 <- 500
    sigma <- 100
    mu <- c(seq(450, 650, 10))  #alternatives
    M <- length(mu)
    power <- numeric(M)
    for (i in 1:M) {
        mu1 <- mu[i]
        pvalues <- replicate(m, expr = {
            #simulate under alternative mu1
            x <- rnorm(n, mean = mu1, sd = sigma)
            ttest <- t.test(x,
                     alternative = "greater", mu = mu0)
            ttest$p.value  } )
        power[i] <- mean(pvalues <= .05)
    }
    
    se <- sqrt(power * (1-power) / m)
    df <- data.frame(mean=mu, power=power, upper=power+2*se, lower=power-2*se)
    
    library(ggplot2)
    ggplot(df, aes(x=mean, y=power)) +
      geom_line() +
      geom_vline(xintercept=500, lty=2) +
      geom_hline(yintercept=c(0,.05), lty=1:2) +
      geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.2, lwd=1.5)
    

### Example 7.10 (Power of the skewness test of normality)

    #set.seed(111)
    alpha <- .1
    n <- 30
    m <- 2500
    epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
    N <- length(epsilon)
    pwr <- numeric(N)
    #critical value for the skewness test
    cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))

    for (j in 1:N) {           #for each epsilon
        e <- epsilon[j]
        sktests <- numeric(m)
        for (i in 1:m) {       #for each replicate
            sigma <- sample(c(1, 10), replace = TRUE,
                size = n, prob = c(1-e, e))
            x <- rnorm(n, 0, sigma)
            sktests[i] <- as.integer(abs(sk(x)) >= cv)
            }
        pwr[j] <- mean(sktests)
        }

    se <- sqrt(pwr * (1-pwr) / m)
    df <- data.frame(epsilon=epsilon, power=pwr, upper=pwr+2*se, lower=pwr-2*se)
    
    #plot power vs epsilon
    library(ggplot2)
    ggplot(df, aes(x=epsilon, y=power)) +
      geom_line() + labs(x=bquote(epsilon)) +
      geom_hline(yintercept=.1, lty=2) +
      geom_pointrange(aes(ymin=lower, ymax=upper))
    

### Example 7.11 (Power comparison of tests of normality)

    #only one loop, for epsilon=0.1, was shown in the text
    #the simulation below takes several minutes to run

    # initialize input and output
    library(energy)
    alpha <- .1
    n <- 30
    m <- 500        #try small m for a trial run
    test1 <- test2 <- test3 <- numeric(m)

    #critical value for the skewness test
    cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
    sim <- matrix(0, 11, 4)

    # estimate power
    for (i in 0:10) {
        epsilon <- i * .1
        for (j in 1:m) {
            e <- epsilon
            sigma <- sample(c(1, 10), replace = TRUE,
                size = n, prob = c(1-e, e))
            x <- rnorm(n, 0, sigma)
            test1[j] <- as.integer(abs(sk(x)) >= cv)
            test2[j] <- as.integer(
                        shapiro.test(x)$p.value <= alpha)
            test3[j] <- as.integer(
                        mvnorm.etest(x, R=200)$p.value <= alpha)
            }
        print(c(epsilon, mean(test1), mean(test2), mean(test3)))
        sim[i+1, ] <- c(epsilon, mean(test1), mean(test2), mean(test3))
    }
    detach(package:energy)

    # plot the empirical estimates of power
    plot(sim[,1], sim[,2], ylim = c(0, 1), type = "l",
        xlab = bquote(epsilon), ylab = "power")
    lines(sim[,1], sim[,3], lty = 2)
    lines(sim[,1], sim[,4], lty = 4)
    abline(h = alpha, lty = 3)
    legend("topright", 1, c("skewness", "S-W", "energy"),
        lty = c(1,2,4), inset = .02)


### Example 7.12 (Count Five test statistic)

    x1 <- rnorm(20, 0, sd = 1)
    x2 <- rnorm(20, 0, sd = 1.5)
    y <- c(x1, x2)

    group <- rep(1:2, each = length(x1))
    boxplot(y ~ group, boxwex = .3, xlim = c(.5, 2.5), main = "")
    points(group, y)

    # now identify the extreme points
    range(x1)
    range(x2)

    i <- which(x1 < min(x2))
    j <- which(x2 > max(x1))

    x1[i]
    x2[j]

    out1 <- sum(x1 > max(x2)) + sum(x1 < min(x2))
    out2 <- sum(x2 > max(x1)) + sum(x2 < min(x1))
    max(c(out1, out2))


### Example 7.13 (Count Five test statistic, cont.)

    maxout <- function(x, y) {
        X <- x - mean(x)
        Y <- y - mean(y)
        outx <- sum(X > max(Y)) + sum(X < min(Y))
        outy <- sum(Y > max(X)) + sum(Y < min(X))
        return(max(c(outx, outy)))
    }

    n1 <- n2 <- 20
    mu1 <- mu2 <- 0
    sigma1 <- sigma2 <- 1
    m <- 1000

    # generate samples under H0
    stat <- replicate(m, expr={
        x <- rnorm(n1, mu1, sigma1)
        y <- rnorm(n2, mu2, sigma2)
        maxout(x, y)
        })
    print(cumsum(table(stat)) / m)
    print(quantile(stat, c(.8, .9, .95)))


### Example 7.14 (Count Five test)

    count5test <- function(x, y) {
        X <- x - mean(x)
        Y <- y - mean(y)
        outx <- sum(X > max(Y)) + sum(X < min(Y))
        outy <- sum(Y > max(X)) + sum(Y < min(X))
        # return 1 (reject) or 0 (do not reject H0)
        return(as.integer(max(c(outx, outy)) > 5))
    }

    n1 <- n2 <- 20
    mu1 <- mu2 <- 0
    sigma1 <-  sigma2 <- 1
    m <- 10000
    tests <- replicate(m, expr = {
        x <- rnorm(n1, mu1, sigma1)
        y <- rnorm(n2, mu2, sigma2)
        x <- x - mean(x)  #centered by sample mean
        y <- y - mean(y)
        count5test(x, y)
        } )

    alphahat <- mean(tests)
    print(alphahat)


### Example 7.15 (Count Five test, cont.)

    n1 <- 20
    n2 <- 30
    mu1 <- mu2 <- 0
    sigma1 <- sigma2 <- 1
    m <- 10000

    alphahat <- mean(replicate(m, expr={
        x <- rnorm(n1, mu1, sigma1)
        y <- rnorm(n2, mu2, sigma2)
        x <- x - mean(x)  #centered by sample mean
        y <- y - mean(y)
        count5test(x, y)
        }))

    print(alphahat)


### Example 7.16 (Count Five, cont.)

    # generate samples under H1 to estimate power
    sigma1 <- 1
    sigma2 <- 1.5

    power <- mean(replicate(m, expr={
        x <- rnorm(20, 0, sigma1)
        y <- rnorm(20, 0, sigma2)
        count5test(x, y)
        }))

    print(power)
