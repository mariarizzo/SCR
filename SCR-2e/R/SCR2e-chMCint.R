#######################################################
###    Statistical Computing with R  2e             ###
###    Maria L. Rizzo                               ###
###    Chapman & Hall/CRC The R Series              ###
###    ISBN 9781466553323 - CAT# K15269             ###
###       March 6, 2019                             ###
###                                                 ###
###  R code for Chapter 6                           ###
###  Monte Carlo Integration and Variance Reduction ###
#######################################################

### Example 6.1 (Simple Monte Carlo integration)

    m <- 10000
    x <- runif(m)
    theta.hat <- mean(exp(-x))
    print(theta.hat)
    print(1 - exp(-1))


### Example 6.2 (Simple Monte Carlo integration, cont.)

    m <- 10000
    x <- runif(m, min=2, max=4)
    theta.hat <- mean(exp(-x)) * 2
    print(theta.hat)
    print(exp(-2) - exp(-4))


### Example 6.3 (Monte Carlo integration, unbounded interval)

    x <- seq(.1, 2.5, length = 10)
    m <- 10000
    u <- runif(m)
    cdf <- numeric(length(x))
    for (i in 1:length(x)) {
        g <- x[i] * exp(-(u * x[i])^2 / 2)
        cdf[i] <- mean(g) / sqrt(2 * pi) + 0.5
    }

    Phi <- pnorm(x)
    print(round(rbind(x, cdf, Phi), 3))


### Example 6.4 (Example 6.3, cont.)

    x <- seq(.1, 2.5, length = 10)
    m <- 10000
    z <- rnorm(m)
    dim(x) <- length(x)
    p <- apply(x, MARGIN = 1,
             FUN = function(x, z) {mean(z < x)}, z = z)

    Phi <- pnorm(x)
    print(round(rbind(x, p, Phi), 3))


### Example 6.5 (Error bounds for MC integration)

    x <- 2
    m <- 10000
    z <- rnorm(m)
    g <- (z < x)  #the indicator function
    v <- mean((g - mean(g))^2) / m
    cdf <- mean(g)
    c(cdf, v)
    c(cdf - 1.96 * sqrt(v), cdf + 1.96 * sqrt(v))


### Example 6.6 (Antithetic variables)

    MC.Phi <- function(x, R = 10000, antithetic = TRUE) {
        u <- runif(R/2)
        if (!antithetic) v <- runif(R/2) else
            v <- 1 - u
        u <- c(u, v)
        cdf <- numeric(length(x))
        for (i in 1:length(x)) {
            g <- x[i] * exp(-(u * x[i])^2 / 2)
            cdf[i] <- mean(g) / sqrt(2 * pi) + 0.5
        }
        cdf
    }


    x <- seq(.1, 2.5, length=5)
    Phi <- pnorm(x)
    set.seed(123)
    MC1 <- MC.Phi(x, anti = FALSE)
    set.seed(123)
    MC2 <- MC.Phi(x)
    print(round(rbind(x, MC1, MC2, Phi), 5))


    m <- 1000
    MC1 <- MC2 <- numeric(m)
    x <- 1.95
    for (i in 1:m) {
        MC1[i] <- MC.Phi(x, R = 1000, anti = FALSE)
        MC2[i] <- MC.Phi(x, R = 1000)
    }

    print(sd(MC1))
    print(sd(MC2))
    print((var(MC1) - var(MC2))/var(MC1))


### Example 6.7 (Control variate)

    m <- 10000
    a <- - 12 + 6 * (exp(1) - 1)
    U <- runif(m)
    T1 <- exp(U)                  #simple MC
    T2 <- exp(U) + a * (U - 1/2)  #controlled

    mean(T1)
    mean(T2)
    (var(T1) - var(T2)) / var(T1)


### Example 6.8 (MC integration using control variates)

    f <- function(u)
        exp(-.5)/(1+u^2)

    g <- function(u)
        exp(-u)/(1+u^2)

    set.seed(510) #needed later
    u <- runif(10000)
    B <- f(u)
    A <- g(u)

    cor(A, B)
    a <- -cov(A,B) / var(B)    #est of c*
    a

    m <- 100000
    u <- runif(m)
    T1 <- g(u)
    T2 <- T1 + a * (f(u) - exp(-.5)*pi/4)

    c(mean(T1), mean(T2))
    c(var(T1), var(T2))
    (var(T1) - var(T2)) / var(T1)


### Example 6.9 (Control variate and regression)

    
    set.seed(510)
    mu <- exp(-.5)*pi/4
    u <- runif(10000)
    f <- exp(-.5)/(1+u^2)
    g <- exp(-u)/(1+u^2)
    L <- lm(g ~ f)
    L
    c.star <- - L$coeff[2]
    c.star 
    
    theta.hat <- sum(L$coeff * c(1, mu))  #pred. value at mu
    theta.hat
    summary(L)$sigma^2
    summary(L)$r.squared

### Example 6.10 (Control variates and multiple regression)
    
    # Example 6.9 continued with a second control variate
    # and multiple regression to estimate vector c*
    
    u <- runif(10000)
    f1 <- exp(-.5) / (1+u^2)
    f2 <- exp(-u) / (1-exp(-1))
    g <- exp(-u) / (1+u^2)
    
    L2 <- lm(g ~ f1 + f2)
    
    L2$coeff
    c.star <-  - L2$coeff[2:3]
    c.star
    mu1 <- exp(-.5)*pi/4
    mu2 <- 1
    mu <- c(mu1, mu2)
    
    # theta.hat is the predicted response at mu
    # alternately can use predict.lm method
    
    theta.hat <- sum(L2$coeff * c(1, mu))  #pred. value at mu
    theta.hat
    
    ## alternately
    df <- data.frame(f1=mu1, f2=mu2)
    theta.hat <-predict(L2, df)
    
    # MSE / n is the est. variance of the control estimator
    MSE <- summary(L2)$sigma^2 
    MSE
    sqrt(MSE / 10000)
    
    
    # compare with the previous estimates using
    # naive MC and control variate f1(u)
    # var1=.060231423 var2=.003124814
    var0 <- 0.060231423  #naive MC
    var1 <- 0.003117644   #controlled estimator with f1
    var2 <- MSE  #new estimator
    
    # percent reduction in variance
    # it is a weighted average of R^2 values
    # so easier to compute directly
    
    100 * (var0 - var1) / var0
    100 * (var1 - var2) / var1
    100 * (var0 - var2) / var0
    
### Example 6.11 (Choice of the importance function)
    #code for plot is at the end of the file

    m <- 10000
    theta.hat <- se <- numeric(5)
    g <- function(x) {
        exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
        }

    x <- runif(m)     #using f0
    fg <- g(x)
    theta.hat[1] <- mean(fg)
    se[1] <- sd(fg)

    x <- rexp(m, 1)   #using f1
    fg <- g(x) / exp(-x)
    theta.hat[2] <- mean(fg)
    se[2] <- sd(fg)

    x <- rcauchy(m)   #using f2
    i <- c(which(x > 1), which(x < 0))
    x[i] <- 2  #to catch overflow errors in g(x)
    fg <- g(x) / dcauchy(x)
    theta.hat[3] <- mean(fg)
    se[3] <- sd(fg)

    u <- runif(m)     #f3, inverse transform method
    x <- - log(1 - u * (1 - exp(-1)))
    fg <- g(x) / (exp(-x) / (1 - exp(-1)))
    theta.hat[4] <- mean(fg)
    se[4] <- sd(fg)

    u <- runif(m)    #f4, inverse transform method
    x <- tan(pi * u / 4)
    fg <- g(x) / (4 / ((1 + x^2) * pi))
    theta.hat[5] <- mean(fg)
    se[5] <- sd(fg)

    rbind(theta.hat, se / sqrt(m))


### Example 6.12 (Example 6.11, cont.)

    M <- 20   #number of replicates
    T2 <- numeric(4)
    estimates <- matrix(0, 10, 2)

    g <- function(x) {
        exp(-x - log(1+x^2)) * (x > 0) * (x < 1) }

    for (i in 1:10) {
        estimates[i, 1] <- mean(g(runif(M)))
        T2[1] <- mean(g(runif(M/4, 0, .25)))
        T2[2] <- mean(g(runif(M/4, .25, .5)))
        T2[3] <- mean(g(runif(M/4, .5, .75)))
        T2[4] <- mean(g(runif(M/4, .75, 1)))
        estimates[i, 2] <- mean(T2)
    }

    estimates
    apply(estimates, 2, mean)
    apply(estimates, 2, var)


### Example 6.13 (Examples 6.11-6.12, cont.)

    M <- 10000  #number of replicates
    k <- 10     #number of strata
    r <- M / k  #replicates per stratum
    N <- 50     #number of times to repeat the estimation
    T2 <- numeric(k)
    estimates <- matrix(0, N, 2)

    g <- function(x) {
        exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
        }

    for (i in 1:N) {
        estimates[i, 1] <- mean(g(runif(M)))
        for (j in 1:k)
            T2[j] <- mean(g(runif(M/k, (j-1)/k, j/k)))
        estimates[i, 2] <- mean(T2)
    }

    apply(estimates, 2, mean)
    apply(estimates, 2, var)



### Plot importance functions in Figures 6.1(a) and 6.1.(b)

    #par(ask = TRUE) #uncomment to pause between graphs

    x <- seq(0, 1, .01)
    w <- 2
    f1 <- exp(-x)
    f2 <- (1 / pi) / (1 + x^2)
    f3 <- exp(-x) / (1 - exp(-1))
    f4 <- 4 / ((1 + x^2) * pi)
    g <- exp(-x) / (1 + x^2)

    #for color change lty to col

    #figure (a)
    plot(x, g, type = "l", main = "", ylab = "",
         ylim = c(0,2), lwd = w)
    lines(x, g/g, lty = 2, lwd = w)
    lines(x, f1, lty = 3, lwd = w)
    lines(x, f2, lty = 4, lwd = w)
    lines(x, f3, lty = 5, lwd = w)
    lines(x, f4, lty = 6, lwd = w)
    legend("topright", legend = c("g", 0:4),
           lty = 1:6, lwd = w, inset = 0.02)

    #figure (b)
    plot(x, g, type = "l", main = "", ylab = "",
        ylim = c(0,3.2), lwd = w, lty = 2)
    lines(x, g/f1, lty = 3, lwd = w)
    lines(x, g/f2, lty = 4, lwd = w)
    lines(x, g/f3, lty = 5, lwd = w)
    lines(x, g/f4, lty = 6, lwd = w)
    legend("topright", legend = c(0:4),
           lty = 2:6, lwd = w, inset = 0.02)

    