#######################################################
###       Statistical Computing with R  2e          ###
###       Maria L. Rizzo                            ###
###       Chapman & Hall/CRC The R Series           ###
###       ISBN 9781466553323 - CAT# K15269          ###
###       March 6, 2019                             ###
###                                                 ###
###       R code for Chapter 11                     ###
###       Markov Chain Monte Carlo                  ###
#######################################################


    
### Example 11.1 (Metropolis-Hastings sampler)

    f <- function(x, sigma) {
        if (any(x < 0)) return (0)
        stopifnot(sigma > 0)
        return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
    }

    m <- 10000
    sigma <- 4
    x <- numeric(m)
    x[1] <- rchisq(1, df=1)
    k <- 0
    u <- runif(m)

    for (i in 2:m) {
        xt <- x[i-1]
        y <- rchisq(1, df = xt)
        num <- f(y, sigma) * dchisq(xt, df = y)
        den <- f(xt, sigma) * dchisq(y, df = xt)
        if (u[i] <= num/den) x[i] <- y else {
             x[i] <- xt
             k <- k+1     #y is rejected
             }
        }

    print(k)

    index <- 5000:5500
    y1 <- x[index]
    plot(index, y1, type="l", main="", ylab="x")


### Example 11.2 (Example 11.1, cont.)

    b <- 2001      #discard the burnin sample
    y <- x[b:m]
    a <- ppoints(100)
    QR <- sigma * sqrt(-2 * log(1 - a))  #quantiles of Rayleigh
    Q <- quantile(y, a)
    
    qqplot(QR, Q, main="", cex=.5,
           xlab="Rayleigh Quantiles", ylab="Sample Quantiles")
    abline(0, 1)
    
    hist(y, breaks="scott", main="", xlab="", freq=FALSE)
    lines(QR, f(QR, 4))
    
    
### Example 11.3 (Expected lifetime beta-binomial model)
    
    f.mu <- function(x) {
      exp(- lbeta(432, 5) + 431 * log(x) -  
            431 * log(1 + x) - 6 * log(1 + x))
    }
    curve(f.mu(x), from=0, to=400, xlab="hours", ylab="")

    fr <- function(x, y) {
      a <- 431 * (log(y) - log(x))
      b <- 437 * (log(1+x) - log(1+y))
      return(exp(a + b))
    }

    set.seed(2016)
    m <- 10000
    x <- numeric(m)
    x[1] <- rchisq(1, df=1)  #initialize chain
    k <- 0
    u <- runif(m)
    for (i in 2:m) {
      xt <- x[i-1]
      y <- rchisq(1, df = xt)
      r <- fr(xt, y) * dchisq(xt, df=y) / dchisq(y, df=xt)
      if (u[i] <= r) x[i] <- y else {
        x[i] <- xt
        k <- k+1     #y is rejected
      }
    }

    k
    plot(acf(x))

### Example 11.4 (Expected lifetime M-H sampler, continued)

    # MCMC Metropolis-Hastings with gamma proposal
    m <- 10000
    x <- numeric(m)
    a <- 4
    x[1] <- rlnorm(1)  #initialize chain
    k <- 0
    u <- runif(m)
    for (i in 2:m) {
      xt <- x[i-1]
      y <- rgamma(1, shape=a, rate=a/xt)
      r <- fr(xt, y) * dgamma(xt, shape=a, rate=a/y) / 
        dgamma(y, shape=a, rate=a/xt)
      if (u[i] <= r) x[i] <- y else {
        x[i] <- xt
        k <- k+1     #y is rejected
      }
    }

    k / m    #proportion rejected
    plot(acf(x))

### Example 11.5 (Learning about the distribution of mean future lifetime)

    burnin <- m/2
    X <- x[-(1:burnin)]
    hist(X, prob=TRUE, breaks="scott", xlab=bquote(psi),
         main="MCMC replicates using gamma proposal") -> h
    curve(f.mu(x), add=TRUE)
    i <- which.max(h$density)  #estimating the mode
    h$mids[i]

    q <- quantile(X, c(.025, .975), type=1)
    round(q, 1)     #hours

    HPDi <- function(x, prob = 0.95) {
      ## HPD interval for a single MCMC chain (a vector)
      x <- sort(x)
      n <- length(x)
      m <- floor(prob * n)
      i <- 1:(n - m)
      L <- x[i + m] - x[i]
      best <- which.min(L)
      return (c(lower=x[best], upper=x[best+m]))
    }
    
    HPDi(X) 

### Example 11.6 (Random walk Metropolis) 

    rw.Metropolis <- function(n, sigma, x0, N) {
        x <- numeric(N)
        x[1] <- x0
        u <- runif(N)
        k <- 0
        for (i in 2:N) {
            y <- rnorm(1, x[i-1], sigma)
                if (u[i] <= (dt(y, n) / dt(x[i-1], n)))
                x[i] <- y  else {
                    x[i] <- x[i-1]
                    k <- k + 1
                }
            }
        return(list(x=x, k=k))
        }

    n <- 4  #degrees of freedom for target Student t dist.
    N <- 2000
    sigma <- c(.05, .5, 2,  16)

    x0 <- 25
    rw1 <- rw.Metropolis(n, sigma[1], x0, N)
    rw2 <- rw.Metropolis(n, sigma[2], x0, N)
    rw3 <- rw.Metropolis(n, sigma[3], x0, N)
    rw4 <- rw.Metropolis(n, sigma[4], x0, N)

    #number of candidate points rejected
    print(c(rw1$k, rw2$k, rw3$k, rw4$k))
    

### Code for Figure 11.6

    par(mfrow=c(2,2))  #display 4 graphs together
    refline <- qt(c(.025, .975), df=n)
    rw <- cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
    for (j in 1:4) {
        plot(rw[,j], type="l",
             xlab=bquote(sigma == .(round(sigma[j],3))),
             ylab="X", ylim=range(rw[,j]))
        abline(h=refline)
    }
    par(mfrow=c(1,1)) #reset to default


### Example 11.7 (Example 11.6, cont.)

    a <- c(.05, seq(.1, .9, .1), .95)
    Q <- qt(a, n)
    rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)
    mc <- rw[501:N, ]
    Qrw <- apply(mc, 2, function(x) quantile(x, a))
    print(round(cbind(Q, Qrw), 3))
    xtable::xtable(round(cbind(Q, Qrw), 3)) #latex format

### Code for Figures 11.7(a) and 11.7(b)

    plot(x, type="l")
    abline(h=b, v=501, lty=3)
    xb <- x[- (1:501)]
    hist(xb, prob=TRUE, xlab=bquote(beta), ylab="X", main="")
    z <- seq(min(xb), max(xb), length=100)
    lines(z, dnorm(z, mean(xb), sd(xb)))


### Example 11.8 (Bayesian inference: A simple investment model)

    b <- .2          #actual value of beta
    w <- .25         #width of the uniform support set
    m <- 5000        #length of the chain
    burn <- 1000     #burn-in time
    days <- 250
    x <- numeric(m)  #the chain

    # generate the observed frequencies of winners
    i <- sample(1:5, size=days, replace=TRUE,
            prob=c(1, 1-b, 1-2*b, 2*b, b))
    win <- tabulate(i)
    print(win)

    prob <- function(y, win) {
        # computes (without the constant) the target density
        if (y < 0 || y >= 0.5)
            return (0)
        return((1/3)^win[1] *
            ((1-y)/3)^win[2] * ((1-2*y)/3)^win[3] *
                ((2*y)/3)^win[4] * (y/3)^win[5])
    }

    u <- runif(m)         #for accept/reject step
    v <- runif(m, -w, w)  #proposal distribution
    x[1] <- .25
    for (i in 2:m) {
        y <- x[i-1] + v[i]
        if (u[i] <= prob(y, win) / prob(x[i-1], win))
            x[i] <- y  else
                x[i] <- x[i-1]
    }

    print(win)
    print(round(win/days, 3))
    print(round(c(1, 1-b, 1-2*b, 2*b, b)/3, 3))
    xb <- x[(burn+1):m]
    print(mean(xb))
    

### Example 11.9 (Independence sampler)

    m <- 5000 #length of chain
    xt <- numeric(m)
    a <- 1             #parameter of Beta(a,b) proposal dist.
    b <- 1             #parameter of Beta(a,b) proposal dist.
    p <- .2            #mixing parameter
    n <- 30            #sample size
    mu <- c(0, 5)      #parameters of the normal densities
    sigma <- c(1, 1)

    # generate the observed sample
    i <- sample(1:2, size=n, replace=TRUE, prob=c(p, 1-p))
    x <- rnorm(n, mu[i], sigma[i])

    # generate the independence sampler chain
    u <- runif(m)
    y <- rbeta(m, a, b)      #proposal distribution
    xt[1] <- .5

    for (i in 2:m) {
        fy <- y[i] * dnorm(x, mu[1], sigma[1]) +
                (1-y[i]) * dnorm(x, mu[2], sigma[2])
        fx <- xt[i-1] * dnorm(x, mu[1], sigma[1]) +
                (1-xt[i-1]) * dnorm(x, mu[2], sigma[2])

        r <- prod(fy / fx) *
               (xt[i-1]^(a-1) * (1-xt[i-1])^(b-1)) /
                 (y[i]^(a-1) * (1-y[i])^(b-1))

        if (u[i] <= r) xt[i] <- y[i] else
            xt[i] <- xt[i-1]
        }

    plot(xt, type="l", ylab="p")
    hist(xt[101:m], main="", xlab="p", prob=TRUE)
    print(mean(xt[101:m]))
    

### Example 11.10 (Gibbs sampler: Bivariate distribution)

    #initialize constants and parameters
    N <- 5000               #length of chain
    burn <- 1000            #burn-in length
    X <- matrix(0, N, 2)    #the chain, a bivariate sample

    rho <- -.75             #correlation
    mu1 <- 0
    mu2 <- 2
    sigma1 <- 1
    sigma2 <- .5
    s1 <- sqrt(1-rho^2)*sigma1
    s2 <- sqrt(1-rho^2)*sigma2

    ###### generate the chain #####

    X[1, ] <- c(mu1, mu2)            #initialize

    for (i in 2:N) {
        x2 <- X[i-1, 2]
        m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
        X[i, 1] <- rnorm(1, m1, s1)
        x1 <- X[i, 1]
        m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
        X[i, 2] <- rnorm(1, m2, s2)
    }

    b <- burn + 1
    x <- X[b:N, ]

    # compare sample statistics to parameters
    colMeans(x)
    cov(x)
    cor(x)

    plot(x, main="", cex=.5, xlab=bquote(X[1]),
         ylab=bquote(X[2]), ylim=range(x[,2]))
         
### Example 11.11 (Gelman-Rubin method of monitoring convergence)

    Gelman.Rubin <- function(psi) {
        # psi[i,j] is the statistic psi(X[i,1:j])
        # for chain in i-th row of X
        psi <- as.matrix(psi)
        n <- ncol(psi)
        k <- nrow(psi)

        psi.means <- rowMeans(psi)     #row means
        B <- n * var(psi.means)        #between variance est.
        psi.w <- apply(psi, 1, "var")  #within variances
        W <- mean(psi.w)               #within est.
        v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
        r.hat <- v.hat / W             #G-R statistic
        return(r.hat)
        }

    normal.chain <- function(sigma, N, X1) {
        #generates a Metropolis chain for Normal(0,1)
        #with Normal(X[t], sigma) proposal distribution
        #and starting value X1
        x <- rep(0, N)
        x[1] <- X1
        u <- runif(N)

        for (i in 2:N) {
            xt <- x[i-1]
            y <- rnorm(1, xt, sigma)     #candidate point
            r1 <- dnorm(y, 0, 1) * dnorm(xt, y, sigma)
            r2 <- dnorm(xt, 0, 1) * dnorm(y, xt, sigma)
            r <- r1 / r2
            if (u[i] <= r) x[i] <- y else
                 x[i] <- xt
            }
        return(x)
        }

    sigma <- .2     #parameter of proposal distribution
    k <- 4          #number of chains to generate
    n <- 15000      #length of chains
    b <- 1000       #burn-in length

    #choose overdispersed initial values
    x0 <- c(-10, -5, 5, 10)

    #generate the chains
    X <- matrix(0, nrow=k, ncol=n)
    for (i in 1:k)
        X[i, ] <- normal.chain(sigma, n, x0[i])

    #compute diagnostic statistics
    psi <- t(apply(X, 1, cumsum))
    for (i in 1:nrow(psi))
        psi[i,] <- psi[i,] / (1:ncol(psi))
    print(Gelman.Rubin(psi))

    #plot psi for the four chains
    par(mfrow=c(2,2))
    for (i in 1:k)
        plot(psi[i, (b+1):n], type="l",
            xlab=i, ylab=bquote(psi))
    par(mfrow=c(1,1)) #restore default

    #plot the sequence of R-hat statistics
    rhat <- rep(0, n)
    for (j in (b+1):n)
        rhat[j] <- Gelman.Rubin(psi[,1:j])
    plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
    abline(h=1.1, lty=2)


### Example 11.12 (Coal mining disasters)

    library(boot)     #for coal data
    data(coal)
    year <- floor(coal)
    y <- table(year)
    plot(y)  #a time plot

    y <- floor(coal[[1]])
    y <- tabulate(y)
    y <- y[1851:length(y)]

    # Gibbs sampler for the coal mining change point

    # initialization
    n <- length(y)    #length of the data
    m <- 1000         #length of the chain
    mu <- lambda <- k <- numeric(m)
    L <- numeric(n)
    k[1] <- sample(1:n, 1)
    mu[1] <- 1
    lambda[1] <- 1
    b1 <- 1
    b2 <- 1

    # run the Gibbs sampler
    for (i in 2:m) {
        kt <- k[i-1]

        #generate mu
        r <- .5 + sum(y[1:kt])
        mu[i] <- rgamma(1, shape = r, rate = kt + b1)

        #generate lambda
        if (kt + 1 > n) r <- .5 + sum(y) else
            r <- .5 + sum(y[(kt+1):n])
        lambda[i] <- rgamma(1, shape = r, rate = n - kt + b2)

        #generate b1 and b2
        b1 <- rgamma(1, shape = .5, rate = mu[i]+1)
        b2 <- rgamma(1, shape = .5, rate = lambda[i]+1)

        for (j in 1:n) {
            L[j] <- exp((lambda[i] - mu[i]) * j) *
                          (mu[i] / lambda[i])^sum(y[1:j])
            }
        L <- L / sum(L)

        #generate k from discrete distribution L on 1:n
        k[i] <- sample(1:n, prob=L, size=1)
    }

    b <- 201
    j <- k[b:m]
    print(mean(k[b:m]))
    print(mean(lambda[b:m]))
    print(mean(mu[b:m]))


### Code for Figure 11.14

    # plots of the chains for Gibbs sampler output

    par(mfcol=c(3,1), ask=TRUE)
    plot(mu, type="l", ylab="mu")
    plot(lambda, type="l", ylab="lambda")
    plot(k, type="l", ylab="change point = k")


### Code for Figure 11.15

    # histograms from the Gibbs sampler output

    par(mfrow=c(2,3))
    labelk <- "changepoint"
    label1 <- paste("mu", round(mean(mu[b:m]), 1))
    label2 <- paste("lambda", round(mean(lambda[b:m]), 1))

    hist(mu[b:m], main="", xlab=label1,
         breaks = "scott", prob=TRUE) #mu posterior
    hist(lambda[b:m], main="", xlab=label2,
         breaks = "scott", prob=TRUE) #lambda posterior
    hist(j, breaks=min(j):max(j), prob=TRUE, main="",
        xlab = labelk)
    par(mfcol=c(1,1), ask=FALSE)  #restore display

    