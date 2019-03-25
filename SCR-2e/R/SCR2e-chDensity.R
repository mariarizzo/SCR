#######################################################
###       Statistical Computing with R  2e          ###
###       Maria L. Rizzo                            ###
###       Chapman & Hall/CRC The R Series           ###
###       ISBN 9781466553323 - CAT# K15269          ###
###       March 6, 2019                             ###
###                                                 ###
###       R code for Chapter 10                     ###
###       Permutation Tests                         ###
#######################################################


# packages to install: ash, ggplot2 

### Example 12.1 (Histogram density estimates using Sturges' Rule)

    set.seed(12345) 
    n <- 25
    x <- rnorm(n)
    # calc breaks according to Sturges' Rule
    nclass <- ceiling(1 + log2(n))
    cwidth <- diff(range(x) / nclass)
    breaks <- min(x) + cwidth * 0:nclass
    h.default <- hist(x, freq = FALSE, xlab = "default",
        main = "hist: default")
    z <- qnorm(ppoints(1000))
    lines(z, dnorm(z))
    h.sturges <- hist(x, breaks = breaks, freq = FALSE,
        main = "hist: Sturges")
    lines(z, dnorm(z))

    print(h.default$breaks)
    print(h.default$counts)
    print(round(h.sturges$breaks, 1))
    print(h.sturges$counts)
    print(cwidth)

    print(h.default$density[5])
    print(h.sturges$density[4])


### Example 12.2 (Density estimates from a histogram)

    #using histograms from Example 12.1
    x0 <- .1
    b <- which.min(h.default$breaks <= x0) - 1
    print(c(b, h.default$density[b]))
    b <- which.min(h.sturges$breaks <= x0) - 1
    print(c(b, h.sturges$density[b]))

    h.default$counts[7] / (n * 0.5)
    h.sturges$counts[6] / (n * cwidth)


### Example 12.3 (Density estimation for Old Faithful)

    library(MASS)  #for geyser and truehist
    waiting <- geyser$waiting
    n <- length(waiting)
    # rounding the constant in Scott's rule
    # and using sample standard deviation to estimate sigma
    h <- 3.5 * sd(waiting) * n^(-1/3)

    # number of classes is determined by the range and h
    m <- min(waiting)
    M <- max(waiting)
    nclass <- ceiling((M - m) / h)
    breaks <- m + h * 0:nclass

    par(ask = TRUE)  #prompt to see next graph
    h.scott <- hist(waiting, breaks = breaks, freq = FALSE,
        main = "")
    truehist(waiting, nbins = "Scott", x0 = 0, prob=TRUE, col = 0)
    hist(waiting, breaks = "scott", prob=TRUE, density=5,
        add=TRUE)
    par(ask = FALSE)


### Example 12.4 (Frequency polygon density estimate)

    waiting <- geyser$waiting   #in MASS
    n <- length(waiting)
    # freq poly bin width using normal ref rule
    h <- 2.15 * sqrt(var(waiting)) * n^(-1/5)

    # calculate the sequence of breaks and histogram
    br <- pretty(waiting, diff(range(waiting)) / h)
    brplus <- c(min(br)-h, max(br+h))
    histg <- hist(waiting, breaks = br, freq = FALSE,
        main = "", xlim = brplus)

    vx <- histg$mids     #density est at vertices of polygon
    vy <- histg$density
    delta <- diff(vx)[1] # h after pretty is applied
    k <- length(vx)
    vx <- vx + delta     # the bins on the ends
    vx <- c(vx[1] - 2 * delta, vx[1] - delta, vx)
    vy <- c(0, vy, 0)

    # add the polygon to the histogram
    polygon(vx, vy)

    # check estimates by numerical integration
    fpoly <- approxfun(vx, vy)
    print(integrate(fpoly, lower=min(vx), upper=max(vx)))
    
    library(ggplot2)
    ggplot(geyser, aes(waiting)) + geom_frepoly(binsize=h)
    
### Example 12.6 (ASH density estimate)

    library(MASS)
    waiting <- geyser$waiting
    n <- length(waiting)
    m <- 20
    a <- min(waiting) - .5
    b <- max(waiting) + .5
    h <- 7.27037
    delta <- h / m

    #get the bin counts on the delta-width mesh.
    br <- seq(a - delta*m, b + 2*delta*m, delta)
    histg <- hist(waiting, breaks = br, plot = FALSE)
    nk <- histg$counts
    K <- abs((1-m):(m-1))

    fhat <- function(x) {
        # locate the leftmost interval containing x
        i <- max(which(x > br))
        k <- (i - m + 1):(i + m - 1)
        # get the 2m-1 bin counts centered at x
        vk <- nk[k]
        sum((1 - K / m) * vk) / (n * h)   #f.hat
        }

    # density can be computed at any points in range of data
    z <- as.matrix(seq(a, b + h, .1))
    f.ash <- apply(z, 1, fhat)   #density estimates at midpts

    # plot ASH density estimate over histogram
    br2 <- seq(a, b + h, h)
    hist(waiting, breaks = br2, freq = FALSE, main = "",
        ylim = c(0, max(f.ash)))
    lines(z, f.ash, xlab = "waiting")
    
    
### Example 12.7 (KDE of Old Faithful waiting time)

    library(MASS)
    waiting <- geyser$waiting
    n <- length(waiting)

    h1 <- 1.06 * sd(waiting) * n^(-1/5)
    h2 <- .9 * min(c(IQR(waiting)/1.34, sd(waiting))) * n^(-1/5)
    plot(density(waiting))

    print(density(waiting))

    sdK <- density(kernel = "gaussian", give.Rkern = TRUE)
    print(c(sdK, sdK * sd(waiting)))
    print(c(sd(waiting), IQR(waiting)))
    print(c(h1, h2))


### Example 12.8 (KDE of precipitation data)

    n <- length(precip)
    h1 <- 1.06 * sd(precip) * n^(-1/5)
    h2 <- .9 * min(c(IQR(precip)/1.34, sd(precip))) * n^(-1/5)
    h0 <- bw.nrd0(precip)

    par(mfrow = c(2, 2))
    plot(density(precip))           #default Gaussian (h0)
    plot(density(precip, bw = h1))  #Gaussian, bandwidth h1
    plot(density(precip, bw = h2))  #Gaussian, bandwidth h2
    plot(density(precip, kernel = "cosine"))
    par(mfrow = c(1,1))

    print(c(h0, h1, h2))
    
    
### Example 12.9 (Computing $\hat f(x)$ for arbitrary x)

    d <- density(precip)
    xnew <- seq(0, 70, 10)
    approx(d$x, d$y, xout = xnew)

    fhat <- approxfun(d$x, d$y)
    fhat(xnew)    


### Example 12.10 (Exponential density)

    x <- rexp(1000, 1)
    plot(density(x), xlim = c(-1, 6), ylim = c(0, 1), main="")
    abline(v = 0)

    # add the true density to compare
    y <- seq(.001, 6, .01)
    lines(y, dexp(y, 1), lty = 2)
    

### Example 12.11 (Reflection boundary technique)

    xx <- c(x, -x)
    g <- density(xx, bw = bw.nrd0(x))
    a <- seq(0, 6, .01)

    ghat <- approx(g$x, g$y, xout = a)
    fhat <- 2 * ghat$y       # density estimate along a

    bw <- paste("Bandwidth = ", round(g$bw, 5))
    plot(a, fhat, type="l", xlim=c(-1, 6), ylim=c(0, 1),
        main = "", xlab = bw, ylab = "Density")
    abline(v = 0)

    # add the true density to compare
    y <- seq(.001, 6, .01)
    lines(y, dexp(y, 1), lty = 2)
    

### Example 12.12 (Bivariate frequency table: bin2d)

    bin2d <-
      function(x, breaks1 = "Sturges", breaks2 = "Sturges"){
      # Data matrix x is n by 2
      # breaks1, breaks2: any valid breaks for hist function
      # using same defaults as hist
      histg1 <- hist(x[,1], breaks = breaks1, plot = FALSE)
      histg2 <- hist(x[,2], breaks = breaks2, plot = FALSE)
      brx <- histg1$breaks
      bry <- histg2$breaks

      # bin frequencies
      freq <- table(cut(x[,1], brx),  cut(x[,2], bry))

      return(list(call = match.call(), freq = freq,
                breaks1 = brx, breaks2 = bry,
                mids1 = histg1$mids, mids2 = histg2$mids))
      }
    
    bin2d(iris[1:50,1:2])


### Example 12.13 (Bivariate density polygon)

    #generate standard bivariate normal random sample
    n <- 2000;   d <- 2
    x <- matrix(rnorm(n*d), n, d)

    # compute the frequency table and density estimates
    # using bin2d function from the previous example
    b <- bin2d(x)
    h1 <- diff(b$breaks1)
    h2 <- diff(b$breaks2)

    # matrix h contains the areas of the bins in b
    h <- outer(h1, h2, "*")

    Z <- b$freq / (n * h)  # the density estimate

    persp(x=b$mids1, y=b$mids2, z=Z, shade=TRUE,
          xlab="X", ylab="Y", main="",
          theta=45, phi=30, ltheta=60)


### Example 12.14 (Bivariate ASH density estimate)

    library(ash)  # for bivariate ASH density est.
    # generate N_2(0,Sigma) data
    n <- 2000
    d <- 2
    nbin <- c(30, 30)          # number of bins
    m <- c(5, 5)               # smoothing parameters

    # First example with positive correlation
    Sigma <- matrix(c(1, .9, .9, 1), 2, 2)
    set.seed(345)
    
    #rmvn.eigen from Chapter 3 used to generate data
    #alternately mvrnorm (MASS) can be used here
    
    x <- rmvn.eigen(n, c(0, 0), Sigma=Sigma)
    #x <- MASS::mvrnorm(n, c(0, 0), Sigma)
    b <- bin2(x, nbin = nbin)
    # kopt is the kernel type, here triangular
    est <- ash2(b, m = m, kopt = c(1,0))

    persp(x = est$x, y = est$y, z = est$z, shade=TRUE,
          xlab = "X", ylab = "Y", zlab = "", main="",
          theta = 30, phi = 75, ltheta = 30, box = FALSE)
    contour(x = est$x, y = est$y, z = est$z, main="")

    # Second example with negative correlation
    Sigma <- matrix(c(1, -.9, -.9, 1), 2, 2)
    set.seed(345)
#   x <- rmvn.eigen(n, c(0, 0), Sigma=Sigma)  #source from ch 3 or use mvrnorm
    x <- MASS::mvrnorm(n, c(0, 0), Sigma)
    b <- bin2(x, nbin = nbin)
    est <- ash2(b, m = m, kopt = c(1,0))

    persp(x = est$x, y = est$y, z = est$z, shade=TRUE,
          xlab = "X", ylab = "Y", zlab = "", main="",
          theta = 30, phi = 75, ltheta = 30, box = FALSE)
    contour(x = est$x, y = est$y, z = est$z, main="")
    par(ask = FALSE)
    detach(package:ash)
    
        
### Example 12.15 (Product kernel estimate of a BVN mixture)

    library(MASS)  #for mvrnorm and kde2d
    #generate the normal mixture data
    n <- 2000
    p <- c(.2, .3, .5)
    mu <- matrix(c(0, 1, 4, 0, 3, -1), 3, 2)
    Sigma <- diag(2)
    i <- sample(1:3, replace = TRUE, prob = p, size = n)
    k <- table(i)

    x1 <- mvrnorm(k[1], mu = mu[1,], Sigma)
    x2 <- mvrnorm(k[2], mu = mu[2,], Sigma)
    x3 <- mvrnorm(k[3], mu = mu[3,], Sigma)
    X <-  rbind(x1, x2, x3)   #the mixture data
    x <- X[,1]
    y <- X[,2]
    
    print(c(bandwidth.nrd(x), bandwidth.nrd(y)))

    # accepting the default normal reference bandwidth
    fhat <- kde2d(x, y)
    contour(fhat)
    persp(fhat, phi = 30, theta = 20, d = 5, xlab = "x")

    # select bandwidth by unbiased cross-validation
    h = c(ucv(x), ucv(y))
    fhat <- kde2d(x, y, h = h)
    contour(fhat)
    persp(fhat, phi = 30, theta = 20, d = 5, xlab = "x")



### Code to generate data as shown in Table 10.1

    N <- c(10, 20, 30, 50, 100, 200, 500, 1000, 5000, 10000)
    m <- length(N)
    out <- matrix(0, nrow = m, ncol = 8)
    out[ ,1] <- N
    out[ ,5] <- N
    for (i in 1:m) {
        x <- rnorm(N[i])
        out[i, 2:4] <- c(nclass.Sturges(x),
            nclass.scott(x), nclass.FD(x))
        x <- rexp(N[i])
        out[i, 6:8] <- c(nclass.Sturges(x),
            nclass.scott(x), nclass.FD(x))
    }
    print(out)


### Code to plot the histograms in Figure 12.4

    library(MASS)  #for truehist
    par(mfrow = c(2, 2))
    x <- sort(rnorm(1000))
    y <- dnorm(x)
    o <- (1:4) / 4
    h <- .35
    for (i in 1:4) {
        truehist(x, prob = TRUE, h = .35, x0 = o[i],
            xlim = c(-3.5, 3.5), ylim = c(0, 0.45),
            ylab = "Density", main = "")
        lines(x, y)
    }
    par(mfrow = c(1, 1))


### Code to plot Figure 12.6

    #set.seed(7555)
    n <- 10
    y <- rnorm(n)

    par(mfrow = c(2, 2))
    for (h in c(.25, .4, .6, 1)) {
        x <- seq(-4, 4, .01)
        fhat <- rep(0, length(x))
        # set up the plot window first
        plot(x, fhat, type="n", xlab="", ylab="",
            main=paste("h=",h), xlim=c(-4,4), ylim=c(0, .5))
        for (i in 1:n) {
            # plot a normal density at each sample pt
            z <- (x - y[i]) / h
            f <- dnorm(z)
            lines(x, f / (n * h))
            # sum the densities to get the estimates
            fhat <- fhat + f / (n * h)
        }
        lines(x, fhat, lwd=2) # add density estimate to plot
    }

    par(mfrow = c(1, 1))


### Code to plot kernels in Figure 12.7

    #see examples for density, kernels in S parametrization
    (kernels <- eval(formals(density.default)$kernel))

    plot(density(0, from=-1.2, to=1.2, width=2,
        kern="gaussian"), type="l", ylim=c(0, 1),
        xlab="", main="")
    for(i in 2:5)
        lines(density(0, width=2, kern=kernels[i]), lty=i)
    legend("topright", legend=kernels[1:5],
        lty=1:5, inset=.02)
