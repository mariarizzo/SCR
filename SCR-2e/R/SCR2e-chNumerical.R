#######################################################
###       Statistical Computing with R  2e          ###
###       Maria L. Rizzo                            ###
###       Chapman & Hall/CRC The R Series           ###
###       ISBN 9781466553323 - CAT# K15269          ###
###       January 2019                              ###
###                                                 ###
###       R code for Chapter 13                     ###
###       Introduction to Numerical Methods in R    ###
#######################################################

# packages to install: gsl

### Example 13.1 (Identical and nearly equal)

    isTRUE(all.equal(.2, .3 - .1))
    all.equal(.2, .3)          #not a logical value
    isTRUE(all.equal(.2, .3))  #always a logical value

    x <- 1:4
    y <- 2
    y == 2
    x == y  #not necessarily a single logical value
    identical(x, y)  #always a single logical value
    identical(y, 2)


### Example 13.2 (Ratio of two large numbers)

    n <- 400
    (gamma((n-1)/2) / (sqrt(pi) * gamma((n-2)/2)))
    exp(lgamma((n-1)/2) - lgamma((n-2)/2)) / sqrt(pi)


### Example 13.3 (Taylor expansion)

    system.time({
        for (i in 1:1000) {
            a <- rep(0, 24)
            a0 <- pi / 6
            a2 <- a0 * a0
            a[1] <- -a0^3 / 6
            for (i in 2:24)
                a[i] <- - a2 * a[i-1] / ((2*i+1)*(2*i))
            a0 + sum(a)}
        })

    system.time({
        for (i in 1:1000) {
            K <- 2 * (0:24) + 1
            i <- rep(c(1, -1), length=25)
            sum(i * (pi/6)^K / factorial(K))}
        })


### Example 13.4 (Derivative of zeta function)

    zeta.deriv <- function(a) {
        z <- a - 1
        # Stieltjes constants gamma_k for k=1:5
        g <- c(
            -.7281584548367672e-1,
            -.9690363192872318e-2,
             .2053834420303346e-2,
             .2325370065467300e-2,
             .7933238173010627e-3)
        i <- c(-1, 1, -1, 1, -1)
        n <- 0:4
        -1/z^2 + sum(i * g * z^n / factorial(n))
    }


### Example 13.5 (Derivative of zeta function, cont.)

    library(gsl)   #for zeta function
    z <- c(1.001, 1.01, 1.5, 2, 3, 5)
    h <- .Machine$double.eps^0.5
    dz <- dq <- rep(0, length(z))
    for (i in 1:length(z)) {
        v <- z[i] + h
        h <- v - z[i]
        a0 <- z[i] - h
        if (a0 < 1) a0 <- (1 + z[i])/2
        a1 <- z[i] + h
        dq[i] <- (zeta(a1) - zeta(a0)) / (a1 - a0)
        dz[i] <- zeta.deriv(z[i])
    }

    h

    cbind(z, dz, dq)

    
### Example 13.6 (Solving f(x)=0)

    f <- function(y, a, n) {
        a^2 + y^2 + 2*a*y/(n-1) - (n-2)
    }

    a <- 0.5
    n <- 20
    b0 <- 0
    b1 <- 5*n

    #solve using bisection
    it <- 0
    eps <- .Machine$double.eps^0.25
    r <- seq(b0, b1, length=3)
    y <- c(f(r[1], a, n), f(r[2], a, n), f(r[3], a, n))
    if (y[1] * y[3] > 0)
        stop("f does not have opposite sign at endpoints")

    while(it < 1000 && abs(y[2]) > eps) {
        it <- it + 1
        if (y[1]*y[2] < 0) {
            r[3] <- r[2]
            y[3] <- y[2]
        } else {
            r[1] <- r[2]
            y[1] <- y[2]
        }
        r[2] <- (r[1] + r[3]) / 2
        y[2] <- f(r[2], a=a, n=n)
        print(c(r[1], y[1], y[3]-y[2]))
    }


### Example 13.7 (Solving f(x)=0 with Brent's method: uniroot)

    a <- 0.5
    n <- 20
    out <- uniroot(function(y) {
               a^2 + y^2 + 2*a*y/(n-1) - (n-2) },
               lower = 0, upper = n*5)
    unlist(out)
    uniroot(function(y) {a^2 + y^2 + 2*a*y/(n-1) - (n-2)},
            interval = c(-n*5, 0))$root


### Example 13.8 (Numerical integration with integrate)

    f <- function(y, N, r, rho) {
        (cosh(y) - rho * r)^(1 - N)
    }
    integrate(f, lower=0, upper=Inf,
              rel.tol=.Machine$double.eps^0.25,
              N=10, r=0.5, rho=0.2)

    ro <- seq(-.99, .99, .01)
    v <- rep(0, length(ro))
    for (i in 1:length(ro)) {
        v[i] <- integrate(f, lower=0, upper=Inf,
                  rel.tol=.Machine$double.eps^0.25,
                  N=10, r=0.5, rho=ro[i])$value
        }
    plot(ro, v, type="l", xlab=expression(rho),
         ylab="Integral Value (n=10, r=0.5)")


### Example 13.9 (Density of sample correlation coefficient)

    .dcorr <- function(r, N, rho=0) {
        # compute the density function of sample correlation
        if (abs(r) > 1 || abs(rho) > 1) return (0)
        if (N < 4) return (NA)

        if (isTRUE(all.equal(rho, 0.0))) {
            a <- exp(lgamma((N - 1)/2) - lgamma((N - 2)/2)) /
                     sqrt(pi)
            return (a * (1 - r^2)^((N - 4)/2))
        }

        # if rho not 0, need to integrate
        f <- function(w, R, N, rho)
            (cosh(w) - rho * R)^(1 - N)

        #need to insert some error checking here
        i <- integrate(f, lower=0, upper=Inf,
                R=r, N=N, rho=rho)$value
        c1 <- (N - 2) * (1 - rho^2)^((N - 1)/2)
        c2 <- (1 - r^2)^((N - 4) / 2) / pi
        return(c1 * c2 * i)
    }

    r <- as.matrix(seq(-1, 1, .01))
    d1 <- apply(r, 1, .dcorr, N=10, rho=.0)
    d2 <- apply(r, 1, .dcorr, N=10, rho=.5)
    d3 <- apply(r, 1, .dcorr, N=10, rho=-.5)
    plot(r, d2, type="l", lty=2, lwd=2, ylab="density")
    lines(r, d1, lwd=2)
    lines(r, d3, lty=4, lwd=2)
    legend("top", inset=.02,
           c("rho = 0", "rho = 0.5", "rho = -0.5"), lty=c(1,2,4), lwd=2)


### Example 13.10 (MLE using mle)

    #the observed sample
    y <- c(0.04304550, 0.50263474)

    mlogL <- function(theta=1) {
        #minus log-likelihood of exp. density, rate theta
        return( - (length(y) * log(theta) - theta * sum(y)))
    }

    library(stats4)
    fit <- mle(mlogL)
    summary(fit)

    # Alternately, the initial value for the optimizer could
    # be supplied in the call to mle; two examples are

    mle(mlogL, start=list(theta=1))
    mle(mlogL, start=list(theta=mean(y)))

### Application: Evaluating an Expected Value

    y <- VGAM::rpareto(1000, scale=1, shape=3)
    hist(y, prob=TRUE, breaks="scott", main="", ylim=c(0,3))
    curve(VGAM::dpareto(x, scale=1, shape=3), add=TRUE)

    a <- 1
    s <- 1
    b <- 0.5
    
### Example 13.11 (Numerical integration)

    f1 <- function(x, y, s, a, b) {
      # the integrand function: |y-x|^b f(x)
      (abs(y-x))^b * a * s^a / x^(a+1)
    }
    
    integrate(f1, lower=s, upper=Inf, y=2, s=1, a=1, b=b)
    
### Example 13.12 (Direct evaluation)

    g2 <- function(y, s, b)   {
      # Compute E|y-X|^b for Pareto(I) with a=1
      y0 <- (y - s)/y
      c0 <-  gsl::hyperg_2F1(b+1, 1, b+2, y0)
      ((y - s)^b) - s*b*(y^(b-1)) *
      ((y0^b)/b + c0*(y0^(b+1))/(b+1)) +
      s*y^(b-1) * beta(b+1,1-b)
      }

    g2(2, s, b)

### Example 13.13 (Plot of expected distance function)

    p <- c(ppoints(50)/100, ppoints(50))
    x <- VGAM::qpareto(p, scale=s, shape=1)
    ex <- g2(x, s=s, b=b)
    plot(log(x), log(ex), cex=.25, type="l")
    for (i in 1:length(x)) {
      y <- x[i]
      zi <- integrate(f1, lower=s, upper=Inf,
      subdivisions=200, rel.tol=.Machine$double.eps^.5,
      stop.on.error=FALSE, y=y, s=s, a=1, b=b)
      if (zi$message == "OK")
      points(log(y), log(zi$value), col=2, cex=.25)
      else print(paste("y=", y, zi$Message))
    }

    integrate(f1, lower=s, upper=Inf, y=s, s=s, a=1, b=b)
    g2(s, s, b)
