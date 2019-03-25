#######################################################
###       Statistical Computing with R  2e          ###
###       Maria L. Rizzo                            ###
###       Chapman & Hall/CRC The R Series           ###
###       ISBN 9781466553323 - CAT# K15269          ###
###       March 6, 2019                             ###
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
        #minus log-likelihood of exp. density, rate 1/theta
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
      zi <- integrate(f, lower=s, upper=Inf,
      subdivisions=200, rel.tol=.Machine$double.eps^.5,
      stop.on.error=FALSE, y=y, s=s, a=1, b=b)
      if (zi$message == "OK")
      points(log(y), log(zi$value), col=2, cex=.25)
      else print(paste("y=", y, zi$Message))
    }

    integrate(f1, lower=s, upper=Inf, y=s, s=s, a=1, b=b)
    g2(s, s, b)

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
