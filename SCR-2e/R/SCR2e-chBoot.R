    #######################################################
    ###       Statistical Computing with R  2e          ###
    ###       Maria L. Rizzo                            ###
    ###       Chapman & Hall/CRC The R Series           ###
    ###       ISBN 9781466553323 - CAT# K15269          ###
    ###       March 6, 2019                             ###
    ###                                                 ###
    ###       R code for Chapter 8                      ###
    ###       Bootstrap and Jackknife                   ###
    #######################################################

    
# packages to install: bootstrap, DAAG

    
### Example 8.2 (Bootstrap estimate of standard error)

    library(bootstrap)    #for the law data
    print(cor(law$LSAT, law$GPA))
    print(cor(law82$LSAT, law82$GPA))

    #set up the bootstrap
    B <- 200            #number of replicates
    n <- nrow(law)      #sample size
    R <- numeric(B)     #storage for replicates

    #bootstrap estimate of standard error of R
    for (b in 1:B) {
        #randomly select the indices
        i <- sample(1:n, size = n, replace = TRUE)
        LSAT <- law$LSAT[i]       #i is a vector of indices
        GPA <- law$GPA[i]
        R[b] <- cor(LSAT, GPA)
    }
    #output
    print(se.R <- sd(R))
    hist(R, prob = TRUE)


### Example 8.3 (Bootstrap estimate of standard error: boot function)

    r <- function(x, i) {
        #want correlation of columns 1 and 2
        cor(x[i,1], x[i,2])
    }

    library(boot)       #for boot function
    obj <- boot(data = law, statistic = r, R = 2000)
    obj
    y <- obj$t
    sd(y)


### Example 8.4 (Bootstrap estimate of bias)

    #sample estimate for n=15
    theta.hat <- cor(law$LSAT, law$GPA)

    #bootstrap estimate of bias
    B <- 2000   #larger for estimating bias
    n <- nrow(law)
    theta.b <- numeric(B)

    for (b in 1:B) {
        i <- sample(1:n, size = n, replace = TRUE)
        LSAT <- law$LSAT[i]
        GPA <- law$GPA[i]
        theta.b[b] <- cor(LSAT, GPA)
    }
    bias <- mean(theta.b - theta.hat)
    bias

    
### Example 8.5 (Bootstrap estimate of bias of a ratio estimate)


    data(patch, package = "bootstrap")
    patch

    n <- nrow(patch)  #in bootstrap package
    B <- 2000
    theta.b <- numeric(B)
    theta.hat <- mean(patch$y) / mean(patch$z)

    #bootstrap
    for (b in 1:B) {
        i <- sample(1:n, size = n, replace = TRUE)
        y <- patch$y[i]
        z <- patch$z[i]
        theta.b[b] <- mean(y) / mean(z)
        }
    bias <- mean(theta.b) - theta.hat
    se <- sd(theta.b)
    print(list(est=theta.hat, bias = bias,
               se = se, cv = bias/se))


### Example 8.6 (Jackknife estimate of bias)

    data(patch, package = "bootstrap")
    n <- nrow(patch)
    y <- patch$y
    z <- patch$z
    theta.hat <- mean(y) / mean(z)
    print (theta.hat)

    #compute the jackknife replicates, leave-one-out estimates
    theta.jack <- numeric(n)
    for (i in 1:n)
        theta.jack[i] <- mean(y[-i]) / mean(z[-i])
    bias <- (n - 1) * (mean(theta.jack) - theta.hat)

    print(bias)  #jackknife estimate of bias
    
    
### Example 8.7 (Jackknife estimate of standard error)

    se <- sqrt((n-1) *
        mean((theta.jack - mean(theta.jack))^2))
    print(se)
    

### Example 8.8 (Failure of jackknife)

    set.seed(123) #for the specific example given
    #change the seed to see other examples

    n <- 10
    x <- sample(1:100, size = n)

    #jackknife estimate of se
    M <- numeric(n)
    for (i in 1:n) {        #leave one out
        y <- x[-i]
        M[i] <- median(y)
    }
    Mbar <- mean(M)
    print(sqrt((n-1)/n * sum((M - Mbar)^2)))

    #bootstrap estimate of se
    Mb <- replicate(1000, expr = {
            y <- sample(x, size = n, replace = TRUE)
            median(y) })
    print(sd(Mb)) 
    print(x)
    print(M)
    print(Mb)


### Example 8.9 (Bootstrap confidence intervals for patch ratio statistic)

    library(boot)       #for boot and boot.ci
    data(patch, package = "bootstrap")

    theta.boot <- function(dat, ind) {
        #function to compute the statistic
        y <- dat[ind, 1]
        z <- dat[ind, 2]
        mean(y) / mean(z)
    }

    y <- patch$y
    z <- patch$z
    dat <- cbind(y, z)
    boot.obj <- boot(dat, statistic = theta.boot, R = 2000)

    print(boot.obj)
    print(boot.ci(boot.obj,
                  type = c("basic", "norm", "perc")))


    #calculations for bootstrap confidence intervals
    alpha <- c(.025, .975)

    #normal
    print(boot.obj$t0 + qnorm(alpha) * sd(boot.obj$t))

    #basic
    print(2*boot.obj$t0 - 
        quantile(boot.obj$t, rev(alpha), type=1))

    #percentile
    print(quantile(boot.obj$t, alpha, type=6))


### Example 8.10 (Bootstrap confidence intervals for the correlation statistic)

    library(boot)
    data(law, package = "bootstrap")
    boot.obj <- boot(law, R = 2000,
             statistic = function(x, i){cor(x[i,1], x[i,2])})
    print(boot.ci(boot.obj, type=c("basic","norm","perc")))


### Example 8.11 (Bootstrap t confidence interval)

    boot.t.ci <-
    function(x, B = 500, R = 100, level = .95, statistic){
        #compute the bootstrap t CI
        x <- as.matrix(x);  n <- nrow(x)
        stat <- numeric(B); se <- numeric(B)

        boot.se <- function(x, R, f) {
            #local function to compute the bootstrap
            #estimate of standard error for statistic f(x)
            x <- as.matrix(x); m <- nrow(x)
            th <- replicate(R, expr = {
                i <- sample(1:m, size = m, replace = TRUE)
                f(x[i, ])
                })
            return(sd(th))
        }

        for (b in 1:B) {
            j <- sample(1:n, size = n, replace = TRUE)
            y <- x[j, ]
            stat[b] <- statistic(y)
            se[b] <- boot.se(y, R = R, f = statistic)
        }
        stat0 <- statistic(x)
        t.stats <- (stat - stat0) / se
        se0 <- sd(stat)
        alpha <- 1 - level
        Qt <- quantile(t.stats, c(alpha/2, 1-alpha/2), type = 1)
        names(Qt) <- rev(names(Qt))
        CI <- rev(stat0 - Qt * se0)
    }


### Example 8.12 (Bootstrap t confidence interval for patch ratio statistic)

    #boot package and patch data were loaded in Example 8.10
    #library(boot)       #for boot and boot.ci
    #data(patch, package = "bootstrap")
    
    dat <- cbind(patch$y, patch$z)
    stat <- function(dat) {
        mean(dat[, 1]) / mean(dat[, 2]) }
    ci <- boot.t.ci(dat, statistic = stat, B=2000, R=200)
    print(ci)


### Example 8.13 (BCa bootstrap confidence interval)

    boot.BCa <-
    function(x, th0, th, stat, conf = .95) {
        # bootstrap with BCa bootstrap confidence interval
        # th0 is the observed statistic
        # th is the vector of bootstrap replicates
        # stat is the function to compute the statistic

        x <- as.matrix(x)
        n <- nrow(x) #observations in rows
        N <- 1:n
        alpha <- (1 + c(-conf, conf))/2
        zalpha <- qnorm(alpha)

        # the bias correction factor
        z0 <- qnorm(sum(th < th0) / length(th))

        # the acceleration factor (jackknife est.)
        th.jack <- numeric(n)
        for (i in 1:n) {
            J <- N[1:(n-1)]
            th.jack[i] <- stat(x[-i, ], J)
        }
        L <- mean(th.jack) - th.jack
        a <- sum(L^3)/(6 * sum(L^2)^1.5)

        # BCa conf. limits
        adj.alpha <- pnorm(z0 + (z0+zalpha)/(1-a*(z0+zalpha)))
        limits <- quantile(th, adj.alpha, type=6)
        return(list("est"=th0, "BCa"=limits))
    }
    

### Example 8.14 (BCa bootstrap confidence interval)

    #boot package and patch data were loaded in Example 8.10
    #library(boot)       #for boot and boot.ci
    #data(patch, package = "bootstrap")

    n <- nrow(patch)
    B <- 2000
    y <- patch$y
    z <- patch$z
    x <- cbind(y, z)
    theta.b <- numeric(B)
    theta.hat <- mean(y) / mean(z)

    #bootstrap
    for (b in 1:B) {
        i <- sample(1:n, size = n, replace = TRUE)
        y <- patch$y[i]
        z <- patch$z[i]
        theta.b[b] <- mean(y) / mean(z)
        }
    #compute the BCa interval
    stat <- function(dat, index) {
        mean(dat[index, 1]) / mean(dat[index, 2])  }

    boot.BCa(x, th0 = theta.hat, th = theta.b, stat = stat)
    
    

### Example 8.15 (BCa bootstrap confidence interval using boot.ci)

    #using x from Example 8.15
    boot.obj <- boot(x, statistic = stat, R=2000)
    boot.ci(boot.obj, type=c("perc", "bca"))

### Example 8.16 (Model selection)

    #to prompt for next graph, uncomment line below
    #par(ask = TRUE)   
    
    library(DAAG); attach(ironslag)
    a <- seq(10, 40, .1)     #sequence for plotting fits

    L1 <- lm(magnetic ~ chemical)
    plot(chemical, magnetic, main="Linear", pch=16)
    yhat1 <- L1$coef[1] + L1$coef[2] * a
    lines(a, yhat1, lwd=2)

    L2 <- lm(magnetic ~ chemical + I(chemical^2))
    plot(chemical, magnetic, main="Quadratic", pch=16)
    yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2
    lines(a, yhat2, lwd=2)

    L3 <- lm(log(magnetic) ~ chemical)
    plot(chemical, magnetic, main="Exponential", pch=16)
    logyhat3 <- L3$coef[1] + L3$coef[2] * a
    yhat3 <- exp(logyhat3)
    lines(a, yhat3, lwd=2)

    L4 <- lm(log(magnetic) ~ log(chemical))
    plot(log(chemical), log(magnetic), main="Log-Log", pch=16)
    logyhat4 <- L4$coef[1] + L4$coef[2] * log(a)
    lines(log(a), logyhat4, lwd=2)

### Example 8.17 (Model selection: Cross validation)

    # Example 8.16, cont.
    n <- length(magnetic)   #in DAAG ironslag
    e1 <- e2 <- e3 <- e4 <- numeric(n)

    # for n-fold cross validation
    # fit models on leave-one-out samples
    for (k in 1:n) {
        y <- magnetic[-k]
        x <- chemical[-k]

        J1 <- lm(y ~ x)
        yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
        e1[k] <- magnetic[k] - yhat1

        J2 <- lm(y ~ x + I(x^2))
        yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
                J2$coef[3] * chemical[k]^2
        e2[k] <- magnetic[k] - yhat2

        J3 <- lm(log(y) ~ x)
        logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
        yhat3 <- exp(logyhat3)
        e3[k] <- magnetic[k] - yhat3

        J4 <- lm(log(y) ~ log(x))
        logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
        yhat4 <- exp(logyhat4)
        e4[k] <- magnetic[k] - yhat4
    }


    c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

    #selected model, fitted in Example 8.16
    L2

    par(mfrow = c(2, 2))    #layout for graphs
    plot(L2$fit, L2$res)    #residuals vs fitted values
    abline(0, 0)            #reference line
    qqnorm(L2$res)          #normal probability plot
    qqline(L2$res)          #reference line
    par(mfrow = c(1, 1))    #restore display
    
