#######################################################
###       Statistical Computing with R  2e          ###
###       Maria L. Rizzo                            ###
###       Chapman & Hall/CRC The R Series           ###
###       ISBN 9781466553323 - CAT# K15269          ###
###       January 2019                              ###
###                                                 ###
###       R code for Chapter 9                      ###
###       Resampling Applications                   ###
#######################################################


# packages to install: bootstrap, DAAG, ggplot2


### Example 9.1 (Jackknife-after-bootstrap)

library(boot)
library(bootstrap)
set.seed(1111)

theta.boot <- function(patch, i) {
  # function to compute the patch ratio statistic
  y <- patch[i, "y"]
  z <- patch[i, "z"]
  mean(y) / mean(z)
}

boot.out <- boot(bootstrap::patch,
                 statistic = theta.boot, R=2000)
A <- boot.array(boot.out)
head(A, 3)
mean(A[, 1] == 0)

# jackknife-after-bootstrap to est. se(se)
A <- boot.array(boot.out)
theta.b <- boot.out$t
n <- NROW(patch)
jack.se <- numeric(n)

for (i in 1:n) {
  #in i-th replicate omit all samples with x[i]
  keep <- which(A[, i] == 0)
  jack.se[i] <- sd(theta.b[keep])
}

print(boot.out)  #for se_boot
se.bar <- mean(jack.se)
se.se <- sqrt((n-1) * mean((jack.se - se.bar)^2))
print(paste("Jackknife-after-bootstrap est. se(se)=", se.se))


### Example 9.2 (Jackknife-after-bootstrap)

# initialize
data(patch, package = "bootstrap")
y <- patch$y
z <- patch$z
dat <- cbind(y, z)
n <- NROW(dat)
B <- 2000

# jackknife-after-bootstrap step 1: run the bootstrap
theta_boot <- function(dat, ind) {
  # function to compute the statistic
  y <- dat[ind, 1]
  z <- dat[ind, 2]
  mean(y) / mean(z)
}

boot.obj <- boot(dat, statistic = theta_boot, R=2000)
theta.hat <- boot.obj$t0
theta.b <- boot.obj$t
se.boot <- sd(theta.b)

# jackknife-after-bootstrap to est. se(se)
sample.freq <- boot.array(boot.obj)
se.se.reps <- numeric(n)
N <- 1:n

for (i in N) {
  # jackknife-after-bootstrap
  # omit all bootstrap samples that contain obs i
  keep <- which(sample.freq[ ,i] == 0)
  se.se.reps[i] <- sd(theta.b[keep])
}

print(boot.obj)
se.bar <- mean(se.se.reps)
se.se <- sqrt((n-1) * mean((se.se.reps - se.bar)^2))
se.se

### Example 9.3 (ironslag linear model)

library(ggplot2)
library(DAAG)
L1 <- lm(magnetic ~ chemical, data=ironslag)
cf3 <- round(L1$coeff, 3)
cap <- paste("Fit: magnetic =", cf3[1], "+", cf3[2], "chemical")

ggplot(data=ironslag, aes(chemical, magnetic)) +
geom_point() + geom_smooth(method="lm") +
ggtitle(cap)

plot(L1, which=1:2, ask=FALSE)  #residual plots

### Example 9.4 (mammals data)

library(MASS)
cor(log(mammals$body), log(mammals$brain))
summary(mammals)

y <- log(mammals$brain)
x <- log(mammals$body)
L <- lm(y ~ x)
L

cap <- paste("Fit: log(brain) =", round(L$coeff[1],3),
             "+", round(L$coeff[2],3), "log(body)")
ggplot(data=mammals, aes(x, y)) +
geom_point() + geom_smooth(method="lm") +
labs(x = "log(body)", y = "log(brain)", title = cap)

summary(L)$r.squared

### Example 9.5 (ironslag data, resampling cases)

x <- ironslag$chemical
y <- ironslag$magnetic
m <- 2000
n <- NROW(x)
L1 <- lm(y ~ x)  #estimate the model
b0 <- L1$coeff[1]; b1 <- L1$coeff[2]

## run bootstrap of cases
out <- replicate(m, expr={
i <- sample(1:n, replace=TRUE, size=n)
xstar <- x[i]
ystar <- y[i]
Lb <- lm(ystar ~ xstar)
s <- summary(Lb)$sigma
c(Lb$coeff[1], slope=Lb$coeff[2], s=s)
})

bootCases <- t(out)
meanCases <- colMeans(bootCases)
sdCases <- apply(bootCases, 2, "sd")
meanCases
sdCases

biasInt <- mean(bootCases[,1] - b0)  #bias for intercept
biasSlope <- mean(bootCases[,2] - b1)  #bias for slope

rbind(estimate=c(b0, b1), bias=c(biasInt, biasSlope),
      se=sdCases[1:2], cv=c(biasInt, cv=biasSlope)/sdCases[1:2])

### Example 9.6 (Resampling cases using the boot function)

# set.seed(1104)
library(boot)
m <- 2000
stats <- function(dat, i) {
x <- dat$chemical[i]
y <- dat$magnetic[i]
Lb <- lm(y ~ x)
s <- summary(Lb)$sigma
c(Lb$coeff[1], slope=Lb$coeff[2], s=s)
}

boot.out <- boot(ironslag, statistic=stats, R=2000)
boot.out
boot.out$t0

sd(boot.out$t[,2])
boottbl <- broom::tidy(boot.out)
boottbl$std.error[2]

MASS::truehist(boot.out$t[ ,2], main="", xlab="slopes")
abline(v = boot.out$t0[2], lwd=2)

boot.ci(boot.out, index=2, type=c("norm","perc","basic","bca"))

### Example 9.7 (Resampling errors: mammals data)

# Insert and run Example 9.4 code here

m.resid <- rstandard(L, sd = 1)
r <- m.resid - mean(m.resid)
m <- 1000; n <- NROW(x)
estsErr <- replicate(m, expr={
  estar <- sample(r, replace=TRUE, size=n)
  ystar <- L$fitted.values + estar
  Lb <- lm(ystar ~ x)
  s <- summary(Lb)$sigma
  c(b0=Lb$coeff[1], b1=Lb$coeff[2], s=s)
})
ests <- t(estsErr)
summary(ests)

### Example 9.8 (Resampling errors, continued)
sd(ests[,2])
s <- summary(L)$sigma
SSx <- (n - 1) * var(x)
se.beta1 <- sqrt(s^2 / SSx)
se.beta1
s * sqrt(1/n + mean(x)^2 / SSx)
sd(ests[,1])

betas <- summary(L)$coeff
betas
betas[, "Std. Error"]

broom::tidy(summary(L))
broom::tidy(summary(L))$std.error

### Example 9.9 (Model based resampling with the boot function)
regstats <- function(dat, i) {
  #dat is a data frame (r, x, yhat)
  #r are the modified centered residuals, yhat are the fits
  ystar <- dat$yhat[i] + dat$r[i]
  xstar <- dat$x[i]
  Lnew <- lm(ystar ~ xstar)
  Lnew$coefficients
}

y <- log(mammals$brain)
x <- log(mammals$body)
L <- lm(y ~ x)
r <- rstandard(L, sd=1)
r <- r - mean(r)
df <- data.frame(r=r, x=x, yhat=L$fitted)
head(df)
boot.obj <- boot(data=df, statistic=regstats, R=2000)
broom::tidy(boot.obj)


### Example 9.10 (Empirical influence values for the patch ratio statistic)

library(boot)
library(bootstrap)
theta_boot <- function(dat, ind) {
  # function to compute the patch ratio statistic
  mean(dat[ind, ]$y) / mean(dat[ind, ]$z)
}
boot.out <- boot(patch, theta_boot, R = 2000)
infl <- empinf(boot.out, type = "jack")
theta.hat <- boot.out$t0
jack <- theta.hat - infl / (nrow(patch) - 1)
rbind(infl, jack)

### Example 9.11 (Jackknife-after-bootstrap plot)

jack.after.boot(boot.out, useJ=TRUE, stinf=FALSE)

n <- NROW(patch) 
J <- numeric(n)
b.freq <- boot.array(boot.out)
theta.b <- boot.out$t

for (i in 1:n) {
  keep <- which(b.freq[ ,i] == 0)
  J[i] <- mean(theta.b[keep])
}

# the jackknife influence values
(n - 1) * (mean(J) - J)

jack.after.boot(boot.out, useJ=TRUE, stinf=TRUE)

