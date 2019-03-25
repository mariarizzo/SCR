#######################################################
###       Statistical Computing with R  2e          ###
###       Maria L. Rizzo                            ###
###       Chapman & Hall/CRC The R Series           ###
###       ISBN 9781466553323 - CAT# K15269          ###
###       January 2019                              ###
###                                                 ###
###       R code for Chapter 5                      ###
###       Visualization of Multivariate Data        ###
#######################################################


# packages to install:
# bootstrap, corrplot, DAAG, FactoMineR, ggplot2, hexbin


### Example 5.1 (Scatterplot matrix)

    data(iris)
    #virginica data in first 4 columns of the last 50 obs.
    
    # not shown in text
    pairs(iris[101:150, 1:4])

    panel.d <- function(x, ...) {
        usr <- par("usr")
        on.exit(par(usr))
        par(usr = c(usr[1:2], 0, .5))
        lines(density(x))
     }
    
    # Fig. 5.1
    x <- scale(iris[101:150, 1:4])
    r <- range(x)
    pairs(x, diag.panel = panel.d, xlim = r, ylim = r)

    library(lattice)
    splom(iris[101:150, 1:4])    #plot 1

    #for all 3 at once, in color, plot 2
    splom(iris[,1:4], groups = iris$Species)

    # Fig. 5.2
    #for all 3 at once, black and white, plot 3
    splom(~iris[1:4], groups = Species, data = iris,
        col = 1, pch = c(1, 2, 3),  cex = c(.5,.5,.5))

### Example 5.2 (Correlation plots: Decathlon data)

    library(FactoMineR)  #decathlon data
    library(corrplot)
    data("decathlon")
    str(decathlon)
    
    corrMat <- cor(decathlon[, 1:10])
    corrplot(corrMat, type="upper", tl.col="black", tl.srt=45)
    corrplot(corrMat, type = "upper", method = "square",
         addCoef.col = "black", diag=FALSE)

### Example 5.3 (Plot bivariate normal density)

    #the standard BVN density
    f <- function(x,y) {
        z <- (1/(2*pi)) * exp(-.5 * (x^2 + y^2))
        }

    y <- x <- seq(-3, 3, length= 50)
    z <- outer(x, y, f)   #compute density for all (x,y)

    persp(x, y, z)        #the default plot

    persp(x, y, z, theta = 45, phi = 30, expand = 0.6,
          ltheta = 120, shade = 0.75, ticktype = "detailed",
          xlab = "X", ylab = "Y", zlab = "f(x, y)")


### Example 5.4 (Add elements to perspective plot)

     #store viewing transformation in M
     persp(x, y, z, theta = 45, phi = 30,
           expand = .4, box = FALSE) -> M

     #add some points along a circle
     a <- seq(-pi, pi, pi/16)
     newpts <- cbind(cos(a), sin(a)) * 2
     newpts <- cbind(newpts, 0, 1)  #z=0, t=1
     N <- newpts %*% M
     points(N[,1]/N[,4], N[,2]/N[,4], col=2)

     #add lines
     x2 <- seq(-3, 3, .1)
     y2 <- -x2^2 / 3
     z2 <- dnorm(x2) * dnorm(y2)
     N <- cbind(x2, y2, z2, 1) %*% M
     lines(N[,1]/N[,4], N[,2]/N[,4], col=4)

     #add text
     x3 <- c(0, 3.1)
     y3 <- c(0, -3.1)
     z3 <- dnorm(x3) * dnorm(y3) * 1.1
     N <- cbind(x3, y3, z3, 1) %*% M
     text(N[1,1]/N[1,4], N[1,2]/N[1,4], "f(x,y)")
     text(N[2,1]/N[2,4], N[2,2]/N[2,4], bquote(y==-x^2/3))


### Example 5.5 (Surface plot using wireframe(lattice))

    library(lattice)
    x <- y <- seq(-3, 3, length= 50)

    xy <- expand.grid(x, y)
    z <- (1/(2*pi)) * exp(-.5 * (xy[,1]^2 + xy[,2]^2))
    wireframe(z ~ xy[,1] * xy[,2])


### Example 5.6 (3D scatterplot)

    library(lattice)
    attach(iris)
    #basic 3 color plot with arrows along axes
    print(cloud(Petal.Length ~ Sepal.Length * Sepal.Width,
          data=iris, groups=Species))

    print(cloud(Sepal.Length ~ Petal.Length * Petal.Width,
        data = iris, groups = Species, main = "1", pch=1:3,
        scales = list(draw = FALSE), zlab = "SL",
        screen = list(z = 30, x = -75, y = 0)),
        split = c(1, 1, 2, 2), more = TRUE)

    print(cloud(Sepal.Width ~ Petal.Length * Petal.Width,
        data = iris, groups = Species, main = "2", pch=1:3,
        scales = list(draw = FALSE), zlab = "SW",
        screen = list(z = 30, x = -75, y = 0)),
        split = c(2, 1, 2, 2), more = TRUE)

    print(cloud(Petal.Length ~ Sepal.Length * Sepal.Width,
        data = iris, groups = Species, main = "3", pch=1:3,
        scales = list(draw = FALSE), zlab = "PL",
        screen = list(z = 30, x = -55, y = 0)),
        split = c(1, 2, 2, 2), more = TRUE)

    print(cloud(Petal.Width ~ Sepal.Length * Sepal.Width,
        data = iris, groups = Species, main = "4", pch=1:3,
        scales = list(draw = FALSE), zlab = "PW",
        screen = list(z = 30, x = -55, y = 0)),
        split = c(2, 2, 2, 2))
    detach(iris)


### Example 5.7 (Contour plot)

     #contour plot with labels
     contour(volcano, asp = 1, labcex = 1)

     #another version from lattice package
     library(lattice)
     contourplot(volcano) #similar to above


### Example 5.8 (Filled contour plots)

     image(volcano, col = terrain.colors(100), axes = FALSE)
     contour(volcano, levels = seq(100,200,by = 10), add = TRUE)

     filled.contour(volcano, color = terrain.colors, asp = 1)
     levelplot(volcano, scales = list(draw = FALSE),
               xlab = "", ylab = "")
               

### Example 5.9 (2D histogram)

    library(hexbin)
    x <- matrix(rnorm(4000), 2000, 2)
    plot(hexbin(x[,1], x[,2]))

    # ggplot version 
    library(ggplot2)
    x <- data.frame(x)
    ggplot(x, aes(x[,1], x[,2])) + geom_hex()

### Example 5.10 (Andrews curves)

    library(DAAG)
    attach(leafshape17)

    f <- function(a, v) {
        #Andrews curve f(a) for a data vector v in R^3
        v[1]/sqrt(2) + v[2]*sin(a) + v[3]*cos(a)
    }


    #scale data to range [-1, 1]
    x <- cbind(bladelen, petiole, bladewid)
    n <- nrow(x)
    mins <- apply(x, 2, min)  #column minimums
    maxs <- apply(x, 2, max)  #column maximums
    r <- maxs - mins          #column ranges
    y <- sweep(x, 2, mins)    #subtract column mins
    y <- sweep(y, 2, r, "/")  #divide by range
    x <- 2 * y - 1            #now has range [-1, 1]

    #set up plot window, but plot nothing yet
    plot(0, 0, xlim = c(-pi, pi), ylim = c(-3,3),
        xlab = "t", ylab = "Andrews Curves",
        main = "", type = "n")

    #now add the Andrews curves for each observation
    #line type corresponds to leaf architecture
    #0=orthotropic, 1=plagiotropic
    a <- seq(-pi, pi, len=101)
    dim(a) <- length(a)
    for (i in 1:n) {
        g <- arch[i] + 1
        y <- apply(a, MARGIN = 1, FUN = f, v = x[i,])
        lines(a, y, lty = g)
    }
    legend(3, c("Orthotropic", "Plagiotropic"), lty = 1:2)
    detach(leafshape17)
    

### Example 5.11 (Parallel coordinates)

    library(MASS)
    library(lattice)
    #trellis.device(color = FALSE) #black and white display
    x <- crabs[seq(5, 200, 5), ]  #get every fifth obs.
    parallelplot(~x[4:8] | sp*sex, x)
    
    #trellis.device(color = FALSE)    #black and white display
    x <- crabs[seq(5, 200, 5), ]     #get every fifth obs.
    a <- x$CW * x$CL                 #area of carapace
    x[4:8] <- x[4:8] / sqrt(a)       #adjust for size
    parallelplot(~x[4:8] | sp*sex, x)

### Example 5.12 (Segment plot)

    #segment plot
    x <- MASS::crabs[seq(5, 200, 5), ]         #get every fifth obs.
    x <- subset(x, sex == "M")           #keep just the males
    a <- x$CW * x$CL                     #area of carapace
    x[4:8] <- x[4:8] / sqrt(a)           #adjust for size
    
    #use default color palette or other colors
    #palette(gray(seq(.4, .95, len = 5))) #use gray scale
    palette(rainbow(6))                 #or use color
    stars(x[4:8], draw.segments = TRUE,
          labels =levels(x$sp), nrow = 4,
          ylim = c(-2,10), key.loc = c(3,-1))
    
    #after viewing, restore the default colors
    palette("default")
    
    
### Example 5.13 (PCA for open and closed book exams)
    
    library(bootstrap)
    str(scor)
    pairs(scor)   
    cor(scor)
    
    n <- nrow(scor)
    x <- scale(scor)  #center and scale
    s <- cov(x)
    e <- eigen(s)
    lam <- e$values   #vector of eigenvalues
    P <- e$vectors    #matrix of eigenvectors
    
    plot(lam, type = "b", xlab = "eigenvalues", main = "")
    barplot(lam, xlab = "eigenvalues")
    
    tab <- rbind(lam / sum(lam), cumsum(lam) / sum(lam))
    tab
    
    z <- x %*% P
    dim(z)
    head(z)
    
    pc <- prcomp(scor, center = TRUE, scale = TRUE)
    summary(pc)
    
    df <- scor[1:5, ]
    predict(pc, newdata = df)  #same as z above
    
    head(x %*% pc$rotation, 5)
    head(pc$rotation)
    head(P)
    

### Example 5.14 (PC Biplot)
    
    ## plot scor data in the (PC1, PC2) coordinate system
    biplot(pc, pc.biplot = TRUE) 
    round(cor(x, z), 3)

### Example 5.15 (PCA: Decathlon data)
    
    library(FactoMineR)
    data(decathlon)
    pc <- princomp(decathlon[,1:10], cor = TRUE, scores = TRUE)
    plot(pc)     # screeplot
    biplot(pc)
    summary(pc)    
    
    