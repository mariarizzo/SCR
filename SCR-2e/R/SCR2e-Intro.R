#######################################################
###       Statistical Computing with R  2e          ###
###       Maria L. Rizzo                            ###
###       Chapman & Hall/CRC The R Series           ###
###       ISBN 9781466553323 - CAT# K15269          ###
###       March 6, 2019                             ###
###                                                 ###
###       R code for Chapter 1                      ###
###       Introduction                              ###
#######################################################


# packages to install: ggplot2


### Example 1.1 

    sumdice <- function(n) {
        k <- sample(1:6, size=n, replace=TRUE)
        return(sum(k))
    }

    sumdice(2)

    #to store the result rather than print it
    a <- sumdice(100)

    #we expect the mean for 100 dice to be close to 3.5
    a / 100

    sumdice <- function(n)
        sum(sample(1:6, size=n, replace=TRUE))

    sumdice <- function(n, sides = 6) {
        if (sides < 1) return (0)
        k <- sample(1:sides, size=n, replace=TRUE)
        return(sum(k))
    }

    sumdice(5)      #default 6 sides
    sumdice(n=5, sides=4)  #4 sides

### Example 1.2 (iris data)
    
  names(iris)
  table(iris$Species)
  w <- iris[[2]]   #Sepal.Width
  mean(w)

  attach(iris)
  summary(Petal.Length[51:100]) #versicolor petal length

  with(iris, summary(Petal.Length[51:100]))

  out <- with(iris, summary(Petal.Length[51:100]))

  by(iris[,1:4], Species, colMeans)
  detach(iris)

### Example 1.3 (Arrays)

    x <- 1:24                      # vector
    dim(x) <- length(x)            # 1 dimensional array
    matrix(1:24, nrow=4, ncol=6)   # 4 by 6 matrix
    x <- array(1:24, c(3, 4, 2))   # 3 by 4 by 2 array

### Example 1.4 (Matrices)

    A <- matrix(0, nrow=2, ncol=2)
    A <- matrix(c(0, 0, 0, 0), nrow=2, ncol=2)
    A <- matrix(0, 2, 2)

    A <- matrix(1:8, nrow=2, ncol=4)

### Example 1.5 (Iris data cont.)
    
   x <- as.matrix(iris[,1:4]) #all rows of columns 1 to 4

   mean(x[,2])         #mean of sepal width, all species
   mean(x[51:100,3])   #mean of petal length, versicolor

   y <- array(x, dim=c(50, 3, 4))
   mean(y[,,2])  #mean of sepal width, all species
   mean(y[,2,3]) #mean of petal length, versicolor

   y <- array(c(x[1:50,], x[51:100,], x[101:150,]),
              dim=c(50, 4, 3))
   mean(y[,2,])  #mean of sepal width, all species
   mean(y[,3,2]) #mean of petal length, versicolor

### Example 1.6 (Run length encoding)
   
 
    n <- 1000
    x <- rbinom(n, size = 1, prob = .5)
    table(x)
    head(x, 30)

    r <- rle(x)
    str(r)

    head(r$lengths)
    head(r[[1]])

    max(r$lengths)
    log2(length(x))

### Example 1.7 (Named list)

    w <- wilcox.test(rnorm(10), rnorm(10, 2))
    w    #print the summary

    w$statistic       #stored in object w
    w$p.value
    unlist(w)
    unclass(w)

### Example 1.8 (A list of names)

    a <- matrix(runif(8), 4, 2)   #a 4x2 matrix
    dimnames(a) <- list(NULL, c("x", "y"))

    # if we want row names
    dimnames(a) <- list(letters[1:4], c("x", "y"))
    a

    # another way to assign row names
    row.names(a) <- list("NE", "NW", "SW", "SE")
    a

### Example 1.9 (Parallel boxplots)

    boxplot(iris$Sepal.Length ~ iris$Species)
 
    boxplot(iris$Sepal.Length ~ iris$Species,
            ylab = "Sepal Length", boxwex = .4)
### Example 1.10 (Plotting characters and colors)

    plot(0:25, rep(1, 26), pch = 0:25)
    text(0:25, 0.9, 0:25)


### Example 1.11 (Barplot for run lengths)

  barplot(table(r$lengths))  #R graphics version

  ## ggplot version
  library(ggplot2)
  df <- data.frame(lengths = factor(r$lengths))
  ggplot(df, aes(lengths)) + geom_bar()

### Example 1.12 (Scatterplots)  

  ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width)) +
        geom_point()

  ggplot(iris, aes(Sepal.Length, Sepal.Width,
         color = Species, shape = Species)) + geom_point(size = 2)

### Example 1.13 (ggplot: parallel boxplots and violin plots)

  ggplot(iris, aes(Species, Sepal.Length)) + geom_boxplot()
  ggplot(iris, aes(Species, Sepal.Length)) + geom_violin()

  ggplot(iris, aes(Species, Sepal.Length)) +
     geom_boxplot() + coord_flip()
  ggplot(iris, aes(Species, Sepal.Length)) +
     geom_violin() + coord_flip()

### Example 1.14 (MPG by engine displacement)

  ggplot(mpg, aes(displ, hwy)) +
         geom_point() +
         facet_wrap(~ class)

### Example 1.15 (Import data from a local text file)
  
  # copy "FOREARM.DAT" into working directory 
  forearm <- scan(file = "FOREARM.DAT") #a vector
  
  # use pathname if file is not in working directory
  # forearm <- scan(file = "./DATASETS/FOREARM.DAT") #a vector
  
  head(forearm)

### Example 1.16 (Importing data from a web page)

fileloc <- "https://archive.ics.uci.edu/ml/machine-learning-databases/auto-mpg/auto-mpg.data"

df <- read.table(file = fileloc, na.strings = "?", as.is = TRUE)
str(df)
names(df) <- c("mpg", "cyl", "displ", "hp", "wt", "accel",
    "year", "origin", "name")
summary(df)

### Example 1.17 (Importing/exporting .csv files)

    #create a data frame
    dates <- c("3/27/1995", "4/3/1995",
               "4/10/1995", "4/18/1995")
    prices <- c(11.1, 7.9, 1.9, 7.3)
    d <- data.frame(dates=dates, prices=prices)

    #create the .csv file
    filename <- "temp.csv"
    write.table(d, file = filename, sep = ",",
                row.names = FALSE)

    #read the .csv file
    read.table(file = filename, sep = ",", header = TRUE)
    read.csv(file = filename) #same thing

    