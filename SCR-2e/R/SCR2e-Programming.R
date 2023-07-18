#######################################################
###       Statistical Computing with R  2e          ###
###       Maria L. Rizzo                            ###
###       Chapman & Hall/CRC The R Series           ###
###       ISBN 9781466553323 - CAT# K15269          ###
###       January 2019                              ###
###                                                 ###
###       R code for Chapter 15                     ###
###       Programming Topics                        ###
#######################################################

# packages to install:
# energy, ggplot2, profvis, pryr, microbenchmark, mvtnorm, rbenchmark
# for the last section:  dplyr, Lahman, Rcpp 

### Example 15.1 (Benchmarking methods to generate a sequence)

    s1 <- 1:10
    s2 <- seq(1, 10, 1)
    s3 <- seq.int(1, 10, 1)
    df <- data.frame(s1=s1, s2=s2, s3=s3)
    str(df)

    library(microbenchmark)
    library(ggplot2)
    
    n <- 1000
    mb <- microbenchmark(
      seq(1, n, 1),
      seq.int(1, n, 1),
      1:n
    )
    
    print(mb)
    autoplot(mb)  # display a violin plot

### Example 15.2 (Benchmarking methods to initialize a vector)

    n <- 100
    mb2 <- microbenchmark(
      numeric = numeric(n) + 1,
      rep = rep(1, n),
      seq = seq(from=1, to=1, length=n),
      ones = matrix(1, nrow=n, ncol=1),
      as.ones = as.matrix(rep(1, n))
    )
    
    print(mb2)

### Example 15.3 (Timings of two multivariate normal generators)

  library(rbenchmark)
  library(MASS)
  library(mvtnorm)
  n <- 100          #sample size
  d <- 30           #dimension
  N <- 2000         #iterations
  mu <- numeric(d)

  benchmark(
      columns = c("test", "replications", "elapsed", "relative"),
      replications = 2000,
      cov = {S <- cov(matrix(rnorm(n*d), n, d))},
      mvrnorm = mvrnorm(n, mu, S),
      rmvnorm = rmvnorm(n, mu, S)
  )

### Example 15.4 (Profiling with Rprof)

  x <- rnorm(1000)
  y <- rnorm(1000)

  Rprof("pr.out", line.profiling = TRUE)
  energy::dcor(x, y)
  Rprof(NULL)
  summaryRprof("pr.out")

### Example 15.5 (profvis interactive visualization)

  library(profvis)
  profvis(energy::dcor(x, y))

### Example 15.6 (Object size)

  x <- matrix(rnorm(5000), 1000, 5) #1000 obs in R^5
  DF <- as.data.frame(x)
  object.size(x)
  object.size(DF)
  pryr::object_size(x)
  pryr::object_size(DF)
  pryr::compare_size(x)

  listTwo <- list(x, x)
  pryr::compare_size(listTwo)
 
### Example 15.7 (Comparing objects and attributes)
  
  str(x)
  str(DF)
  all.equal(x, DF)
  names(attributes(x))
  names(attributes(DF))
  all.equal(x, DF, check.attributes = FALSE)

### Example 15.8 (Comparing objects for equality)

  try(ifelse(all.equal(x, DF), "T", "F"))  # error
  ifelse(isTRUE(all.equal(x, DF)), "T", "F") # correct

  x <- 1 - 10e-4
  y <- x + 2
  x == (y - 2)   # equal mathematically but
  isTRUE(all.equal(x, y - 2))  #gives expected result

  ## does not necessarily evaluate to TRUE or FALSE
  try(ifelse(all.equal(x, y), "T", "F"))
  ## returns TRUE or FALSE
  ifelse(isTRUE(all.equal(x, y)), "T", "F")

### Example 15.9 (Display R function code)

  nclass.scott

### Example 15.10 (RSiteSearch)

  if (interactive())
    RSiteSearch("ggcorr")

### Example 15.11 (UseMethod)

  body(density)
  args(density.default)
  body(density.default)

### Example 15.12 (Show methods)
  
  methods(t.test)
  getAnywhere(t.test.formula)
  body(stats:::t.test.formula)

### Example 15.13 (Object not found or not an exported object)

  try(perc.ci)
  try(boot::perc.ci)
  
  getAnywhere(perc.ci)
  args(boot:::perc.ci)
  body(boot:::perc.ci)
  boot:::perc.ci
  getFromNamespace("perc.ci", "boot")

### Example 15.14 (getS3method)

  library(microbenchmark)
  library(ggplot2)
  getAnywhere(autoplot)
  getS3method("autoplot", class = "microbenchmark")
  getAnywhere(autoplot.microbenchmark)

### Example 15.15 (.Primitive or .Internal)

  if (interactive())
    pryr::show_c_source(.Primitive(cumsum(x)))

### Example 15.16 (.Call, .External, .C or .Fortran)

  #dist  is implemented by a .Call to C_Cdist

  body(dist)

### Example 15.17 (A first Rcpp experiment)

# In RStudio use the menu
# File > New File > C++ file to display this code

### Example 15.18 (cppFunction)

    library(Rcpp)
    
    set.seed(1)
    x <- matrix(rnorm(20), nrow = 5, ncol = 4)
    
    cppFunction('double vecnorm(NumericVector v) {
      // compute the Euclidean norm of vector v
      int d = v.size();
      double s = 0.0;
      for (int i = 0; i < d; i++)
        s += v(i) * v(i);
      return sqrt(s);
    }')
    
    print(vecnorm(x[, 1]))
    print(apply(x, MARGIN = 2, FUN = "vecnorm"))


### Example 15.19 (sourceCpp)

# Create the C++ source file "printme.cpp" by editing the template
# in RStudio menu: File > New File > C++

  library(Rcpp)
  sourceCpp("..//examples-cpp//printme.cpp")
  x <- sample(1:5)
  print_me(x)

### Example 15.20 (Lahman baseball data)

  library(Lahman)
  str(Batting)
  # method 1
  S <- subset(Batting, Batting$yearID == 1999,
             select = c("playerID", "AB", "H"))

  # method 2 (dplyr)
  library(dplyr)
  Batting %>% filter(yearID == 1999) -> b

  # method 1
  AB <- as.vector(by(S$AB, S$playerID, FUN = sum))
  H <- as.vector(by(S$H, S$playerID, FUN = sum))
  S <- data.frame(playerID = unique(S$playerID), 
                  AB = AB, H = H, AVG = round(H / AB, 3),
                  stringsAsFactors = FALSE)
  S400 <- S[S$AB >= 400, ]

  # method 2 (dplyr)
  b %>% group_by(playerID) %>%
    summarize(AB = sum(AB), H = sum(H)) -> S
  S %>% mutate(AVG = round(H / AB, 3)) -> S
  S %>% filter(AB >= 400) -> S400

  # method 1
  o <- order(S400$AVG, decreasing = TRUE)
  S400 <- S400[o, ]
  top <- S400[1:10, ]

  # method 2 (dplyr)
  S400 %>% arrange(desc(AVG)) -> S400
  slice(S400, 1:10) -> top
  top

  People %>% select(playerID, nameFirst, nameLast) -> m
  top %>% inner_join(m) %>%
    select(nameFirst, nameLast, AVG)

### Example 15.21 (Comparison with microbenchmark)

# The source for this example is in R Markdown format
# Open "Lahman.Rmd" in RStudio and knit to see report
