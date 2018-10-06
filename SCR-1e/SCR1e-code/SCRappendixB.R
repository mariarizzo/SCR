    #######################################################
    ###       Statistical Computing with R              ###
    ###       Maria L. Rizzo                            ###
    ###       Chapman & Hall / CRC                      ###
    ###       ISBN 9781584885450                        ###
    ###                                                 ###
    ###       R code for Appendix B Examples            ###
    #######################################################


### Example B.2 (Extracting rows from a matrix)

    A <- matrix(1:25, 5, 5)
    A[-(2:3), ]
    A
    
    A[-(2:3), 4]
    A

### Example B.3 (Subsetting data frames)

    # versicolor petal length
    y <- subset(iris, Species == "versicolor",
                select = Petal.Length)

    summary(y)

    # sepal width, all species
    y <- subset(iris, select = c(Sepal.Length, Sepal.Width))
    mean(y)


### Example B.4 (InsectSprays data)

    attach(InsectSprays)
    InsectSprays
    unstack(count, count ~ spray)


### Example B.5 (Merge by ID)

    data1 <- cbind(1:5, rbinom(5, 20, .5))
    data2 <- cbind(3:7, rbinom(5, 20, .5))
    data1
    data2

    # keep only the common ID's
    merge(data1, data2, by=c(1,1))

    #keep all observations
    merge(data1, data2, by=c(1,1), all=TRUE)


### Example B.6 (Reshape)

    #keep all observations
    a <- merge(data1, data2, by=c(1,1), all=TRUE)
    reshape(a, idvar="ID", varying=c(2,3),
        direction="long", v.names="Scores")


### Example B.7 (Recode)

    #store the previous result into b
    b <- reshape(a, idvar="ID", varying=c(2,3),
                  direction="long", v.names="Scores")
    i <- which(is.na(b$Scores))   #these are missing

    b$Scores[i] <- 0   #replace NA with 0
    b

    m <- merge(data1, data2, by=c(1,1), all=TRUE)
    i <- which(is.na(m), arr.ind=TRUE)   #these are missing
    i


### Example B.8 (Date formats)

    d <- "3/27/1995"
    thedate <- as.Date(d, "%m/%d/%Y")
    print(thedate)
    print(format(thedate, "%Y%m%d"))
    print(format(thedate, "%B %d, %Y"))
    print(format(thedate, "%y-%b-%d"))


### Example B.9 (Date-time class)

    pdate <- as.POSIXlt(thedate)
    print(pdate$year)
    print(pdate$mon)
    print(pdate$mday)


### Example B.10 (Importing/exporting .csv files)

    #create a data frame
    dates <- c("3/27/1995", "4/3/1995",
               "4/10/1995", "4/18/1995")
    prices <- c(11.1, 7.9, 1.9, 7.3)
    d <- data.frame(dates=dates, prices=prices)

    #create the .csv file
    filename <- "/Rfiles/temp.csv"
    write.table(d, file = filename, sep = ",",
                row.names = FALSE)

    #read the .csv file
    read.table(file = filename, sep = ",", header = TRUE)
    read.csv(file = filename) #same thing


### Example B.11 (One-way ANOVA)

    # One-way ANOVA example
    # Completely randomized design
    ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
    trt1 <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
    trt2 <- c(5.19,3.33,3.20,3.13,6.46,5.36,6.95,4.19,3.16,4.95)

    group <- factor(rep(1:3, each=10)) #factor
    weight <- c(ctl, trt1, trt2)       #response
    a <- lm(weight ~ group)

    anova(a)                           #brief summary

    summary(a)                         #more detailed summary


### Example B.12 (Extract p-values and statistics from ANOVA)

    A <- anova(a)
    names(A)
    A$"F value"
    A$"F value"[1]

    B <- summary(a)
    names(B)
    B$sigma
    B$r.squared
    B$df[2]

### Example B.13 (Stacked data entry)

    P <- c(29.6, 24.3, 28.5, 32)
    T <- c(27.3, 32.6, 30.8, 34.8)
    S <- c(5.8,  6.2, 11, 8.3)
    E <- c(21.6, 17.4, 18.3, 19)
    C <- c(29.2, 32.8, 25,  24.2)

    #glue the columns together in a data frame
    x <- data.frame(P, T, S, E, C)

    #now stack the data for ANOVA
    y <- stack(x)
    names(y) <- c("Binding", "Antibiotic")

    #check the default formula
    print(formula(y))  #default formula is right one

    lm(y)
    lm(formula = y)
    anova(lm(y))


### Example B.14 (Two-way ANOVA)

    data(leafshape, package = "DAAG")
    attach(leafshape)
    anova(lm(petiole ~ location * arch))


### Example B.15 (Multiple comparisons)

    qtukey(p = .95, nmeans = 5, df = 15)

    #alternately: Tukey Honest Significant Difference

    a <- aov(formula(y), data = y)

    TukeyHSD(a, conf.level=.95)


### Example B.16 (Regression)

    library(DAAG)
    attach(ironslag)

    # simple linear regression model
    lm(magnetic ~ chemical)

    # quadratic regression model
    lm(magnetic ~ chemical + I(chemical^2))

    # exponential regression model
    lm(log(magnetic) ~ chemical)

    # log-log model
    lm(log(magnetic) ~ log(chemical))

    # cubic polynomial model
    lm(magnetic ~ poly(chemical, degree = 3))

    detach(ironslag)
    detach(package:DAAG)
    
