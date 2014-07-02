################################################################
# Exercise script: Statistical programming for the life sciences
#
# Matthias Duschl & Daniel Lee, 2014
################################################################

## Practice exercise 1: Apply what you've learned to series of random numbers.
# 1. Generate a series of 100 normally distributed numbers centered at 0 with a
#    standard deviance of 200.
series <- rnorm(100, sd = 200)
# 2. What's their standard deviance? Is it really what you specified?
sd(series)
# 3. What's their variance?
var(series)
# 4. Sort them in ascending order.
sort(series)
# 5. Sort them in descending order.
sort(series, decreasing = TRUE)
# 6. Filter out negative values.
series[series > 0]
# 7. Do the exercises again for a series of exponentially distributed values,
#    this time filtering for values under 1.
series <- rexp(100)
sd(series)
var(series)
sort(series)
sort(series, decreasing = TRUE)
series[series > 1]

## Practice exercise 2: Do some simple statistics on real data
# 1. Record the occurences of influenza-like illnesses per 100,000 (ILI) for
#    all Scandinavian and central European countries.
# Scandinavian countries
denmark <- 12.9
finland <- NA
iceland <- NA
norway <- 19.5
sweden <- 1.8
# Central Europe
belgium <- 13.9
germany <- NA
luxembourg <- 0.3
netherlands <- 16.5
switzerland <- 2.7
# 2. Make a vector containing the ILI for each country in both regions
scandinavia <- c(denmark, finland, iceland, norway, sweden)
central.europe <- c(belgium, germany, luxembourg, netherlands, switzerland)
# 3. Compute the average ILI for each region
# Watch out for NAs!
avg.scandinavia <- mean(scandinavia)  # Wrong: this will produce NA
avg.scandinavia <- mean(scandinavia, na.rm = TRUE)
avg.central.europe <- mean(central.europe, na.rm = TRUE)
# 4. Compute the ratio of ILIs between the Netherlands and Belgium
netherlands / belgium
# 5. Compute the ratio of ILIs between the Netherlands and every other country
#    in our sample
# Combine both vectors
sample.group <- c(scandinavia, central.europe)
# "netherlands" vector is automatically recycled to match length of
# "sample.group"
netherlands / sample.group

## Practice exercise 3 (1/3): Do some simple statistics on data frame columns
# 1. Find the maximum of all numeric columns
apply(sample.table.numbers, 2, max, na.rm = TRUE)
# 2. Find the minimum of all numeric columns
apply(sample.table.numbers, 2, min, na.rm = TRUE)
# 3. Find the variance of all numeric columns
apply(sample.table.numbers, 2, var, na.rm = TRUE)
# 4. Find the standard deviance of all numeric columns
apply(sample.table.numbers, 2, sd, na.rm = TRUE)

## Practice exercise 3 (2/3): Add columns to a data frame
## In this exercise, we will duplicate some of the columns in the file
## "experiment_data.xlsx" in the sheet "Plier-Algorithm". The results will
## be added to our "sample.table".
# 1. Compute column G
sample.table$MW_Hang_30min <- (sample.table$Hang_30min_Set1 +
                               sample.table$Hang_30min_Set2) / 2
# 2. Compute column H
sample.table$MW_Biof_30min <- (sample.table$Biof_30min_Set1 +
                               sample.table$Biof_30min_Set2) / 2
# 3. Compute column I
sample.table$MW_Hang_minus_Biof_30min <- abs(sample.table$MW_Hang_30min -
                                             sample.table$MW_Biof_30min) 
# 4. Compute column J
sample.table$log2_hang_biof_30min <- log(sample.table$MW_Biof_30min /
                                         sample.table$MW_Hang_30min, base = 2)
# 5. Compute column M
sample.table$fc.hang.biof <- sample.table$MW_Biof_30min / sample.table$MW_Hang_30min
# 6. Compute column O
# Digression: cbind, rbind
# There are multiple ways of solving each problem in R. Here, we use "cbind" to
# attach two columns to each other as a small "mini-dataframe" and then use
# "apply" to compute their maximum. "cbind" binds vectors of equal length
# together as columns, whereas its relative "rbind" does the same for rows.
g.h <- sample.table$MW_Hang_30min / sample.table$MW_Biof_30min
h.g <- sample.table$MW_Biof_30min / sample.table$MW_Hang_30min
sample.table$fc.biof.hang <- apply(cbind(g.h, h.g), 1, max)
# 7. Save the completed data frame
write.csv(sample.table, "sample_table.csv")

## Practice exercise 3 (3/3): Sorting data frames
## Sorting data frames is the same as sorting any vector. It is important,
## however, that the sort vector is used as the column coordinate, not the row
## coordinate. See above for a refresher on how to do this using the function
## "order".
# 1. Sort the data frame by "Strain"
auswertung[order(auswertung$Strain),]
# 2. Sort by "Gene"
auswertung[order(auswertung$Gene),]
# 3. Sort by "CT"
auswertung[order(auswertung$CT),]
# 4. Sort by "reference.target"
auswertung[order(auswertung$reference.target),]
# 5. Sort by "delta.CT"
auswertung[order(auswertung$delta.CT),]


## Practice exercise 4 (1/2): Correlations
## Compute and compare the correlation coefficients between Hang and Biof 
## after 30min of exposure with the one after 1 hour. To get more precise
## estimates of the correlation, use the average between Set1 and Set2. 
## Can you observe any difference in the correlation estimate?
# 1. Calculate the averages between Set1 and Set2
Hang_30min <- apply(sample.table[,c(3,4)],1,mean)
Hang_1h <- apply(sample.table[,c(5,6)],1,mean)
Biof_30min <- apply(sample.table[,c(7,8)],1,mean)
Biof_1h <- apply(sample.table[,c(9,10)],1,mean)
# 2. Correlation analysis after 30 min
cor.test(Hang_30min, Biof_30min, method="spearman")
# 3. Correlation analysis after 1 hour
cor.test(Hang_1h, Biof_1h, method="spearman")



## Practice exercise 4 (1/2): heatmaps
## Heatmap with dendogram from cluster algorithm are implemented in heatmap()
?heatmap
# The examples in the documentation uses mtcars dataset, which comes with R
# This dataframe must the transformed into a matrix
x  <- as.matrix(mtcars)
# T
# Make a heatmap (scaling the columns for better visual comparisons)
heatmap(x, col = cm.colors(256), scale = "column",
        xlab = "specification variables", ylab =  "Car Models")
# See the help file for further options, like ColSideColors or differetnt 
# cluster analysis specifications



############
# Challenges
############

## These are special challenges for the adventurous among you. They all work
## with data from the file "experimental_data.csv". You might still have the
## raw data loaded in the "sample.table" data frame. They all use techniques
## that aren't explicitly covered in the workshop slides. Good luck!

# 1. Compute column N
# We use boolean values rather than a string
# Note that comparisons are lowest in the order of precedence
sample.table$hang.over.biof <-
    sample.table$MW_Biof_30min / sample.table$MW_Hang_30min > 1

# 2. Compute column L
average.bc.bd <- apply(sample.table[, 3:4], 1, mean)
average.bi.bj <- apply(sample.table[, 7:8], 1, mean)
abs.bc.bd <- abs(sample.table$Hang_30min_Set1 - sample.table$Hang_30min_Set2)
abs.bi.bj <- abs(sample.table$Biof_30min_Set1 - sample.table$Biof_30min_Set2)
# The functions "cbind" and "rbind" bind columns or rows together, respectively,
# as a matrix
max.abs <- apply(cbind(abs.bc.bd, abs.bi.bj), 1, max)
sample.table$signifikanz <- abs(average.bc.bd - average.bi.bj) >= max.abs * 1.5

# 3. Compute column K
# The t-test is a bit complex. We will need to work with the object produced by
# the "t.test" function.
?t.test
# These are the right columns...
i <- t.test(sample.table[, c(3:4)], sample.table[, c(7:8)])
# The result is a t.test object. Examine its structure.
str(i)
# The value you want is the p-value. You can access it using list notation, as
# if it were a column in a data frame.
i$p.value
# We see that this works, but we can the values for the entire series, not just
# for a subsection. The "apply" function is well-suited to this job.
# First, extract the columns you're interested in.
t.test.columns <- sample.table[, c(3:4, 7:8)]
# "apply" requires a function that is called on the data it creates for each
# iteration. The data created is passed to the function as the first and only
# argument. You can either pass function defined by R, or your own function
# that you defined yourself. If you want to use your own function, you can
# even define it for one-time use right inside the function.
sample.table$t.test <-
  apply(t.test.columns, 1,
        # This is our own custom function. The name of the argument it receives
        # is irrelevant. It computes the t-test for the first 2 samples against
        # the second 2 and returns the p-value.
        function(x) {t.test(x[1:2], x[3:4], paired = TRUE)$p.value})
