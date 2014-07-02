###############################################################
# Lecture script: Statistical programming for the life sciences
#
# Matthias Duschl & Daniel Lee, 2014
###############################################################

########################################
# Section II: Fundamentals
########################################

Working on the command line

## Expressions
## Every expression you run produces a result. If you don't "capture" it by
## saving it to a variable, it is reported back on the command line. Some
## examples:
1 + 2  # Addition
8 %% 3  # Remainder division
1:10  # Generate a series
mean(c(1, 2, 3))  # Use a function
## There are different data types. The most important are:
0  # Numeric
"zero"  # Character
FALSE  # Logical
## Incompatible types can't be used together in ways R won't expect
"a string" + 2

## Vectors
## Vectors are collections of singular elements that have the same type.
## Elements are combined using the "c" function, which stands for "concatenate"
## or "combine".
# Here I save a numeric vector to the variable "numbers".
numbers <- c(5, 79, 234, 150, 0, -986)
# The result of a variable name is its value.
numbers
# You can access individual elements of a vector by index number. R starts
# counting a 1, so this will return element 1 of "numbers".
numbers[1]
# Making that negative negates the query, so that all elements *except* the
# named position are returned.
numbers[-1]
# So if I want all but the 5th number it's like this
numbers[-5]
# You can also specify a range of index positions.
numbers[2:4]
# For higher granularity, provide a vector of precisely the positions you want.
numbers[c(1:3, 5)]
# Vectors with the same type can be concatenated to produce a new vector.
numbers2 <- -2:-10
new.numbers <- c(numbers, numbers2)
# Of course, the same goes for a vector of characters or logical values
some.strings <- c("This", "is", "a", "bunch", "of", "strings")
some.strings[2:4]
some.strings[-(2:4)]

## Functions
## Functions are collections of instructions. They keep you from having to
## repeat the same statements again and again. They can be saved to named
## variables like any other value. They are called by using the call operator
## "()" containing any arguments the function requires. The arguments are used
## inside the function, so that its results depend on the data passed to it.
## For information on any function, use the "?" operator. You'll be shown a
## description of the function. There are many standard functions that you can
## use, and if you like you can even write your own.
?mean  # How do I use the function "mean"?
mean(numbers)  # Mean
sd(numbers)  # Standard deviation
var(numbers)  # Variance

## NAs
## "NA" stands for "not available" and is used in R to signify a missing value.
## The value "NA" can be coerced to any type and can't produce a valid value in
## expressions.
## NAs are often used to stand for values that can't be computed or are not
## present in input data.
# NAs can be recorded like this
not.available <- NA
# NAs propogate when used in expressions.
1 + NA
mean(c(numbers, not.available))
# Often you can instruct functions to ignore them by setting the argument
# na.rm = TRUE
mean(c(numbers, not.available), na.rm = TRUE)

## Filtering and sorting
# Sort "numbers"
sort(numbers)
sort(numbers, decreasing = TRUE)
# If you want to replace the old vector with its sorted version, overwrite it
numbers <- sort(numbers)
# Filter "numbers"
# Filtering uses expressions that produce logical vectors to tell R which
# values from a vector to keep and which to throw away
numbers[numbers < 100]
numbers[numbers > 0]
# Logical statements can be combined with "&" for "and" and "|" for "or".
numbers[numbers < 100 & numbers > 0]
numbers[numbers > 100 | numbers < 0]
# Find out the length of filtered vector with "length"
length(numbers[numbers %% 2 == 0])

## Random numbers
## R has lots of standard functions for random number generation. There are
## several distribution types you can use - see "distributions" for details.
?distributions
# Generate a series of 100 uniformly distributed random numbers in the range
# from -500 to 500.
random.numbers <- runif(min = -500, max = 500, n = 100)


## SEE EXERCISE xx

############################################################################
# Section 2: Working with real data
# This section uses data from EuroFlu: The WHO European Influenza Network.
# You can download it from http://www.euroflu.org, following these links:
# Display bulletin --> Tables and graphs
# The relevant data is the column "ILI per 100 000".
# Type the relevant data manually.
############################################################################

## Filtering and sorting real data
## Suppose we want to not only produce numbers, but also identify candidates
## from our sample that fulfill certain criteria. In order to do that, we need
## to explicitly associate the country data with the country names.
# First, make a vector of the relevant countries
scandinavian.countries <- c("denmark", "finland", "iceland", "norway",
                            "sweden")
central.european.countries <- c("belgium", "germany", "luxembourg",
                                "netherlands", "switzerland")
sample.names <- c(scandinavian.countries, central.european.countries)
# First you can filter "sample.group" to find countries with ILI > 10
sample.group[sample.group > 10]
# In order to understand what's happening here, it's useful to look at the
# statement inside the brackets - it produces a logical vector showing which
# elements fulfil the given criteria. In the previous statement, we used that
# logical vector to select the elements that fit our critieria. However, the
# logical vector doesn't know which vector it was originally associated with -
# it can be used to filter any other vector as well.
sample.group > 10
# Use this principle to filter "sample.names" based on criteria that can only
# be discovered by analyzing the numbers stored in "sample.group"
sample.names[sample.group > 10]  # Country names with ILI > 10
# We can also order the names by ILI. This works using the same principles as
# filtering: First, we produce a vector that stores the rank of all elements
# from "sample.group" on the position that they were found. Order puts NAs
# last, so you might want to filter beforehand.
order(sample.group)
# We can apply that ordering to extract values from another vector of the same
# length
sample.names[order(sample.group)]  # Ascending order
sample.names[order(sample.group, decreasing = TRUE)]  # Descending order

# Now that we're finished, you might have noticed that it's a bit of work to
# track multiple vectors of associated data. In the next exercise, you'll learn
# to associate them intrinsically by storing them as a data frame - R's
# table-like structure.

##########################################
# Section 3: Working with imported data
##########################################

## The working directory
## R always works in a given folder. If you read from or write to any files,
## you must either specify their absolute path or a relative path. Relative
## paths are relative to the current working directory. Often, you can save a
## lot of typing by navigating into a working directory and using relative
## paths from there.
getwd()  # Find out current working directory
work.dir <- "~/Dropbox/Bio-Workshop"  # Set this to directory with your data
setwd(work.dir)

## Importing data frames
## Your first data frame will be a sheet from "experiment_data.xlsx". Open
## it up and export the sheet named "Tabelle1" as a comma separated value (CSV)
## file.
# Then produce the same results that are in the tables
csv.source <- "experiment_data.csv"  # Set this to the name you chose
# Often you can use "read.csv" to read (CSVs). The function relies on the data
# being formated the way R expects it.
sample.table <- read.csv(csv.source)
# If reading with "read.csv" doesn't work the way you want it to, read the help
# and try other arguments. If it still doesn't work, you can try using the more
# general function that "read.csv" uses internally: "read.table".
sample.table <- read.table(csv.source)
# This function is *very* general so it probably won't work the way you want it
# to without giving it any arguments. If we examine what it does with the
# default settings, it turns out to be wrong in our case:
head(sample.table)
# A glance at the help can give us some hints...
?read.table
# One problem here was that it wasn't recognizing the ';' as a separator.
sample.table <- read.table(csv.source, sep = ";")
head(sample.table)
# Another problem was that it wasn't recognizing the first line as the header.
sample.table <- read.table(csv.source, sep = ";", header = TRUE)
head(sample.table)
# If your numbers exported in English, that last command will have worked. If
# your computer is set to German, you might need to use the "dec" argument,
# which stands for the symbol used to separate decimal places.

## Examining data frames
# You can view the table like any other variable
sample.table
# Or if you don't want to see all of the data, just look at its head or tail
head(sample.table)
tail(sample.table)
# You can also get a quick summary of the table with "summary"
summary(sample.table)
# And you can get very precise information about the data types in your data
# frame with "str" (for "structure").
str(sample.table)
# Some other cool functions:
colnames(sample.table)  # Vector of data frame's column names
nrow(sample.table)  # Number of rows in data frame

## Working with subsets of data
## Often, you don't want to work with all the data in a data frame. R has many
## efficient ways of filtering out individual columns, rows, values or
## arbitrary subsections.
## Access a single column like this:
# By variable name
sample.table$Hang_30min_Set2
# By name string
sample.table[["Hang_30min_Set2"]]
# By column number
sample.table[,4]
sample.table[,-4]
## In fact, you can access any value you want with coordinates like that. The
## notation is ${ROW}, ${COLUMN}. If you leave out a dimension it
## gives you all values along that coordinate.
sample.table[42,]  # Row 42
sample.table[37, 7]  # Row 37, column 7
## Columns are just vectors, so we can pass them to a function without having
## to preprocess them.
mean(sample.table$Hang_30min_Set1, na.rm = TRUE)

## Operating on multiple columns / rows
## The "apply" family of functions allow us to repeat actions for several
## objects, e.g. across all rows or columns of a data frame
?apply
# The help tells us that we tell "apply" to do the same thing either for all
# rows or for all columns of a table-like structure. The arguments are passed
# like this:
# apply(${data}, ${number to row or column}, ${function}, ${other arguments})
apply(sample.table, 2, summary)  # Apply "summary" to each column
summary(sample.table)  # For comparison
# Watch out: "apply" tries to coerce data to having the same type. This doesn't
# always work. For example, this will produce an error:
apply(sample.table, 2, mean, na.rm = TRUE)
# This is because the first column is of a type that doesn't automatically
# convert to a number, so all the other columns are also assumed to not be
# numbers. We can prevent this from happening by removing the non-number
# columns.
apply(sample.table[, -c(1, 2)], 2, mean, na.rm = TRUE)
# If you think that looks hard to read, you're right. To be more clear about
# what you're doing, and to not have to type as much, you can assign the
# columns you chose to a variable.
sample.table.numbers <- sample.table[, -c(1, 2)]
apply(sample.table.numbers, 2, mean, na.rm = TRUE)

## Adding data to data frames
# In order to add a column to a data frame, just write to it as it were already
# there, using any of the ways of accessing data frames listed above.
sample.table$id <- 1:nrow(sample.table)  # Assign an ID number to each row
head(sample.table)

## Alternative ways to repeat tasks: for loops
## In this section we will use a loop to repeat the same commands on several
## subsections of a data frame
# We will be working on a table exported from the file
# "Auswertung Real-time Gesamt.xlsx". If you open it you'll see that it's
# formated in a way that will difficult for R to understand. We've prepared a
# version of the table that R can read easily by repeating the cells that were
# spread over a range into the range they used to occupy. You can load that
# without modification from the file "Auswertung Real-time Gesamt.csv".
# Remember to use the "as.is" flag to keep R from converting characters to
# factors.
auswertung <- read.csv("Auswertung Real-time Gesamt.csv", sep=";", as.is = TRUE)
head(auswertung)
# Now we will compute averages for each strain. First, we make a data frame
# containing only the gene that interests us.
cdc28 <- auswertung[auswertung$Gene == "CDC28",]
# Next, we make a vector of the unique strains with that gene.
strains <- unique(cdc28$Strain)
# The basic syntax for a for-loop looks like this:
for (strain in strains) print(strain)
# If we need to do more complex things inside the for loop we can enclose an
# entire code block in curvy brackets. Here we make a vector called
# "reference.target" and then loop over each strain in our strain vector. For
# each strain, we select the entries from our data frame that correspond to
# that strain, compute mean CT, and append the deviance from the mean for each
# sample to our "reference.target" vector. When we're finished, we attach the
# "reference.target" vector as a column to our data frame.
reference.target <- c()
for (strain in strains) {
  selection <- cdc28[cdc28$Strain == strain,]
  strain.mean <- mean(selection$CT)
  reference.target <- c(reference.target, strain.mean -
                        auswertung[auswertung$Strain == strain, 5])
}
auswertung$reference.target <- reference.target
head(auswertung)
# To complete the exercise, compute the squared change in CT.
auswertung$delta.CT <- auswertung$reference.target ^ 2


######################################
# Section 4: Data analysis: Statistics and graphics
######################################

## In this section you will learn how to use your data to perform some 
## statistical tests and to create some graphics. The most important tests
## (chi-squared, t-test, correlation, regression, anova) and graphic types
## (scatterplot, histogram, barchart, boxplot) are covered. You will also 
## export your visualizations for further processing. 


## First, we start by analysing the data from "experiment_data.csv". If you need
## to read them into R again, here is the code:
csv.source <- "experiment_data.csv"
sample.table <- read.table(csv.source, sep = ";", header = TRUE)
# A good approach for data analysis is to study the distributions of the 
# variables. Therefore, we can draw a histogram of one our variable:
hist(sample.table$Biof_30min_Set1)
# The data looks fairly skewed. Changing the number of bars makes the plot 
# more precise
hist(sample.table$Biof_30min_Set1, breaks=30)
# ... and different color more appealing
hist(sample.table$Biof_30min_Set1, breaks=30, col="blue")
# There are many ways to specify the color
hist(sample.table$Biof_30min_Set1, breaks=30, col=3)
hist(sample.table$Biof_30min_Set1, breaks=30, col="#756bb1")
hist(sample.table$Biof_30min_Set1, breaks=30, col=rgb(117,107,177, maxColorValue = 255))
# We can also change the title and axes labels
hist(sample.table$Biof_30min_Set1, breaks=30, col="#756bb1", 
     xlab="Biof 30min Set1", main="Histogram")
# ... and add some low-level graphics. Low-level graphics, like lines or rugs, 
# are added by calling the correspondung function after plotting the high-level
# graphic. Also the low-level graphics can be customized using paramters.
abline(v=mean(sample.table$Biof_30min_Set1), col="red", lwd=2)
abline(v=median(sample.table$Biof_30min_Set1), col="blue", lwd=2)
rug(sample.table$Biof_30min_Set1)
# Finally, to save the graphic, use the export functions in RStudio. 
# As a command-line alternative, try
dev.copy2pdf(file="histogram.pdf", width=7, height=5)

## Mostly, we want to study differences and relationships between more
## than one variable, for instance, between "Hang" and "Biof". 
# Again, we start by plotting the data, using a scatterplot
plot(sample.table$Hang_30min_Set1, sample.table$Biof_30min_Set1)
# The same plot with some customizations
plot(sample.table$Hang_30min_Set1, sample.table$Biof_30min_Set1,
     pch=19, cex=0.5, col="darkgrey", xlab="Hang 30min", ylab="Biof 30min")
# Adding a linear regression line
abline(lm(sample.table$Biof_30min_Set1 ~ sample.table$Hang_30min_Set1), col="blue", lwd=2, lty=2)

## At this point, we often want to make some correlation analysis
# cor() gives us the correlation coefficient as a descriptive statistic
cor(sample.table$Hang_30min_Set1, sample.table$Biof_30min_Set1)
# There are also some correlation tests. The Pearson test is called by
cor.test(sample.table$Hang_30min_Set1, sample.table$Biof_30min_Set1)
# However, we know that our data are not normally distribued. Hence, we
# have to use the Spearman rank correlation:
cor.test(sample.table$Hang_30min_Set1, sample.table$Biof_30min_Set1, method="spearman")
# The results from this test can be stored in new variables
cor_results <- cor.test(sample.table$Hang_30min_Set1, sample.table$Biof_30min_Set1, method="spearman")
str(cor_results)
cor_results$p.value
cor_results$estimate

## Now, we want to compare also the mean values. Have they increased or 
## decreased from 30min to 1h?
# First, we simply can look at the values
mean_samples <- apply(sample.table[,c(3,5,7,9)], 2, mean)
# Are the observed changes significant? Therefore, we can use the t-test:
t.test(sample.table$Hang_30min_Set1, sample.table$Hang_1h_Set1)
t.test(sample.table$Biof_30min_Set1, sample.table$Biof_1h_Set1)
# However, the data is not normally distributed. 
# This could be tested by using one of the normality tests,
# like the Kolmogorov-Smirnov test. Here, we specify our test distribution 
# (the mean and sd) from the empirical values
ks.test(sample.table$Hang_30min_Set1, "pnorm", mean(sample.table$Hang_30min_Set1), sd=sd(sample.table$Hang_30min_Set1) )
# Hence, the results from the t-test are not reliable. Therefore, we can 
# refer to non-parametric alternative, the Wilcoxon rank-sum test:
wilcox.test(sample.table$Hang_30min_Set1, sample.table$Hang_1h_Set1)
wilcox.test(sample.table$Biof_30min_Set1, sample.table$Biof_1h_Set1)

# We still have some additional information in our data, which can be
# used for testing, namely the pairing of the data (genes). Instead of comparing 
# simply the mean values, we could also look at the differences between each pair
# and test if they are different from zero. In our case, this would mean that some change
# has occured over time. This test can be performed simply by adding the paramter
# paired=T
wilcox.test(sample.table$Hang_30min_Set1, sample.table$Hang_1h_Set1, paired=T)
wilcox.test(sample.table$Biof_30min_Set1, sample.table$Biof_1h_Set1, paired=T)



## Now, we use the data from "Auswertung Real-time Gesamt.csv" again. Make sure
## that these data are available in your R workspace:
auswertung <- read.csv("Auswertung Real-time Gesamt.csv", sep=";", as.is = TRUE)

# Before proceeding to the analysis, we exclude the control group data
auswertung_sub <- auswertung[auswertung$Strain != "Kontrolle (Wasser)", ]
# The distributions of CT in the three strains can be compared with boxplots
boxplot(auswertung_sub$CT ~ auswertung_sub$Strain)
# The strains should also be named differently
names_strains <- c("YHUM1694", "YHUM1700", "YHUM1900")
boxplot(auswertung_sub$CT ~ auswertung_sub$Strain, col="gold", names=names_strains)

# To display the single CT values, we can use a barplot 
# The row-based organization of the CT values must the transformed to contingency table  
CT_table <- xtabs(auswertung_sub$CT ~ auswertung_sub$Strain + auswertung_sub$Gene + auswertung_sub$Replicate)
# We get two tables, for replicate1 and replicate2. These can be addressed by []
CT_table[,,1]
# Caclulate the average between both tables
CT_table.replicate <- (CT_table[,,1] + CT_table[,,2] ) / 2
# Make a barplot. The grouping will be stacked
barplot(CT_table.replicate)
# Add paramter besides=T to get a dodged version. Besides, we change colors and 
# add a legend
barplot(CT_table.replicate, beside=T, col=heat.colors(3), legend=names_strains)
# Finally, legend position might be customized
barplot(CT_table.replicate, beside=T, col=heat.colors(3), legend=names_strains,
        args.legend=list(x=30, y=28))
# To change the grouping, simply transpose the table. We also adapt colors, names and legend
barplot(t(CT_table.replicate), beside=T, col=rainbow(8), names=names_strains, legend=T,
        args.legend=list(x="bottomright"))
       
## To test differences in the mean value of CT among the strains, t-test can be used
# First, we select two different strains
auswertung_WT_tec1 <- auswertung[auswertung$Strain == "WT (YHUM1694)" | auswertung$Strain == "tec1? (YHUM1700)", ]
# T-test can now be performed. Note that we have a different data structure as in the example above
# Therefore, we can use the formula notation "~" instead of "," within t.test() function
t.test(auswertung_WT_tec1$CT ~ auswertung_WT_tec1$Strain)
# Are the data normal distributed
ks.test(auswertung_WT_tec1$CT, "pnorm", mean=mean(auswertung_WT_tec1$CT), sd=sd(auswertung_WT_tec1$CT))

## anova analysis can be used to detect differences among the groups
# First, we write down the model to be estimated
anova_model <- aov(auswertung_sub$CT ~ auswertung_sub$Strain)
# With summary, we get the results of our estimation
summary(anova_model)
# We can estimate more complex models
summary(aov(auswertung_sub$CT ~ auswertung_sub$Strain + auswertung_sub$Gene))
summary(aov(auswertung_sub$CT ~ auswertung_sub$Strain + auswertung_sub$Gene + auswertung_sub$Strain * auswertung_sub$Gene))


## R provides many different types of plots. Here are some examples

## Line charts can be generated from a scatterplot using type="l" 
## and some index variable (e.g., 1:10 or 1999:2008) for the x-axis.
x <- 1:10
y <- rnorm(10, 3, 2)
plot(x,y, type="l", lwd=2, col=2)




## There are also many interactive visualization techniques
# This example shows the googleVis implementation (you need R3.0 to run this code)
install.packages("googleVis") 
library(googleVis)
demo(WorldBank)


## Curve fitting in R
# http://www.itc.nl/~rossiter/teach/R/R_CurveFit.pdf
# http://davetang.org/muse/2013/05/09/on-curve-fitting/
# http://www.walkingrandomly.com/?p=5254


######################################
# Section 5: Selected packages for biology
######################################

## In this section we will discover some packages useful for biology/genetics.
## First, genoPlotR is a package "to produce reproducible, publication-grade 
## graphics of gene and genome maps. It allows the user to read from usual format such as 
## protein table files and blast results, as well as home-made tabular files."
## Secondly, we explore the Bioconductor project, which provides many packages for
## biostatistics and the analysis of high-throughput genomic data. As an example,
## we use the affy package to read and analyse probe level data (.CEL files). 


## genoPlotR
# Install the package on your computer 
install.packages("genoPlotR")
# Now, package can be loaded
library(genoPlotR)

## First, download two gene sequences from http://www.ncbi.nlm.nih.gov/
## Display the data as GenBank (full) and "send" them to a file
## For the example, we use the Genes NC_005956 and NC_005955
## Additionally, have process and download the comparisons between the two sequences.
## This can be done on this page: http://blast.ncbi.nlm.nih.gov/Blast.cgi and then 
## going to nucleodite_blast and entering both sequence numbers. After "Blasting", 
## the comparison can be downloaded as "Hit Table(text)"
## All files for this example are stored in the folder "gene_data"

# Read the gene sequence data, which are in a subfolder of the working directory
gene1 <- read_dna_seg_from_file("gene_data/NC_005956.gb")
gene2 <- read_dna_seg_from_file("gene_data/NC_005955.gb")
# To read  gene data from different sources and in different formats, see also
?read_dna_seg_from_file
# Now, read the comparison file of both gene sequences
comp <- read_comparison_from_blast("gene_data/V2UK5JJ0114-Alignment.txt")

# The package contains one core function to plot the genes and comparisons
# In this function, you input your genes and comparisons as list objects
# It might some computational time if the entire sequence is plotted
plot_gene_map(dna_segs=list(gene1,gene2), comparisons=list(comp))
# To plot just a sample of the sequence, define the limits first
xlims <- list(c(1, 50000), c(1, 50000))
# The same plot using the defined limits and some further plotting options
plot_gene_map(dna_segs = list(gene1,gene2), comparisons = list(comp),
              xlims = xlims, main = "Gene_1 vs Gene_2, comparison of the first 50 kb",
              gene_type = "side_blocks", dna_seg_scale = TRUE, scale = FALSE)
# You find a lot of paramters to customize your plot
?plot_gene_map

# For a more detailed description on the data structure, plotting options, etc.
# read carefully the package Vignette: 
# http://genoplotr.r-forge.r-project.org/pdfs/genoPlotR.pdf


## Bioconductor
## Bioconductor is collection of packages for the analysis of high-throughput genomic data

# To install bioconductor, see also guide on http://www.bioconductor.org/install/
# Install the basic packages
source("http://bioconductor.org/biocLite.R")
biocLite()
# This takes some time. You are also asked to update many of your pre-installed packages 
# (not all of them are part of bioconductor). 

# Then you can install the specific packages you want to use. For an overview on the packages
# see http://www.bioconductor.org/packages/release/BiocViews.html#___Software 
# and http://www.bioconductor.org/help/workflows/  
# Here, we want to use the packages "affy", which contains functions for exploratory 
# oligonucleotide array analysis. 
# The installation is different from the installation of packages from R-Cran.
biocLite("affy")
# For this package, we have to install another package in addition
biocLite("tkWidgets")
# Load package affy
library(affy)
# After loading your package, a lot of functions will be masked from standard R

# Now, we can read probe level data (.CEL files) using the interactive screen 
# Let us select all four .CEL files
Data <- ReadAffy(widget=TRUE)
# To check for the quality, we explore the data first
# Simply calling the data already provides a summary (instead of displaying the original data)
Data
# Of course, one might plot the cel data as images to detect some spatial patterns
# To plot all four images together, adjust the plot device
par(mfrow=c(2,2))
# Now, plot the Cel data
image(Data)
# Don't forget to re-store your plot device 
par(mfrow=c(1,1))

# Data can be also visualized using histograms or boxplots
hist(Data)
hist(Data[,1:2])
boxplot(Data, col=2:5)

# MAplot offers pairwise graphical comparison of intensity data
MAplot(Data, pairs=TRUE, plot.method="smoothScatter")

# The degradation can be visually analyzed
# Calculate the degradation and print the results
deg <- AffyRNAdeg(Data)
summaryAffyRNAdeg(deg)
# Plot the degradation curves
plotAffyRNAdeg(deg)

# From the probe level data (CEL), one often wants to get the expression measures
# Get the expression measures without any transformations
eset <- rma(Data)
# write the measures to a file 
write.exprs(eset, file="mydata.txt")
# Get the expression measures performing some normalizations and corrections 
eset <- expresso(Data, normalize.method="qspline",
                 bgcorrect.method="rma",pmcorrect.method="pmonly",
                 summary.method="liwong")

## For further examples and detailed descriptions, see the Vignette of this package
## http://www.bioconductor.org/packages/release/bioc/vignettes/affy/inst/doc/affy.pdf


######################################
# Section 4: Automated file processing
######################################

## In this section you will learn to process several files at once. In our
## example, you will process 100 text files. You will plot the observed
## intensities along the X axis in a simple line plot and save it to a file with
## a name that corresponds to the name of the text file. You will also produce a
## CSV file containing the mean non-saturated intensities for each file.
## The input data is in the folder "Textfiles, Mannobi., 05.11.13".

# Reading and skipping a few lines
# Putting that in a dataframe
# Time, title and y values change so focus on those
# Variability, % == 999.999, average, plot for each file. How do I write the
# outputs to a file now?

# Change into the directory containing the data
setwd("Mannobi")
# Get a list of all files in that directory
files <- list.files()
# Initialize a data frame of mean intensities with 1 row
mean.intensities <- data.frame(1)
# Loop over all files
for (file in files) {
  # Read each file, skipping the first 18 lines
  read.file <- read.table(file, skip = 18)
  # Assume that observations of 999.999 mean that the sensor was oversaturated
  read.file$V2[read.file$V2 == 999.999] <- NA
  # Construct a file name to save to by concatenating ".png" to the end of the
  # text file's name
  filename <- paste(file, ".png", sep = "")
  # Open the file as a PNG so you can save the plot to it
  png(filename)
  # Create the plot from the two columns in the table we read
  plot(read.file$V1, read.file$V2,
    type = "l",  # Plot the observations as a line
    xlab = "Nanometers",  # Label X axis
    ylab = "Intensity")  # Label Y axis
  dev.off()  # Save and close plot file
  # Record mean intensities to the data frame
  mean.intensities[[file]] <- mean(read.file$V2, na.rm = TRUE)
}
# Save the mean intensity data frame to disk
write.csv(mean.intensities, "mean_intensities.csv")
