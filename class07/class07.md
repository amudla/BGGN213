Class07: R function and packages
================
Anusorn Mudla
4/24/2019

### Today lecture is more on writing R function and will introduce R packages

``` r
source("http://tinyurl.com/rescale-R")
```

Test the **rescale()** function

``` r
rescale(c(1:10))
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
x <- c(1,2,NA)
y <- c(1,3,NA)
# check for NA in the vector
is.na(x)
```

    ## [1] FALSE FALSE  TRUE

``` r
is.na(y)
```

    ## [1] FALSE FALSE  TRUE

``` r
# check if both vector has NA at the same position
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE

``` r
# check how many NA in both vector
sum(is.na(x) & is.na(y))
```

    ## [1] 1

``` r
# turn the snipple into a fuction
both_na <- function(x,y) {
  sum( is.na(x) & is.na(y) )
}
both_na(x,y)
```

    ## [1] 1

``` r
# check the length of the vector before operate on them
both_na2 <- function(x, y) {
  ## Check for NA elements in both input vectors and don't allow re-cycling 
  if(length(x) != length(y)) {
    stop("Input x and y should be vectors of the same length", call.=FALSE)
  }
  sum( is.na(x) & is.na(y) )
}
```

``` r
#adding some messages to the end of function
both_na3(x,y)
```

    ## Found 1 NA's at position(s):3

    ## $number
    ## [1] 1
    ## 
    ## $which
    ## [1] 3

Writing the grade function
--------------------------

``` r
#student 1
std1 <- c(100,100,100,100,100,100,100,90)
std2 <- c(100,NA,90,90,90,90,97,80)

index.remove <- which.min(std1)
std1.new <- std1[-index.remove]
std1.new
```

    ## [1] 100 100 100 100 100 100 100

``` r
mean(std1.new)
```

    ## [1] 100

``` r
grade <- function(score) {
  #remove the na from the score
  score.noNA <- score[!is.na(score)]
  index.remove <- which.min(score.noNA)
  score.new <- score.noNA[-index.remove]
  mean(score.new)
}
```

``` r
# for not turning in homework will not have negative effect
grade2 <- function(score) {
  (sum(score,na.rm = TRUE)-min(score, na.rm = TRUE)) / (length(score[!is.na(score)])-1)
}
```

``` r
# for not turing 
grade3 <- function(score) {
  (sum(score,na.rm = TRUE)-min(score, na.rm = TRUE)) / (length(score)-1)
}
```

``` r
## apply the function to the csv file
url <- "http://tinyurl.com/gradeinput"
class_data <- read.csv(url,row.names = 1)

grade2(class_data[1,])
```

    ## [1] 91.75

``` r
# use the apply function to each row. 1 = by row and 2 = by column
average <- apply(class_data,1,grade3)
# sort the average to rank student
sort(average, decreasing = TRUE)
```

    ##  student-7  student-8 student-13  student-1 student-12 student-16 
    ##      94.00      93.75      92.25      91.75      91.75      89.50 
    ##  student-6  student-5 student-17  student-9 student-14 student-11 
    ##      89.00      88.25      88.00      87.75      87.75      86.00 
    ##  student-3 student-19 student-20  student-2 student-18  student-4 
    ##      84.25      82.75      82.75      82.50      72.75      66.00 
    ## student-15 student-10 
    ##      62.50      61.00

``` r
## add the aveage column to the class data
class_data$average <- apply(class_data,1,grade3)
# how to sort the data frame by any particular column
```

### find intersaction between two vectors

``` r
x <- df1$IDs
y <- df2$IDs
intersect(x,y)
```

    ## [1] "gene2" "gene3"

``` r
x %in% y
```

    ## [1] FALSE  TRUE  TRUE

``` r
gene_intersect <- function(x, y) {
  cbind(x[x %in% y], y[y %in% x])
}
```

``` r
merge(df1,df2,by= "IDs")
```

    ##     IDs exp.x exp.y
    ## 1 gene2     1    -2
    ## 2 gene3     1     1
