---
title: "vars Impulse response fucntion with ggplot2"
author: "Matthieu"
date: "January 10, 2019"
output:
  html_document:
    keep_md: yes
---


## Objective

using package `vars`, the function `irf()` produces Impulse Response Function. I wanted to plot them with `ggplot2`, instead of the built-in `plot()` function. This script provides two functions:

* `as.data.frame.varirf()` converts output intoa  standard long data.frame
* `irfplot()` : do the whole conversion and plotting

## Example



```r
library(tidyverse)
library(devtools)
source_url("https://raw.githubusercontent.com/MatthieuStigler/Misconometrics/master/irf_for_ggplot/as.data.frame.varirf.R")

## standard VAR analysis
library(vars)
data(Canada)
var.2c <- VAR(Canada, p = 2, type = "const")
irf_1 <- irf(var.2c, runs = 1000)
irf_1_rw <- irf(var.2c, impulse = "rw", response =  "rw", boot = TRUE, runs = 500)
```


Plot now all


```r
irfplot(irf_1)
```

![](as.data.frame.varirf_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

Plot when called with just 1:


```r
irfplot(irf_1_rw)
```

![](as.data.frame.varirf_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

Check various options of the as.data.frame method:


```r
as.data.frame(x=irf_1, format = "long") %>%  as_tibble()
```

```
## # A tibble: 528 x 5
##    type  impulse n_ahead response  value
##    <chr> <chr>     <int> <chr>     <dbl>
##  1 IRF   e             1 e        0.363 
##  2 IRF   e             2 e        0.548 
##  3 IRF   e             3 e        0.618 
##  4 IRF   e             4 e        0.611 
##  5 IRF   e             5 e        0.552 
##  6 IRF   e             6 e        0.461 
##  7 IRF   e             7 e        0.354 
##  8 IRF   e             8 e        0.244 
##  9 IRF   e             9 e        0.139 
## 10 IRF   e            10 e        0.0449
## # … with 518 more rows
```

```r
as.data.frame(x=irf_1, format = "wide_resp") %>%  as_tibble()
```

```
## # A tibble: 132 x 7
##    type  impulse n_ahead      e     prod       rw        U
##    <chr> <chr>     <int>  <dbl>    <dbl>    <dbl>    <dbl>
##  1 IRF   e             1 0.363  -0.0206  -0.116   -0.190  
##  2 IRF   e             2 0.548  -0.00120 -0.202   -0.329  
##  3 IRF   e             3 0.618   0.0148  -0.180   -0.369  
##  4 IRF   e             4 0.611  -0.0216  -0.100   -0.353  
##  5 IRF   e             5 0.552  -0.0849   0.00805 -0.301  
##  6 IRF   e             6 0.461  -0.156    0.127   -0.230  
##  7 IRF   e             7 0.354  -0.221    0.242   -0.152  
##  8 IRF   e             8 0.244  -0.275    0.344   -0.0752 
##  9 IRF   e             9 0.139  -0.313    0.427   -0.00584
## 10 IRF   e            10 0.0449 -0.335    0.489    0.0534 
## # … with 122 more rows
```

```r
as.data.frame(x=irf_1, format = "wide_type") %>%  as_tibble()
```

```
## # A tibble: 176 x 6
##    impulse n_ahead response value_IRF value_Lower value_Upper
##    <chr>     <int> <chr>        <dbl>       <dbl>       <dbl>
##  1 e             1 e           0.363       0.284        0.395
##  2 e             2 e           0.548       0.371        0.602
##  3 e             3 e           0.618       0.351        0.710
##  4 e             4 e           0.611       0.267        0.736
##  5 e             5 e           0.552       0.165        0.717
##  6 e             6 e           0.461       0.0319       0.664
##  7 e             7 e           0.354      -0.105        0.583
##  8 e             8 e           0.244      -0.226        0.509
##  9 e             9 e           0.139      -0.323        0.454
## 10 e            10 e           0.0449     -0.413        0.385
## # … with 166 more rows
```

  


