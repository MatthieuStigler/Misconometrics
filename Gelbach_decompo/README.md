---
title: 'R script for Gelbach: `When Do Covariates Matter? And Which Ones, and How Much?`'
author: "Matthieu"
date: "September 4, 2018"
output:
  html_document:
    keep_md: yes
---



## Purpose

This code implements the basics of the paper: Gelbach (2016) *When Do Covariates Matter? And Which Ones, and How Much?* Journal of Labor Economics, 2016, vol. 34(2)

It does just do the decomposition, without inference. It does not do the grouping either: in the standard OLS case, R can do a multivariate OLS, so computing the auxiliary regressions for one or multiple LHS variables is done efficiently. User can simply use `tidyverse::group_by` tools to do the grouping. 


## Use

The function is: `dec_covar(object = , var_main)`. `var_main` is the variable of interest, as character string.

The function works for `object` of class:
 * OLS regressions computed with `lm()` 
 * panel computed with  `plm()` or `felm()` from package *lfe*


## Illustration


Source the script:


```r
library(devtools)
source_url("https://raw.githubusercontent.com/MatthieuStigler/Misconometrics/master/Gelbach_decompo/dec_covar.R")
```

```
## SHA-1 hash of file is cc7ac4d4dcdb89a4df6c70ce66999c555dbd62ad
```


```r
## run a OLS with built-in data freeny
model_full_1 <- lm(y ~ lag.quarterly.revenue + price.index + income.level + market.potential, data=freeny)

## use decomposition:
dec_covar(object = model_full_1, var_main = "lag.quarterly.revenue")
```

```
##                               variable      gamma     beta_K     delta
## 1                          price.index -0.4182537 -0.7542401 0.3154637
## 2                         income.level  0.3747902  0.7674609 0.2876368
## 3                     market.potential  0.2038996  1.3305577 0.2713002
##                                  Total         NA         NA 0.8744008
## lag.quarterly.revenue            Check         NA         NA 0.8744008
##                           perc
## 1                     36.07770
## 2                     32.89531
## 3                     31.02699
##                             NA
## lag.quarterly.revenue       NA
```

