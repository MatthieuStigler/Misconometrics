---
title: "Index insurance metrics"
author: "Matthieu"
date: "2022-12-11"
output:
  html_document:
    keep_md: yes
---



# Intro

This page documents functions used in the paper  Stigler, Matthieu, and Lobell, David (2023) “Optimal Index Insurance and Basis Risk Decomposition: An Application to Kenya.” _American Journal of Agricultural Economics_ 1–24. https://doi.org/10.1111/ajae.12375

The functions contain fast algorithms to extract the first eigenvalue of the covariance/correlation matrix when N>>T (or p>>n). Scripts are stored in Github: [1_0_eigenvalue_metrics.R](https://raw.githubusercontent.com/MatthieuStigler/Misconometrics/master/index_insurance_metrics/)

## Functions description


- `idx_r2_general`: compute the individual $R^2$ for an index Yw (cf Equation 5), specifying weights `w`
- `idx_r2_optimal`: compute the individual $R^2$ based on the optimal index (cf Theorem 1, 1.1). This function is optimized for T <<N. 
- `idx_total_r2_general`: compute the total $R^2$ specifying weights `w`
- `idx_total_r2_optimal`: compute the total $R^2$ based on the optimal index

## Demo

Load the script:


```r
source("1_0_eigenvalue_metrics.R")

## load the testthat package to have pretty test output:
library(testthat)

## alternatively, source the script from github:
# source_url("https://raw.githubusercontent.com/MatthieuStigler/Misconometrics/master/index_insurance_metrics/1_0_eigenvalue_metrics.R")
```


Generate some data using `idx_mvnorm_sim()`:


```r
X_sim <- idx_mvnorm_sim(n=30, lambda = c(200, rep(1, 199)), seed = 123)
dim(X_sim)
```

```
## [1]  30 200
```

The data has 30 rows (years T) and 200 columns (fields N).

Compute the individual $R^2$ using the mean yield ($w=1/N$) as index:


```r
r2_mean <- idx_r2_general(X_sim, w = rep(1/200, 200))
```

This is the same as using `lm()`:


```r
regs_manu <- lm(X_sim~rowMeans(X_sim))
r2_manu <- sapply(summary(regs_manu), \(x) x$r.squared)

test_that("Using 'idx_r2_general' is the same as using lm manually:",
          expect_equal(r2_mean, 
                       r2_manu, check.attributes=FALSE))
```

```
## Test passed 😀
```

Compute the individual $R^2$ using the optimal index with `idx_r2_optimal()`:


```r
r2_opt <- idx_r2_optimal(X_sim)
```

Do the same but manually:


```r
w_opt <- eigen(cor(X_sim))$vectors[,1] * sqrt(1/diag(cov(X_sim))) # Theorem 1, 1.1
r2_opt_manu <- idx_r2_general(X_sim, w = w_opt)
all.equal(r2_opt_manu, r2_opt)
```

```
## [1] TRUE
```

Check Theorem 1: $\bar{R}^2(w^*_{cor}) =\lambda^{cor}_1/\Sigma_i \lambda^{cor}_i$:


```r
## Compute first eigenvalue manually:
eig_values_manu <- eigen(cor(X_sim))$values
eig_values_manu_first <- eig_values_manu[1]/sum(eig_values_manu)

test_that("The mean of the R2 with the optimal vector is the same as the first eigenvalue:",
          expect_equal(mean(r2_opt),
                       eig_values_manu_first))
```

```
## Test passed 🎉
```

### Total R2

Compute the total R2: $\bar{\bar{R}}=1-\Sigma_i SSR_i/\Sigma_i SST_i$ using the mean as index:


```r
R2_total_mean <- idx_total_r2_general(X_sim, w = rep(1/200, 200))
R2_total_mean
```

```
## [1] 0.3200817
```

This is the same as doing manually with `lm()`:


```r
SSR_tot <- sum(deviance(regs_manu))
SST_tot <- sum(deviance( lm(X_sim~1)))
R2_total_mean_manual <- 1-SSR_tot/SST_tot 
  
test_that("Using 'idx_total_r2_general' is the same as using lm manually:",
          expect_equal(R2_total_mean_manual,
                       R2_total_mean))
```

```
## Test passed 🥇
```

Check theorem 2: the total $\bar{\bar{R}}^2(w^*_{cov}) =\lambda^{cov}_1/\Sigma_i \lambda^{cov}_i$:


```r
R2_total_opt <- idx_total_r2_general(X_sim, w = eigen(cov(X_sim))$vectors[,1])
R2_total_opt_2 <- idx_total_r2_optimal(X_sim)
eig_S_values_manu <- eigen(cov(X_sim))$values
eig_S_values_manu_first <- eig_S_values_manu[1]/sum(eig_S_values_manu)

test_that("The Total R2 with the optimal vector is the same as the first eigenvalue",
          expect_equal(R2_total_opt, 
                       eig_S_values_manu_first))
```

```
## Test passed 🎊
```

Check also that this is  weighted sum of the $R^2$ obtained with the optimal vector for S:


```r
r2_opt_S <- idx_r2_general(df_w = X_sim, w = eigen(cov(X_sim))$vectors[,1])
w_var <- apply(X_sim, 2, var)

test_that("The Total R2 with the optimal vector is also the weighted average of R2:",
          expect_equal(weighted.mean(r2_opt_S, w_var),
                       R2_total_opt))
```

```
## Test passed 🥳
```


## Benchmark

Functions to get the first eigenvalue (`idx_r2_optimal` and `idx_total_r2_optimal`) have been optimized for the case where T << N (by using the "inverted cor/cov matrix" with is T x T instead of N x N), and by using `RSpectra::eigs_sym` which takes profit of the symmetry of cor/cov matrix and allows to compute only the first eigenvalue/vector. 

The code illustrates the speed gains with N=1000 and T=30. As the simulation shows below, most of the speed gain is obtained by using the "inverted correlation matrix" with is T x T instead of N x N. 

Prep the data, create a function to do it manually:


```r
## testing function:
fully_manual <- function(df_w, type = c("cor", "cov"), cross=TRUE) {
    type <- match.arg(type)
    if(cross) {
      covr_fo <- switch(type, "cov" = idx_cov_cross, "cor" = idx_cor_cross)
    } else {
      covr_fo <- switch(type, "cov" = cov, "cor" = cor)
    }
    eigs_covr <- eigen(covr_fo(df_w))$value 
    eigs_covr[1]/sum(eigs_covr)
  }

## simulate big data
lambdas_big <- c(1000, rep(1, 1000-1))
lambdas_big[1]/sum(lambdas_big)
```

```
## [1] 0.5002501
```

```r
X_sim_big <- idx_mvnorm_sim(n=30, lambdas_big)
```

Estimate the first eigenvalue of the correlation matrix

```r
microbenchmark::microbenchmark(fo_cross = mean(idx_r2_optimal(X_sim_big, cross = TRUE)),
                               fo_no_cross = mean(idx_r2_optimal(X_sim_big, cross = FALSE)),
                               manual_cross = fully_manual(X_sim_big),
                               manual_no_cross = fully_manual(X_sim_big, cross = FALSE),
                               times = 20,
                               check = "equal")
```

```
## Unit: milliseconds
##             expr        min         lq       mean     median         uq
##         fo_cross   2.996057   3.191487   3.859267   3.636570   3.813746
##      fo_no_cross  32.204291  35.320998  45.238038  40.638125  43.204420
##     manual_cross   3.152854   3.587471   4.418494   3.818208   4.258894
##  manual_no_cross 458.205645 486.035717 545.635268 562.012745 587.490566
##         max neval cld
##    9.652584    20 a  
##  155.452985    20  b 
##   12.463229    20 a  
##  683.575322    20   c
```

Estimate the first eigenvalue of the covariance matrix

```r
microbenchmark::microbenchmark(fo_cross = idx_total_r2_optimal(X_sim_big, cross = TRUE),
                               fo_no_cross = idx_total_r2_optimal(X_sim_big, cross = FALSE),
                               manual_cross = fully_manual(X_sim_big, type = "cov"),
                               manual_no_cross = fully_manual(X_sim_big, cross = FALSE, type = "cov"),
                               times = 20,
                               check = "equal")
```

```
## Unit: microseconds
##             expr        min          lq        mean     median          uq
##         fo_cross    669.723    716.6715    954.4426    761.317    876.6855
##      fo_no_cross  27266.692  28392.1900  30704.4762  31225.996  32478.6295
##     manual_cross    914.739   1053.0610   1291.6271   1138.018   1203.8055
##  manual_no_cross 535457.994 544686.4115 577602.4378 568600.546 606080.9475
##         max neval cld
##    3911.919    20 a  
##   35287.872    20  b 
##    4095.436    20 a  
##  661109.095    20   c
```


