---
title: "Index insurance metrics"
author: "Matthieu"
date: "2022-12-11"
output:
  html_document:
    keep_md: yes
---



# Intro

This page documents functions used in the paper *Optimal index insurance and basis risk decomposition: an application to Kenya*.

The functions contain fast algorithms to extract the first eigenvalue of the covariance/correlation matrix when N>>T (or p>>n). Scripts are stored in: [1_0_eigenvalue_metrics.R`](https://raw.githubusercontent.com/MatthieuStigler/Misconometrics/master/index_insurance_metrics/)

## Functions description


- `idx_r2_general`: compute the individual $R^2$ for an index Yw (cf Equation 5), specifying weights `w`
- `idx_r2_optimal`: compute the individual $R^2$ based on the optimal index (cf Theorem 1, 1.1). This function is optimized for T <<N. 
- `idx_total_r2_general`: compute the total $R^2$ specifying weights `w`
- `idx_total_r2_optimal`: compute the total $R^2$ based on the optimal index

## Demo

Load the script:


```r
source("1_0_eigenvalue_metrics.R")
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

Compute the individual $R^2$ using the mean yield as index:


```r
r2_mean <- idx_r2_general(X_sim)
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

test_that("The mean of the R2 with th eoptimal vector is the same as the first eigenvalue:",
          expect_equal(mean(r2_opt),
                       eig_values_manu_first))
```

```
## Test passed 🎉
```

### Total R2

Compute the total R2: $\bar{\bar{R}}=1-\Sigma_i SSR_i/\Sigma_i SST_i$


```r
R2_total_mean <- idx_total_r2_general(X_sim)
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

Functions to get the first eigenvalue (`idx_r2_optimal` and `idx_total_r2_optimal`) have been optimized for the case where T << N (by using the "inverted cor/cov matrix" with is T x T instead of N x N), and by using `RSpectra::eigs_sym` which takes profit of the symmetry of cor/cov matrix and allows ot compute only the first eigenvalue/vector. 

The code illustrates the speed gains with N=1000 and T=30. As the simulation shows below, most of the speed gain is obtained by using the "inverted correlation matrix" with is T x T instead of N x N. 

Prep the data, create a funciton to do it manually:


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
##         fo_cross   2.986692   3.267720   3.793061   3.447777   3.725261
##      fo_no_cross  31.909704  33.057319  42.670157  35.803299  40.124958
##     manual_cross   3.325531   3.459369   4.163353   3.589930   3.810346
##  manual_no_cross 456.282802 468.088342 504.411868 477.260693 544.825107
##        max neval cld
##   10.04797    20 a  
##  153.56017    20  b 
##   12.29082    20 a  
##  597.92157    20   c
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
##             expr        min          lq        mean      median         uq
##         fo_cross    552.062    637.7365    785.1199    653.3010    721.319
##      fo_no_cross  22477.128  24184.4895  25993.4480  25473.8800  27821.875
##     manual_cross    825.092    956.0210   1126.0622    992.9075   1093.710
##  manual_no_cross 446836.443 458979.7775 499052.2690 482053.4860 537114.317
##         max neval cld
##    3072.081    20 a  
##   31854.594    20  b 
##    3375.411    20 a  
##  590969.294    20   c
```


