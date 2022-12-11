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

The functions contain fast algorithms to extract the firt eigenvalue of the covariance/correlation matrix when N>>T (or p>>n) 
stored in 1_0_eigenvalue_metrics.R used  https://raw.githubusercontent.com/MatthieuStigler/Misconometrics/master/index_insurance_metrics/

## Functions description


- `idx_r2_general`: compute the individual $R^2$ for an index Yw (cf Equation 5), specifying weights `w`
- `idx_r2_optim`: compute the individual $R^2$ based on the optimal index (cf Theorem 1, 1.1). This function is optimized for T <<N. 
- `idx_total_R2`: compute the total $R^2$ 

## Demo

Load the script:


```r
source("1_0_eigenvalue_metrics.R")
```


Generate some data:


```r
X_sim <- idx_mvnorm_sim(n=30, lambda = c(200, rep(1, 199)), seed = 123)
dim(X_sim)
```

```
## [1]  30 200
```

The data has 30 rows (years T) and 200 columns (fields N).

Compute the $R^2$ using the mean yield as index:


```r
r2_mean <- idx_r2_general(X_sim)
```

This is the same as using `lm()`:


```r
regs_manu <- lm(X_sim~rowMeans(X_sim))
r2_manu <- sapply(summary(regs_manu), \(x) x$r.squared)
all.equal(r2_mean, 
          r2_manu, check.attributes=FALSE)
```

```
## [1] TRUE
```

Compute the $R^2$ using the optimal index:


```r
w_opt <- eigen(cor(X_sim))$vectors[,1] * sqrt(1/diag(cov(X_sim))) # Theorem 1, 1.1
r2_opt <- idx_r2_general(X_sim, w = w_opt)
```

Check Theorem 1: $\bar{R}^2(w^*_{cor}) =\lambda^{cor}_1/\Sigma_i \lambda^{cor}_i$:


```r
eig_values_manu <- eigen(cor(X_sim))$values
eig_values_manu_first <- eig_values_manu[1]/sum(eig_values_manu)
all.equal(mean(r2_opt),
          eig_values_manu_first)
```

```
## [1] TRUE
```

### Total R2


```r
R2_total_mean <- idx_total_R2(X_sim)
R2_total_mean
```

```
## [1] 0.3200817
```

This is the same as doing manually $1-\Sigma_i SSR_i/\Sigma_i SST_i$


```r
SSR_tot <- sum(deviance(regs_manu))
SST_tot <- sum(deviance( lm(X_sim~1)))
all.equal(1-SSR_tot/SST_tot,
          R2_total_mean)
```

```
## [1] TRUE
```

Check theorem 2: the total $\bar{\bar{R}}^2(w^*_{cov}) =\lambda^{cov}_1/\Sigma_i \lambda^{cov}_i$:


```r
R2_total_opt <- idx_total_R2(X_sim, w = eigen(cov(X_sim))$vectors[,1])
eig_S_values_manu <- eigen(cov(X_sim))$values
eig_S_values_manu_first <- eig_S_values_manu[1]/sum(eig_S_values_manu)

all.equal(R2_total_opt, 
          eig_S_values_manu_first)
```

```
## [1] TRUE
```

Check also that this is  weighted sum of the $R^2$ obtained with the optimal vector for S:


```r
r2_opt_S <- idx_r2_general(X_sim, w = eigen(cov(X_sim))$vectors[,1])
w_var <- apply(X_sim, 2, var)

all.equal(weighted.mean(r2_opt_S, w_var),
          R2_total_opt)
```

```
## [1] TRUE
```


## Benchmark



