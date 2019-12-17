---
title: "Anova manual"
author: "Matthieu"
date: "December 17, 2019"
output:
  html_document:
    keep_md: yes
---



# Goal


Goal of this quick note is to do a manual ANOVA, and see how to relate different levels of factors. 

Facts learned:

 - When a factor is subset of another, `anova` will drop larger factor (cf region and state)
 - If factor orthogonal, variance is additive
 - With two hierarchical factors, doing an anova on the residuals of the higher factor will give same results. 

## Prepare data and models
 

```r
library(plm)
library(tidyverse)
library(matPkg)

data(Produc)
Produc <- as_tibble(Produc) %>% 
  mutate(zone = if_else(region %in% c("1", "2", "3"), "A", "B"))
```
 
 
Estimate models:


```r
tot_SSR <- var(Produc$unemp)*(nrow(Produc)-1)


lm0 <- lm(unemp~1, data=Produc)
lm1_y <- lm(unemp~year, data=Produc)
lm1_s <- lm(unemp~state, data=Produc)
lm1_r <- lm(unemp~region, data=Produc)
lm1_z <- lm(unemp~zone, data=Produc)
lm2_ys <- lm(unemp~year+state, data=Produc)
lm2_yr <- lm(unemp~year+region, data=Produc)
lm2_sr <- lm(unemp~state+region, data=Produc)
lm3 <- lm(unemp~year+state+region, data=Produc)
```
 
Estimate cumulative models

```r
lm_cum_sr <- lm(unemp~region, data=mutate(Produc, unemp=predict(lm1_s)))
```
 

## Convenience functions


```r
aov_clean <- function(x) {
  nam <- rlang::ensym(x)
  as.data.frame(anova(x)) %>% 
    rownames_to_column(var = "variable") %>% 
    rename(SSR = `Sum Sq`) %>% 
    mutate(SSR_perc = 100* SSR/sum(SSR)) %>% 
    mat_add_total_row() %>% 
    mutate(model = as.character(nam)) %>% 
    select(model, variable, Df, SSR, SSR_perc)
  
}

aov_alter <- function(x) {
  nam <- rlang::ensym(x)
  
  res <- residuals(x)
  tibble(model = as.character(nam),
         variable ="Residuals",
         SSR = var(res) * (nobs(x)-1),
         SSR_perc = 100*SSR/tot_SSR)
}
```

## Check results


```r
aov_clean(x=lm0)
```

```
##   model  variable  Df      SSR SSR_perc
## 1   lm0 Residuals 815 4064.616      100
## 2   lm0     Total 815 4064.616      100
```

```r
aov_clean(x=lm1_y)
```

```
##   model  variable  Df       SSR  SSR_perc
## 1 lm1_y      year   1  627.0127  15.42612
## 2 lm1_y Residuals 814 3437.6033  84.57388
## 3 lm1_y     Total 815 4064.6160 100.00000
```

```r
aov_clean(x=lm1_s)
```

```
##   model  variable  Df      SSR  SSR_perc
## 1 lm1_s     state  47 1298.745  31.95248
## 2 lm1_s Residuals 768 2765.871  68.04752
## 3 lm1_s     Total 815 4064.616 100.00000
```

```r
aov_clean(x=lm1_r)
```

```
##   model  variable  Df       SSR  SSR_perc
## 1 lm1_r    region   8  663.0408  16.31251
## 2 lm1_r Residuals 807 3401.5753  83.68749
## 3 lm1_r     Total 815 4064.6160 100.00000
```

```r
aov_clean(x=lm2_ys)
```

```
##    model  variable  Df       SSR  SSR_perc
## 1 lm2_ys      year   1  627.0127  15.42612
## 2 lm2_ys     state  47 1298.7454  31.95248
## 3 lm2_ys Residuals 767 2138.8578  52.62140
## 4 lm2_ys     Total 815 4064.6160 100.00000
```

```r
aov_clean(x=lm2_yr)
```

```
##    model  variable  Df       SSR  SSR_perc
## 1 lm2_yr      year   1  627.0127  15.42612
## 2 lm2_yr    region   8  663.0408  16.31251
## 3 lm2_yr Residuals 806 2774.5625  68.26137
## 4 lm2_yr     Total 815 4064.6160 100.00000
```

```r
aov_clean(x=lm2_sr)
```

```
##    model  variable  Df      SSR  SSR_perc
## 1 lm2_sr     state  47 1298.745  31.95248
## 2 lm2_sr Residuals 768 2765.871  68.04752
## 3 lm2_sr     Total 815 4064.616 100.00000
```

```r
aov_clean(x=lm3)
```

```
##   model  variable  Df       SSR  SSR_perc
## 1   lm3      year   1  627.0127  15.42612
## 2   lm3     state  47 1298.7454  31.95248
## 3   lm3 Residuals 767 2138.8578  52.62140
## 4   lm3     Total 815 4064.6160 100.00000
```

### Alternative: residuals

Check with the alter method: just var of residuals


```r
#
aov_clean(x=lm2_yr)
```

```
##    model  variable  Df       SSR  SSR_perc
## 1 lm2_yr      year   1  627.0127  15.42612
## 2 lm2_yr    region   8  663.0408  16.31251
## 3 lm2_yr Residuals 806 2774.5625  68.26137
## 4 lm2_yr     Total 815 4064.6160 100.00000
```

```r
aov_alter(x=lm2_yr)
```

```
## # A tibble: 1 x 4
##   model  variable    SSR SSR_perc
##   <chr>  <chr>     <dbl>    <dbl>
## 1 lm2_yr Residuals 2775.     68.3
```

### Decompose hierarchical ones

Check residual approach: if state gives 31%, and region 16%, can I claim that region is generating 16/31 of variance? yes


```r
aov_clean(lm1_s)
```

```
##   model  variable  Df      SSR  SSR_perc
## 1 lm1_s     state  47 1298.745  31.95248
## 2 lm1_s Residuals 768 2765.871  68.04752
## 3 lm1_s     Total 815 4064.616 100.00000
```

```r
aov_clean(lm1_r)
```

```
##   model  variable  Df       SSR  SSR_perc
## 1 lm1_r    region   8  663.0408  16.31251
## 2 lm1_r Residuals 807 3401.5753  83.68749
## 3 lm1_r     Total 815 4064.6160 100.00000
```

```r
aov_clean(lm_cum_sr)
```

```
##       model  variable  Df       SSR SSR_perc
## 1 lm_cum_sr    region   8  663.0408  51.0524
## 2 lm_cum_sr Residuals 807  635.7047  48.9476
## 3 lm_cum_sr     Total 815 1298.7454 100.0000
```

Using residuals only:


```r
cumdiff <- function(x) cumsum(diff(x))

dat <- rbind(aov_alter(lm1_z) %>% 
               mutate(variable="resid_zone"),
             aov_alter(lm1_s) %>% 
               mutate(variable="resid_state"),
             aov_alter(lm1_r)%>% 
               mutate(variable="resid_region")) %>% 
  bind_rows(tibble(model = "all",
                   variable ="resid_total",
                   SSR = tot_SSR, SSR_perc=100))

dat
```

```
## # A tibble: 4 x 4
##   model variable       SSR SSR_perc
##   <chr> <chr>        <dbl>    <dbl>
## 1 lm1_z resid_zone   4009.     98.6
## 2 lm1_s resid_state  2766.     68.0
## 3 lm1_r resid_region 3402.     83.7
## 4 all   resid_total  4065.    100
```

```r
dat[1:3,] %>% 
  mutate(var_take = 100 - SSR_perc) %>%   
  arrange(var_take) %>% 
  mutate(var_use = c(var_take[1], abs(diff(var_take))),
         check = cumsum(var_use)) %>% 
  as.data.frame()
```

```
##   model     variable      SSR SSR_perc  var_take   var_use     check
## 1 lm1_z   resid_zone 4009.177 98.63606  1.363937  1.363937  1.363937
## 2 lm1_r resid_region 3401.575 83.68749 16.312507 14.948570 16.312507
## 3 lm1_s  resid_state 2765.871 68.04752 31.952476 15.639969 31.952476
```

Check: want to have same values as in column `var_take`:


```r
rbind(aov_clean(lm1_z),
      aov_clean(lm1_r),
      aov_clean(lm1_s)) %>% 
  filter(!variable %in% c("Total", "Residuals"))
```

```
##   model variable Df        SSR  SSR_perc
## 1 lm1_z     zone  1   55.43879  1.363937
## 2 lm1_r   region  8  663.04076 16.312507
## 3 lm1_s    state 47 1298.74544 31.952476
```

