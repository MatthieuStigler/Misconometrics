---
title: 'Panel in R: interact time FE with space units'
author: "Matthieu"
abstract: "Goal is to add interactions with a time effects and some group, such as region, county etc"
date: "May 15, 2019"
output:
  html_document:
    toc: true
    keep_md: yes
---




```r
library(plm)
library(lfe)
library(broom)
library(tidyverse)
```


# Summary, refs

Discussed on Stack Overflow:

- [r-plm-and-lm-fixed-effects](https://stackoverflow.com/questions/43636724/r-plm-and-lm-fixed-effects)
- [fixed-effects-regression-with-state-specific-trends](https://stackoverflow.com/questions/34232834/fixed-effects-regression-with-state-specific-trends)

There is a quite complicated story about counting the degrees of freedom: 

- standard FE2: just need to remove one df
- year-region: here interacting time with 9 regions implies we should remove (i.e. normalize the parameters) $n_R$ -1=8 here!

`lm()` and `felm()` give same answer (-8), `plm()` not. Getting the fixed-effects in `felm()` shows indeed 8 FE are set to 0. In `plm()`, noen of them is zero, so it is not clear what the transformation is?



State specific trend: need to convert to integer


```r
fe1_Ttrend_plm <- plm(log(gsp) ~  state : as.integer(year), data = produc_plm)
```


# Preliminary: Load data


```r
data("Produc", package = "plm")
```

Add region:


```r
Produc$year_region <- paste(Produc$year, Produc$region, sep = "_")
```


# With package plm

## Simple reg


```r
zz <- plm(gsp ~ pcap, data = Produc, index = c("state","year"), effect = "twoway")
summary(zz)
```

```
## Twoways effects Within Model
## 
## Call:
## plm(formula = gsp ~ pcap, data = Produc, effect = "twoway", index = c("state", 
##     "year"))
## 
## Balanced Panel: n = 48, T = 17, N = 816
## 
## Residuals:
##      Min.   1st Qu.    Median   3rd Qu.      Max. 
## -69163.02  -1934.11    363.24   2328.71 101946.55 
## 
## Coefficients:
##      Estimate Std. Error t-value  Pr(>|t|)    
## pcap  2.23051    0.15994  13.946 < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Total Sum of Squares:    8.6552e+10
## Residual Sum of Squares: 6.8748e+10
## R-Squared:      0.20571
## Adj. R-Squared: 0.13802
## F-statistic: 194.494 on 1 and 751 DF, p-value: < 2.22e-16
```

## Reg with year_region


```r
zz2 <- plm(gsp ~ pcap, data = Produc, index = c("state","year_region"), effect = "twoway")
zz
```

```
## 
## Model Formula: gsp ~ pcap
## 
## Coefficients:
##   pcap 
## 2.2305
```

```r
zz2
```

```
## 
## Model Formula: gsp ~ pcap
## 
## Coefficients:
##   pcap 
## 2.2845
```


Check FE2:


```r
FE_1 <- fixef(zz, effect = "time")
FE_2 <- fixef(zz2, effect = "time")
FE_df_1 <- tibble(fe = names(FE_1), value = FE_1) %>% 
  mutate(year = as.integer(fe))
FE_df_2 <- tibble(fe = names(FE_2), value = FE_2) %>% 
  separate(fe, c("year", "region"), remove = FALSE, convert = TRUE)
FE_df_2
```

```
## # A tibble: 153 x 4
##    fe      year region   value
##    <chr>  <int>  <int>   <dbl>
##  1 1970_1  1970      1   2520.
##  2 1970_2  1970      2  14314.
##  3 1970_3  1970      3   7649.
##  4 1970_4  1970      4  -2634.
##  5 1970_5  1970      5   1921.
##  6 1970_6  1970      6  -7259.
##  7 1970_7  1970      7  15484.
##  8 1970_8  1970      8  -1846.
##  9 1970_9  1970      9 -17478.
## 10 1971_1  1971      1   1713.
## # … with 143 more rows
```


Average FE


```r
## ave
FE_df_2_ave <- FE_df_2 %>% 
  group_by(year) %>% 
  summarise(value = mean(value)) %>% 
  ungroup()

FE_df_both <- rbind(select(FE_df_1, -fe),
                    FE_df_2_ave) %>% 
  mutate(type = rep(c("FE_1", "FE_2"), each = nrow(FE_df_1)))

FE_df_both
```

```
## # A tibble: 34 x 3
##    value  year type 
##    <dbl> <int> <chr>
##  1 2353.  1970 FE_1 
##  2 1640.  1971 FE_1 
##  3 2663.  1972 FE_1 
##  4 4164.  1973 FE_1 
##  5 2437.  1974 FE_1 
##  6  149.  1975 FE_1 
##  7 1537.  1976 FE_1 
##  8 3314.  1977 FE_1 
##  9 5733.  1978 FE_1 
## 10 6460.  1979 FE_1 
## # … with 24 more rows
```


Plot it

```r
## plot
FE_df_both %>% 
  ggplot(aes(x = year, y = value, colour = type))+
  geom_point()
```

![](README_files/figure-html/unnamed-chunk-9-1.png)<!-- -->




# With lm


```r
lm1 <- lm(gsp ~ pcap -1 + state + factor(year), data = Produc)
lm2 <- lm(gsp ~ pcap -1 + state + year_region, data = Produc)
```

Compare:


```r
all.equal(tidy(lm1) %>% 
            filter(term=="pcap") %>% 
            as.data.frame(), 
          tidy(zz) %>%
            as.data.frame())
```

```
## [1] TRUE
```

```r
all.equal(tidy(lm2) %>% 
            filter(term=="pcap") %>% 
            as.data.frame(), 
          tidy(zz2) %>%
            as.data.frame())
```

```
## [1] "Component \"std.error\": Mean relative difference: 0.00648305" 
## [2] "Component \"statistic\": Mean relative difference: 0.006441291"
```



# With package lfe (felm  function)

### Run simple


```r
## year
felm <- felm(gsp ~ pcap|state+year, data = Produc)
all.equal(tidy(felm) %>% as.data.frame(), 
          tidy(zz)%>% as.data.frame())
```

```
## [1] TRUE
```

### Run with year_region


```r
## year region
felm2 <- felm(gsp ~ pcap|state+year_region, data = Produc, exactDOF="rM")

all.equal(coef(felm2), coef(zz2))
```

```
## [1] TRUE
```

```r
## but SE not the same!
all.equal(tidy(felm2) %>% as.data.frame(), 
          tidy(zz2)%>% as.data.frame())
```

```
## [1] "Component \"std.error\": Mean relative difference: 0.00648305" 
## [2] "Component \"statistic\": Mean relative difference: 0.006441291"
```

```r
## but same as lm!?
all.equal(tidy(lm2) %>% 
            filter(term=="pcap") %>% 
            as.data.frame(), 
          tidy(felm2) %>%
            as.data.frame())
```

```
## [1] TRUE
```

Why not same? Df are different! felm has 8 more!


```r
c(felm$df.residual, zz$df.residual)
```

```
## [1] 751 751
```

```r
c(felm=felm2$df.residual, plm=zz2$df.residual)
```

```
## felm  plm 
##  623  615
```

```r
diff(c(felm2$df.residual, zz2$df.residual))
```

```
## [1] -8
```


```r
Produc %>% 
  summarise_at(c("state", "year", "year_region", "region"), n_distinct)
```

```
##   state year year_region region
## 1    48   17         153      9
```

```r
## so technically, we estimate 48+153 = 201
nrow(Produc) - 201
```

```
## [1] 615
```

```r
## This is the df of plm!

## But many FEs are zero!
getfe(felm) %>% 
  count(effect_zero=effect==0)
```

```
## # A tibble: 2 x 2
##   effect_zero     n
##   <lgl>       <int>
## 1 FALSE          64
## 2 TRUE            1
```

```r
getfe(felm2) %>% 
  count(effect_zero=effect==0)
```

```
## # A tibble: 2 x 2
##   effect_zero     n
##   <lgl>       <int>
## 1 FALSE         192
## 2 TRUE            9
```

```r
## plm 
tibble(n_ind =length(fixef(zz, effect = "individual")),
       n_time =length(fixef(zz, effect = "time")), 
       tot = n_ind+n_time)
```

```
## # A tibble: 1 x 3
##   n_ind n_time   tot
##   <int>  <int> <int>
## 1    48     17    65
```

```r
tibble(n_ind =length(fixef(zz2, effect = "individual")),
       n_time =length(fixef(zz2, effect = "time")), 
       tot = n_ind+n_time)
```

```
## # A tibble: 1 x 3
##   n_ind n_time   tot
##   <int>  <int> <int>
## 1    48    153   201
```



Do with factors:


```r
Produc2 <- Produc %>% 
  mutate(region = as.factor(region),
         year = as.factor(year))
felm2b <- felm(gsp ~ pcap|state+year:region, data = Produc2)

all.equal(tidy(felm2) %>% as.data.frame(), 
          tidy(felm2b)%>% as.data.frame())
```

```
## [1] TRUE
```



