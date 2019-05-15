---
title: 'R panel plm function with big data: benchmark for speed improvement'
author: "Matthieu"
date: "September 15, 2018"
output:
  html_document:
    toc: true
    keep_md: yes
---
  


## Purpose

R package for panel data, `plm`, ie vey powerful, yet slow in presence of large data. In these scripts, I intendo to do:

1. Measure the speed of `plm` verus alternatives such aspackage `lfe` with its funciton `felm`
2. Measure the speed of the raw function Within(), taht does the FE transformation. 
3. Look for alternative function, using package `data.table`. 


## Results for now:

1. Yes, `plm` can be very slow in presence of large datasets
2. Package `lfe`, with its function `felm`, and underlying Within funciton `demeanlist`, is much faster
3. `data.table` seems to offer very interesting speed imporvements, over plm and even over felm (sometimes)

## Plots

![Graphic](figures/benchm_center_dt_lfe.png?raw=true "Title")
