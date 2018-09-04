#' ---
#' title: "Implement Gelbach dec_covarsition"
#' author: "Matthieu"
#' abstract: Replicate method of Gelbach, link: https://www.journals.uchicago.edu/doi/full/10.1086/683668
#' stata file: http://fmwww.bc.edu/repec/bocode/b/b1x2.ado, 
#' date: "2018-09-03"
#' output:
#'    pdf_document:
#'       toc: true
#' ---


#+ echo=FALSE, message=FALSE, results="hide"
if(!interactive()){
  library(knitr)
  knitr::opts_knit$set(root.dir = '..' )
  
}

library(tidyverse)


#'###############################
#'## First test: OVB formula
################################


freeny_4 <- freeny %>%
  mutate(income.level = income.level + market.potential) %>%
  select(-market.potential)

coef(lm(freeny))
coef(lm(freeny_4))


#'###############################
#'## First test: OVB formula
################################

freeny_tb <- as_tibble(freeny)

reg_k2 <- lm(freeny[, 1:3])
reg_k1 <- lm(freeny[, 1:2])

## aux reg:
reg_aux_a <- lm(freeny[, 3:2])
reg_aux_b <- lm(freeny[, 3:2])

##
coef(reg_k2)
coef(reg_k1)

## diff:
coef(reg_k2)[2] - coef(reg_k1)[2]

## relation
coef(reg_aux_a)[2] * coef(reg_k2)[3]

#'###############################
#'## test 2: on freeny
#'###############################

reg_full <- lm(y ~ lag.quarterly.revenue + price.index + income.level + market.potential, data=freeny)
reg_base <- lm(y ~ lag.quarterly.revenue, data=freeny)

freeny_M <- as.matrix(freeny)
reg_aux <- lm(freeny_M[, -c(1:2)]~freeny_M[, 2])

## diff
coef(reg_base)[2] - coef(reg_full)[2]

## dec_covar:
coef(reg_aux)[2,] %*% coef(reg_full)[-c(1,2)]


## in df
X <- data_frame(variable  = names(coef(reg_full))[-c(1,2)],
                gamma = coef(reg_aux)[2,],
                beta_K = coef(reg_full)[-c(1,2)],
                delta = gamma * beta_K,
                perc = 100 * delta/sum(delta)) %>%
  bind_rows(data_frame(variable = "Total", 
                       gamma = NA, 
                       beta_K =NA,
                       delta = sum(.$delta)))
X
  # bind_rows(summarise_all(., funs(if(is.numeric(.)) sum(.) else "Total")))

## 
x <- "lag.quarterly.revenue"
f2 <- reformulate(x, ff[[2]])


parse(eval)
update(reg_full, . ~ lag.quarterly.revenue)
update(reg_full, reformulate(x, formula(reg_full)[[2]]) )

update(reg_full, cbind(price.index, income.level, market.potential) ~ lag.quarterly.revenue)

string_formula <- sprintf("cbind(%s) ~ .", toString(names(iris)[1:2]))
string_formula <- sprintf("cbind(%s) ~ %s", toString(names(iris)[1:2]), x)


################################
#'## From lm
################################



dec_covar <- function(object, var_main) {
  
  
  ## Y:
  formu_obj <- formula(object)
  response_var <- deparse(formu_obj[[2]])
  var_all <- attr(terms(formu_obj), "term.labels")
  if(!var_main %in% var_all) stop(paste("Variable: '", var_main, "' not found in columns: '", 
                                        paste(var_all, collapse = ", "), "'", sep=""))
  var_other <- var_all[var_all!=var_main]
  
  XY <- model.frame(object)
  y <- XY[, response_var]
  x_main <- XY[, var_main]
  x_others <- XY[, var_other]
  
  
  ## reg base
  reg_base <- lm(y~x_main)
  
  ## aux regs
  reg_aux <- lm(as.matrix(x_others) ~ x_main)
  
  ## RES
  X <- data.frame(variable  = var_other,
                  gamma = coef(reg_aux)["x_main",],
                  beta_K = coef(object)[var_other],row.names = NULL)
  X$delta = X$gamma *  X$beta_K
  X$perc <- 100 * X$delta/sum(X$delta)
  X <- rbind(X, data.frame(variable = "Total", 
                           gamma = NA, 
                           beta_K =NA,
                           delta = sum(X$delta),
                           perc=NA))
  
  X
  
}

dec_covar <- function(object, var_main) {
  
  
  ## Y:
  formu_obj <- formula(object)
  response_var <- deparse(formu_obj[[2]])
  # reg_var <- deparse(formu_obj[[3]])
  var_all <- attr(terms(formu_obj), "term.labels")
  var_other <- var_all[var_all!=var_main]
  
  XY <- model.frame(object)
  y <- XY[, response_var]
  x_main <- XY[, var_main]
  x_others <- XY[, var_other]
  
  
  ## reg base
  reg_base <- update(object, reformulate(var_main, formula(object)[[2]]) )
  
  ## aux regs
  if(inherits(object, "lm")) {
    string_formula <- sprintf("cbind(%s) ~ %s", toString(var_other), var_main)
    reg_aux <- update(object, as.formula(string_formula))  
    gamma <- coef(reg_aux)[var_main,]
  } else {
    string_formula2 <- sprintf("%s ~ %s", var_other, var_main)
    reg_aux_all <- lapply(string_formula2, function(x) update(object, as.formula(x))) 
    coefs <- lapply(reg_aux_all, coef)  
    gamma <- unlist(coefs)
  }
  
  
  ## separately
  
  
  ## RES
  X <- data.frame(variable  = var_other,
                  gamma = gamma,
                  beta_K = coef(object)[var_other], row.names = NULL)
  X$delta = X$gamma *  X$beta_K
  X$perc <- 100 * X$delta/sum(X$delta)
  X <- rbind(X, data.frame(variable = c("Total", "Check"),
                           gamma = c(NA, NA), beta_K = c(NA, NA),
                           delta = c(sum(X$delta), coef(reg_base)[var_main] - coef(object)[var_main]),
                           perc= c(NA, NA)))
  
  X
  
}


model_full_1 <- lm(y ~ lag.quarterly.revenue + price.index + income.level + market.potential, data=freeny)
dec_covar_OLD(object = model_full_1, var_main = "lag.quarterly.revenue")
dec_covar(object = model_full_1, var_main = "lag.quarterly.revenue")

smodel_full_2 <- lm(swiss)
dec_covar_OLD(object = model_full_2, var_main = "Catholic")


## plm
library(plm)
data("Produc", package = "plm")
zz <- plm(log(gsp) ~ pcap + log(pc) + log(emp) + unemp,
          data = Produc, index = c("state","year"))

dec_covar(object=zz, var_main = "pcap")

