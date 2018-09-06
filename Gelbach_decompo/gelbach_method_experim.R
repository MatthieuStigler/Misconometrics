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


################################
#'## Not all covariates?
################################


reg_full <- lm(y ~ lag.quarterly.revenue + price.index + income.level + market.potential, data=freeny)
reg_base_2 <- lm(y ~ lag.quarterly.revenue + price.index, data=freeny)

freeny_M <- as.matrix(freeny)
reg_aux_2 <- lm(freeny_M[, -c(1:3)]~freeny_M[, 2:3])

## diff
coef(reg_base_2)[2] - coef(reg_full)[2]
coef(reg_base_2)[3] - coef(reg_full)[3]

## dec_covar:
coef(reg_aux_2)[2,] %*% coef(reg_full)[-c(1,2, 3)]
coef(reg_aux_2)[3,] %*% coef(reg_full)[-c(1,2, 3)]


## in df
X <- data_frame(variable  = names(coef(reg_full))[-c(1,2,3)],
                gamma_1 = coef(reg_aux_2)[2,],
                gamma_2 = coef(reg_aux_2)[3,],
                beta_K = coef(reg_full)[-c(1,2, 3)],
                delta_1 = gamma_1 * beta_K,
                delta_2 = gamma_2 * beta_K) %>%
  bind_rows(data_frame(variable = "Total", 
                       delta_1 = sum(.$delta_1),
                       delta_2 = sum(.$delta_2))) %>%
  bind_rows(data_frame(variable = "Check", 
                       delta_1 = coef(reg_base_2)[2] - coef(reg_full)[2],
                       delta_2 = coef(reg_base_2)[3] - coef(reg_full)[3]))
X


drop_terms <- function(f, dr) {
  f <- as.Formula(f)
  n_rhs <- length(f)[2]
  if(n_rhs==1) {
    form_rem <- paste(". ~ . -", paste(dr, collapse = " - ")) 
  } else {
    terms_by_part <- sapply(1:length(f)[2], function(x) attr(terms(f, rhs=x), "term.labels"))
    has_terms <- sapply(terms_by_part, function(x) any(dr %in% x))
    inner <- sapply(has_terms, function(x) if(x) paste(" . -", paste(dr, collapse = " - ")) else ".")
    form_rem <- paste(". ~ ", paste(inner, collapse = "|")) 
  }
  
  form_rem
}

dec_covar <- function(object, var_main, format = c("wide", "long"),
                      add_coefs = FALSE) {
  
  format <- match.arg(format)
  if(add_coefs & format == "wide") stop("'add_coefs' only for format = 'long'")
  
  ## get var names
  formu_obj <- formula(object)
  response_var <- get_response(object)
  var_all <- attr(terms(object), "term.labels")
  var_other <- var_all[! var_all %in% var_main]
  
  ## base regression: on main variable(s) only
  formu_base <- drop_terms(formu_obj, var_other)
  reg_base <- update(object, as.Formula(formu_base))
  coef_diffs <- tidy(reg_base) %>%
    filter(term !="(Intercept)") %>%
    left_join(tidy(object), by="term", suffix = c("_base", "_full")) %>%
    mutate(diff = estimate_base-estimate_full) %>%
    select(term, diff)
  
  
  ## auxiliary regs: covariates on main regs
  if(inherits(object, "lm")) {
    string_formula <- sprintf("cbind(%s) ~ %s", toString(var_other), paste(var_main, collapse=" + "))
    reg_aux <- update(object, as.formula(string_formula))  
    gamma_df <- coef(reg_aux)[var_main,, drop=FALSE] %>%
      as_data_frame() %>%
      mutate(variable=var_main) %>%
      gather(covariate, gamma, -variable)
  } else {
    string_formula2 <- sprintf("%s ~ %s", var_other, paste(var_main, collapse=" + "))
    reg_aux_all <- lapply(string_formula2, function(x) update(object, as.formula(x))) 
    gamma_df <- map2_dfr(reg_aux_all, var_other, ~tidy(.x) %>% mutate(covariate = .y)) %>%
      filter(term %in% var_main) %>%
      rename(variable = term,
             gamma= estimate) %>%
      select(covariate, variable, gamma)
  }
  
  ## Assemble results
  res_df <- data_frame(covariate  = var_other,
                       beta_K = coef(object)[var_other]) %>%
    left_join(gamma_df, by="covariate") %>%
    mutate(delta = gamma * beta_K)
  
  ## add coefs in case:
  if(add_coefs) {
    betas_main <- data_frame(model = c("full", "base"),
                      reg =list(full = object,
                                base = reg_base)) %>%
      mutate(reg = map(reg, tidy)) %>%
      unnest(reg) %>%
      select(model, term, estimate) %>%
      filter(term %in% var_main) %>%
      mutate(model = paste("beta_var", model, sep="_")) %>%
      spread(model, estimate) %>%
      rename(variable = term)
    
    res_df <- res_df %>%
      left_join(betas_main, by = "variable") 
  }
  
  # wide version in case 
  if(format == "wide") {
    res_df_w <- res_df %>%
      gather(stat, value, gamma, delta) %>%
      mutate(stat = paste(stat, variable, sep="_")) %>%
      select(-variable) %>%
      spread(stat, value) %>%
      select(covariate, beta_K, starts_with("gamma"), starts_with("delta"))
    
    res_df_w_tot <- res_df_w %>%
      bind_rows(summarise_at(res_df_w, vars(starts_with("delta")), sum) %>% 
                  mutate(covariate="Total")) %>%
      bind_rows(coef_diffs %>% mutate(covariate="Check",
                                      term = paste("delta_", term, sep="")) %>% spread(term, diff))
    res_df <- res_df_w_tot
  }
  
  ## Export results
  res_df
  
}





f3 <- drop_terms(f= y ~lag.quarterly.revenue + price.index + income.level + market.potential , 
           dr = c("income.level", "market.potential"))


f4 <- drop_terms(f= y ~lag.quarterly.revenue + price.index + income.level + market.potential| price.index, 
           dr = c("income.level", "market.potential"))

update(model_full_1, f3)

model_full_1 <- lm(y ~ lag.quarterly.revenue + price.index + income.level + market.potential, data=freeny)
dec_covar2(object = model_full_1, var_main = c("lag.quarterly.revenue", "price.index"))
dec_covar2(object = model_full_1, var_main = c("lag.quarterly.revenue"))

################################
#'## Plot
################################


plot_dec <- function(x) {
  n_var <- length(unique(x$variable))
  pl <- x %>%
    mutate(delta_center = beta_var_base - delta) %>%
    ggplot(aes(x = delta_center, y =covariate)) +
    geom_point() +
    geom_segment(aes(x=beta_var_base, xend = delta_center, yend = covariate)) +
    geom_vline(aes(xintercept = c(beta_var_base)), lty=2, colour="blue")+
    geom_vline(aes(xintercept = c(beta_var_full)), lty=2) 
  if(n_var>1) pl <- pl + facet_grid(. ~ variable, scales="free")
  pl
}


dec_1 <- dec_covar(object = model_full_1, var_main = c("lag.quarterly.revenue"), format="long",
                   add_coefs = TRUE)

dec_2 <- dec_covar(object = model_full_1, var_main = c("lag.quarterly.revenue", "price.index"), format="long",
                   add_coefs = TRUE)
plot_dec(x=dec_1)
plot_dec(x=dec_2)

