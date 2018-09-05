#' ---
#' title: "Implement Gelbach dec_covarsition"
#' author: "Matthieu Stigler"
#' date: "2018-09-03"
#' output:
#'    pdf_document:
#'       toc: true
#' ---

library(tidyverse)
library(broom)
library(Formula)


get_response <- function(x)  
  UseMethod("get_response")

get_response.default <- function(x) {
  deparse(formula(x)[[2]])
}

get_response.felm <- function(x) {
  colnames(x$response)
}

drop_terms <- function(f, dr) {
  f <- as.Formula(f)
  n_rhs <- length(f)[2]
  if(n_rhs==1) {
    form_rem <- paste(". ~ . -", paste(dr, collapse = " - ")) 
  } else {
    terms_by_part <- lapply(1:length(f)[2], function(x) attr(terms(f, rhs=x), "term.labels"))
    has_terms <- sapply(terms_by_part, function(x) any(dr %in% x))
    inner <- sapply(has_terms, function(x) if(x) paste(" . -", paste(dr, collapse = " - ")) else ".")
    form_rem <- paste(". ~ ", paste(inner, collapse = "|")) 
  }
  
  form_rem
}

update.ivreg <- function (object, formula., ..., evaluate = TRUE) 
{
  if (is.null(call <- getCall(object))) 
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.)) 
    call$formula <- update(Formula(formula(object)), formula.)
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (evaluate) 
    eval(call, parent.frame())
  else call
}


dec_covar <- function(object, var_main, format = c("wide", "long")) {
  
  format <- match.arg(format)
  
  
  ## get var names
  formu_obj <- formula(object)
  response_var <- get_response(object)
  var_all <- attr(terms(object), "term.labels")
  var_other <- var_all[! var_all %in% var_main]
  
  ## base regression: on main variable(s) only
  formu_base <- drop_terms(f=formu_obj, dr=var_other)
  reg_base <- update(object, as.formula(formu_base))
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
  
  # wide version in case 
  if(format == "wide") {
    res_df_w <- res_df %>%
      gather(stat, value, gamma, delta) %>%
      mutate(stat = paste(stat, variable, sep="_")) %>%
      select(-variable) %>%
      spread(stat, value) %>%
      select(covariate, beta_K, starts_with("gamma"), starts_with("delta"))
    
    ## add totals and check
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

#'###############################
#'## Test
################################

if(FALSE){
  
  ## lm
  model_full_1 <- lm(y ~ lag.quarterly.revenue + price.index + income.level + market.potential, data=freeny)
  dec_covar(object = model_full_1, var_main = "lag.quarterly.revenue")
  
  model_full_2 <- lm(Fertility ~ . , data=swiss)
  dec_covar(object = model_full_2, var_main = "Catholic")

  
  ## panel
  library(plm)
  data("Produc", package = "plm")
  zz <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp,
            data = Produc, index = c("state","year"))
  
  dec_covar(object=zz, var_main = "log(pc)")
  
  library(lfe)
  model_felm <- felm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp |state|0|state, data = Produc)
  dec_covar(object=model_felm, var_main = "log(pc)")
  
  ## iv
  library(AER)
  example(ivreg, echo=FALSE)
  model_iv_over <- ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff + I(tax/cpi),
                       data = CigarettesSW, subset = year == "1995")
  dec_covar(object=model_iv_over, var_main = "log(rprice)")
  
  model_iv_just <- ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff ,
                         data = CigarettesSW, subset = year == "1995")
  dec_covar(object=model_iv_just, var_main = "log(rprice)")
  
  ## get_response
  get_response(model_full_1)
  get_response(model_iv_over)
  get_response(model_felm)
  
}
