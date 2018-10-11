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


dec_covar <- function(object, var_main, format = c("wide", "long"),
                      add_coefs = FALSE, conf.int = FALSE) {
  
  format <- match.arg(format)
  if(add_coefs & format == "wide") stop("'add_coefs' only for format = 'long'")
  
  ## get var names
  formu_obj <- formula(object)
  response_var <- get_response(object)
  var_all <- attr(terms(object), "term.labels")
  var_other <- var_all[! var_all %in% var_main]
  
  ## beta_main: base regression: on main variable(s) only
  formu_base <- drop_terms(f=formu_obj, dr=var_other)
  reg_base <- update(object, as.formula(formu_base))
  coef_diffs <- tidy(reg_base) %>%
    filter(term !="(Intercept)") %>%
    left_join(tidy(object), by="term", suffix = c("_base", "_full")) %>%
    mutate(diff = estimate_base-estimate_full) %>%
    select(term, diff)
  
  
  ## gamma: auxiliary regs: covariates on main regs
  if(inherits(object, c("lm", "felm"))) {
    if(inherits(object, "lm")) {
      string_formula <- sprintf("cbind(%s) ~ %s", toString(var_other), paste(var_main, collapse=" + "))
    } else {
      string_formula <- sprintf("%s ~ %s", paste(var_other, collapse=" + "), paste(var_main, collapse=" + "))  
    }
    reg_aux <- update(object, as.formula(string_formula))  
    gamma_df <- tidy(x=reg_aux, conf.int = conf.int) %>% 
      select(response, term, estimate, contains("conf")) %>% 
      filter(term %in% var_main) %>% 
      rename(variable = term, 
             covariate = response,
             gamma = estimate) %>% 
      rename_at(vars(contains("conf")), funs(str_replace(., "^conf\\.", "gamma_"))) %>% 
      select(variable, covariate, contains("gamma"))
    
    
    # gamma_df <- coef(reg_aux)[var_main,, drop=FALSE] %>%
    #   as_data_frame() %>%
    #   mutate(variable=var_main) %>%
    #   gather(covariate, gamma, -variable)
    
  } else {
    string_formula2 <- sprintf("%s ~ %s", var_other, paste(var_main, collapse=" + "))
    reg_aux_all <- lapply(string_formula2, function(x) update(object, as.formula(x))) 
    gamma_df <- map2_dfr(reg_aux_all, var_other, ~tidy(.x, conf.int = conf.int) %>%
                           mutate(covariate = .y)) %>%
      filter(term %in% var_main) %>%
      rename(variable = term,
             gamma= estimate) %>%
      rename_at(vars(contains("conf")), funs(str_replace(., "^conf\\.", "gamma_"))) %>% 
      select(covariate, variable, gamma, contains("gamma"))
  }
  
  ## Assemble results
  res_df <- tidy(object, conf.int = conf.int) %>% 
    select(term, estimate, contains("conf")) %>% 
    rename(covariate = term,
           beta_K = estimate) %>% 
    rename_at(vars(contains("conf")), funs(str_replace(., "^conf\\.", "beta_K_"))) %>% 
    filter(covariate %in% var_other) %>%
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

plot_dec <- function(x) {
  
  ## get df of main beta
  beta_orig <- x %>% 
    select(variable, starts_with("beta_var")) %>% 
    distinct() %>% 
    gather(Model, value, starts_with("beta_var")) %>% 
    mutate(Model = case_when(Model == "beta_var_base" ~ "model base",
                             Model == "beta_var_full" ~ "model full"))
  
  n_var <- length(unique(x$variable))
  pl <- x %>%
    mutate(delta_center = beta_var_base - delta) %>%
    ggplot(aes(x = delta_center, y =covariate)) +
    geom_point() +
    geom_segment(aes(x=beta_var_base, xend = delta_center, yend = covariate)) +
    geom_vline(aes(xintercept = value, colour = Model, linetype = Model), data = beta_orig) +
    scale_colour_manual(values = c("model base"= "black", "model full" = "blue"))+
    theme(legend.position = "bottom") +
    xlab("Delta ( = gamma * beta)")
  
  
  if(n_var>1) pl <- pl + facet_grid(. ~ variable, scales="free")
  pl
}

plot_gamma_beta <-  function(x, colour = covariate, size = abs(delta), legend_size = FALSE,
                             add_CI = FALSE) {
  
  colour_var <-  enquo(colour)
  size_var <-  enquo(size)
  
  ## make sure data is in long?
  if(!"gamma" %in% colnames(x)) stop("Data not in 'long' format?")
  
  res <- x %>% 
    ggplot(aes(x = beta_K, y = gamma, colour = !!colour_var))+
    geom_point(aes(size = !!size_var)) +
    geom_hline(yintercept = 0, lty=2) +
    geom_vline(xintercept = 0, lty=2) 
  
  if(!legend_size) res <- res +
    scale_size_continuous(guide = "none")
  
  if(add_CI) res <-  res +
    geom_errorbar(aes(ymin = gamma_low, ymax = gamma_high), lty = 2) +
    geom_errorbarh(aes(xmin = beta_K_low, xmax = beta_K_high)) 
  res
}


## should be removed in near future, see: https://github.com/tidymodels/broom/issues/510
process_lm_mine <- function(ret, x, conf.int = FALSE, conf.level = .95,
                       exponentiate = FALSE) {
  if (exponentiate) {
    # save transformation function for use on confidence interval
    if (is.null(x$family) ||
        (x$family$link != "logit" && x$family$link != "log")) {
      warning(paste(
        "Exponentiating coefficients, but model did not use",
        "a log or logit link function."
      ))
    }
    trans <- exp
  } else {
    trans <- identity
  }
  
  if (conf.int) {
    # avoid "Waiting for profiling to be done..." message
    if(inherits(x, "mlm")) {
      .format.perc <- function (probs, digits) {
        paste(
          format(
            100 * probs,
            trim = TRUE,
            scientific = FALSE,
            digits = digits
          ),
          "%"
        )
      }
      confint.mlm <- function (object, level = 0.95, ...) {
        cf <- coef(object)
        ncfs <- as.numeric(cf)
        a <- (1 - level)/2
        a <- c(a, 1 - a)
        fac <- qt(a, object$df.residual)
        pct <- .format.perc(a, 3)
        ses <- sqrt(diag(stats::vcov(object)))
        ci <- ncfs + ses %o% fac
        setNames(data.frame(ci), pct)
      }
      CI <- confint.mlm(x, level = conf.level)
    } else {
      CI <- suppressMessages(stats::confint(x, level = conf.level))
    }
    # Handle case if regression is rank deficient
    p <- x$rank
    if (!is.null(p) && !is.null(x$qr) & !inherits(x, "mlm")) {
      piv <- x$qr$pivot[seq_len(p)]
      CI <- CI[piv, , drop = FALSE]
    }
    colnames(CI) <- c("conf.low", "conf.high")
    ret <- cbind(ret, trans(broom:::unrowname(CI)))
  }
  ret$estimate <- trans(ret$estimate)
  
  as_tibble(ret)
}

assignInNamespace("process_lm", process_lm_mine, "broom")

#'###############################
#'## Test
################################

if(FALSE){
  library(devtools)
  devtools::install_github("tidymodels/broom")
  
  ## lm
  model_full_1 <- lm(y ~ lag.quarterly.revenue + price.index + income.level + market.potential, data=freeny)
  dec_covar(object = model_full_1, var_main = "lag.quarterly.revenue")
  dec_lm1_l <- dec_covar(object = model_full_1, var_main = "lag.quarterly.revenue", format="long", add_coefs=TRUE,
                         conf.int = TRUE)
  dec_lm1_l
  plot_dec(dec_lm1_l)
  plot_gamma_beta(x=dec_lm1_l)
  
  dec_lm1_k2_l <- dec_covar(object = model_full_1, var_main = c("lag.quarterly.revenue", "price.index"),
                            format="long", add_coefs=TRUE)
  plot_dec(dec_lm1_k2_l)
  plot_gamma_beta(dec_lm1_k2_l)
  
  
  model_full_2 <- lm(Fertility ~ . , data=swiss)
  dec_covar(object = model_full_2, var_main = "Catholic")
  dec_lm2_l <- dec_covar(object = model_full_2, var_main = "Catholic", format="long", add_coefs=TRUE)
  plot_dec(dec_lm2_l)
  plot_gamma_beta(dec_lm2_l)

  ## panel
  library(plm)
  data("Produc", package = "plm")
  zz <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp,
            data = Produc, index = c("state","year"))
  
  dec_covar(object=zz, var_main = "log(pc)")
  dec_plm_l <- dec_covar(object=zz, var_main = "log(pc)", format="long", add_coefs=TRUE)
  plot_dec(dec_plm_l)
  plot_gamma_beta(dec_plm_l)
  
  library(lfe)
  model_felm <- felm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp |state|0|state, data = Produc)
  dec_covar(object=model_felm, var_main = "log(pc)")
  dec_lfe_l <- dec_covar(object=model_felm, var_main = "log(pc)", format="long", add_coefs=TRUE)
  plot_dec(dec_lfe_l)

  dec_lfe_k2_l <- dec_covar(object=model_felm, var_main = c("log(pc)","log(pcap)"), format="long", add_coefs=TRUE)
  plot_dec(dec_lfe_k2_l)
  plot_gamma_beta(dec_lfe_k2_l)

  ## iv
  library(AER)
  example(ivreg, echo=FALSE)
  model_iv_over <- ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff + I(tax/cpi),
                       data = CigarettesSW, subset = year == "1995")
  dec_covar(object=model_iv_over, var_main = "log(rprice)")
  
  model_iv_just <- ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff ,
                         data = CigarettesSW, subset = year == "1995")
  dec_iv <- dec_covar(object=model_iv_just, var_main = "log(rprice)")
  plot_dec(dec_iv)
  
  ## in map calls
  df <- data_frame(data = list(freeny)) %>%
    mutate(reg = map(data, ~lm(y ~ lag.quarterly.revenue + price.index + income.level + market.potential, data=.)))#,
           # dec = map(reg, ~update(., y ~ . - lag.quarterly.revenue))) 
  
  df %>%
    mutate(reg2 = map2(reg, data,  function(x,y) {x$call$data <- y ; update(x, y ~  lag.quarterly.revenue)})) %>%
    .$reg2
  
  ## get_response
  get_response(model_full_1)
  get_response(model_iv_over)
  get_response(model_felm)
  
}
