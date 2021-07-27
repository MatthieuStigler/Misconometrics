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

## Utility function get_response
get_response <- function(x)  
  UseMethod("get_response")

get_response.default <- function(x) {
  deparse(formula(x)[[2]])
}

get_response.felm <- function(x) {
  colnames(x$response)
}

drop_terms <- function(f, dr) {
  f <- Formula::as.Formula(f)
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

## Complicated function to update ivreg, just updating code
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

##### fast auxiliary regression

#' Auxilary regression of covariates on main regressor
#' 
#' @param object regression object
#' @param var_main Main variable
#' @param var_controls Subset of covariates if desired (NULL as default)
#' 
#' @return An extended object of class \code{reg_aux}, on top of the original class
#' @examples 
#' model_full_1 <- lm(y ~ lag.quarterly.revenue + price.index + income.level + market.potential, data=freeny)
#' reg_aux.lm(model_full_1, var_main="lag.quarterly.revenue")
#' reg_aux.lm(model_full_1, var_main="lag.quarterly.revenue", var_controls="price.index")
#' 
#' ## add vcov
#' vcov(reg_aux.lm(model_full_1, var_main="lag.quarterly.revenue", add_vcov=TRUE))


reg_aux <- function(x)  UseMethod("reg_aux")

## not working!! :-()
## need
reg_aux.default <- function(object, var_main, var_controls = NULL) {
  warning("Actually not working...")
  if(is.null(var_controls)) var_controls <- attr(terms(object), "term.labels")
  
  ## exclude main
  if(var_main %in%var_controls) var_controls <- var_controls[-which(var_controls ==var_main)]
  
  ## wrong code!
  string_formula <- sprintf("%s ~ %s", var_controls, paste(var_controls, collapse=" + "))
  print(string_formula)
  reg_aux_all <- lapply(string_formula, function(x) update(object, as.formula(x)))  
}

#' @param method Method to use, update_lmfit is meant to be faster
#' @param add_vcov Compute also variance-covariance matrix?
reg_aux.lm <- function(object, var_main, var_controls = NULL, method = c("update", "update_lmfit",  "sweep"),
                       add_vcov = FALSE) {
  
  method <-  match.arg(method)
  
  if(is.null(var_controls)) {
    var_all <- attr(terms(object), "term.labels")
    var_controls <-  var_all[var_all!= var_main]
  }
  
  if(method == "update") {
    string_formula <- sprintf("cbind(%s) ~ %s", toString(var_controls), paste(var_main, collapse=" + "))
    res <- update(object, as.formula(string_formula)) 
    old_class <- class(res)
    class(res) <-  c("reg_aux_lm", "reg_aux", old_class)
  } else  if(method=="update_lmfit") {
    MM <-  model.matrix(object)
    X <-  MM[, c("(Intercept)", var_main)]
    Y <- MM[, var_controls]
    
    res <- lm.fit(X, Y)
    class(res) <-  c("reg_aux_lm", "reg_aux")
    
    if(add_vcov) {
      rss <- crossprod(res$residuals)
      resvar <- rss/res$df.residual
      p1 <- 1L:res$rank
      R <- chol2inv(res$qr$qr[p1, p1, drop = FALSE])
      # Rinv <- diag(rowSums(backsolve(res$qr$qr, diag(res$rank))^2)) not faster, unlike: 
      # https://stackoverflow.com/questions/39568978/how-to-calculate-variance-of-least-squares-estimator-using-qr-decomposition-in-r
      VC <-  resvar %x% R
      VC_names <-  paste(rep(var_controls, each = length(var_main)+1),
                         rep(c("(Intercept)", var_main), times = length(var_controls)), sep=":")
      colnames(VC) <- rownames(VC) <- VC_names
      res$vcov <-   VC
    }
    
  } else {
  
    require(ISR3)
    XX <-  crossprod(qr.R(object$qr))
    which_main <- which(var_main == colnames(XX))
    var_regs <- c(1, which_main)
    sweep_lm <- ISR3::SWP(XX, var_regs)
    coef <-  sweep_lm[var_regs, - var_regs]
    res <-  list(coefficients = coef)
    if(add_vcov) {
      N <- object$df.residual + object$rank
      df.residual <- N - 2
      S <- sweep_lm[-var_regs, -var_regs]
      VC <-  (-S/df.residual) %x% sweep_lm[var_regs, var_regs]
      VC_names <-  paste(rep(var_controls, each = length(var_main)+1),
                         rep(c("(Intercept)", var_main), times = length(var_controls)), sep=":")
      colnames(VC) <- rownames(VC) <- VC_names
      res$vcov <- VC
      res$df.residual <-  df.residual
    } 
    class(res) <-  c("reg_aux_lm", "reg_aux")
  }
  attr(res, "method") <-  method
  res
}




reg_aux.felm <-  function(object, var_main, var_controls = NULL, method = c("update", "sweep"),
         add_vcov = FALSE) {
  
  method <-  match.arg(method)
  
  if(is.null(var_controls)) {
    var_all <- attr(terms(object), "term.labels")
    var_controls <-  var_all[var_all!= var_main]
  }
  
  if(method == "update") {
    string_formula <- sprintf("%s ~ %s", paste(var_controls, collapse=" + "), paste(var_main, collapse=" + "))  
    res <- update(object, as.formula(string_formula)) 
    old_class <- class(res)
    class(res) <-  c("reg_aux", old_class)
  } else  if(method=="sweep") {
    require(ISR3)
    vc_get_raw <-  function(x, type="iid") vcov(x, type = type) / summary(x)$rse^2
    vc_raw <- vc_get_raw(object)
    which_main <- which(var_main == colnames(vc_raw))
    which_controls <- which(colnames(vc_raw) %in% var_controls)
    
    sweep_lm <- ISR3::RSWP(vc_raw, which_controls)
    coef <-  sweep_lm[which_main, which_controls]
    res <-  list(coefficients = coef)
    if(add_vcov) {
      # N <- object$df.residual + object$rank
      df.residual <- object$df.residual - 1
      S <- sweep_lm[-which_main, -which_main, drop = FALSE]
      VC <-  (-S/df.residual) %x% sweep_lm[which_main, which_main, drop = FALSE]
      VC_names <-  paste(rep(var_controls, each = length(var_main)),
                         rep(var_main, times = length(var_controls)), sep=":")
      colnames(VC) <- rownames(VC) <- VC_names
      res$vcov <- VC
      res$df.residual <-  df.residual
      class(res) <-  c("reg_aux")
    } 
  
  }
  attr(res, "method") <-  method
  res
}

summary.reg_aux_lm <- function(object, ...) {
  
  
  method <-  attr(object, "method")
  if(method %in% c("sweep", "update_lmfit")) {
    se_all <- sqrt(diag(object$vcov))
    se <-  se_all #[-seq(1, by =2, length.out = length(se_all)/2)]
    est <- c(object$coefficients)
    tval <- est/se
    rdf <-  object$df.residual
    object$coefficients <- cbind(Estimate = est, 
                                 `Std. Error` = se, 
                                 `t value` = tval, 
                                 `Pr(>|t|)` = 2 * pt(abs(tval), rdf, lower.tail = FALSE))  
    return(object)
  } else {
    return(stats:::summary.mlm(object))
  }
  
  
}

vcov.reg_aux_lm <- function(object, ...) {
  
  
  method <-  attr(object, "method")
  if(method %in% c("sweep", "update_lmfit")) {
    object$vcov
  } else {
    NextMethod(object)
  }
  
}


## test
if(FALSE) {
  library(ggplot2)

  model_full_1 <- lm(y ~ lag.quarterly.revenue + price.index + income.level + market.potential, data=freeny)
  model_heavy <- lm(carat ~ depth + table + price + x +y+z, data = diamonds)
  
  res_lm <- reg_aux.lm(object=model_full_1, var_main = "lag.quarterly.revenue")
  res_lmf <- reg_aux.lm(object=model_full_1, var_main = "lag.quarterly.revenue", method = "update_lmfit", add_vcov = TRUE)
  res_lmf_heavy <- reg_aux.lm(object=model_heavy, var_main = "depth", method = "update_lmfit", add_vcov = TRUE)
  res_sweep <- reg_aux.lm(model_full_1, var_main = "lag.quarterly.revenue", method = "sweep", add_vcov = TRUE)
  
  res_li <-  list(res_lm= res_lm, res_lmf = res_lmf, 
                  res_sweep = res_sweep)
  
  
  lapply(res_li, coef)
  
  all.equal(coef(res_lmf), coef(res_sweep), coef(res_lm))
  
  
  coef(summary(object = res_lm)) 
  coef(summary(object = res_lmf))
  coef(summary(object = res_sweep))
  all.equal(coef(summary(res_lmf)), coef(summary(res_sweep)))

  
  all.equal(vcov(res_lm), vcov(res_lmf))
  all.equal(vcov(res_lm), vcov(res_sweep))
  
  vcov(res_lm)
  vcov(res_lmf)
  vcov(res_sweep)
  
  confint(res_lm)
  confint(res_sweep)
    
  
  ## felm
  library(lfe)
  data("Produc", package = "plm")
  model_felm <- felm(gsp ~ pcap + pc + emp + unemp |state, data = Produc)
  model_felm_clust <- felm(gsp ~ pcap + pc + emp + unemp |state|0|state, data = Produc)

  res_felm_upd <- reg_aux.felm(object=model_felm, var_main = "pcap")
  res_felm_swp <- reg_aux.felm(object=model_felm, var_main = "pcap", method = "sweep", add_vcov = TRUE)
  
  res_felm_upd_clust <- reg_aux.felm(object=model_felm_clust, var_main = "pcap")
  res_felm_swp_clust <- reg_aux.felm(object=model_felm_clust, var_main = "pcap", method = "sweep", add_vcov = TRUE)
  
  
  map_dfr(list("pc", "emp", "unemp"), ~coef(summary(res_felm_upd, lhs = .)) %>%  as.data.frame)
  coef(summary.reg_aux_lm(res_felm_swp))
  
  
  map_dfr(list("pc", "emp", "unemp"), ~coef(summary(res_felm_upd, lhs = .)) %>%  as.data.frame)
  coef(summary.reg_aux_lm(res_felm_swp))
  
  ## compare
  all.equal(coef(res_lm)[2,], coef(res_sweep))
  all.equal(vcov(res_lm), res_sweep$vcov, check.attributes = FALSE)

  
  coef(summary(object = res_lm))
  coef(summary(object = res_sweep))

  library(broom)
  tidy(res)  %>% 
    filter(term !="(Intercept)") %>% 
    as.matrix()
  
  coef(summary(object = res_sweep))
  
  library(microbenchmark)
  microbenchmark(update  = reg_aux.lm(model_full_1, var_main = "lag.quarterly.revenue") %>%  summary,
                 update_lmf  = reg_aux.lm(model_full_1, var_main = "lag.quarterly.revenue", method = "update_lmfit", add_vcov = TRUE)%>%  summary,
                 sweep  = reg_aux.lm(model_full_1, var_main = "lag.quarterly.revenue", method = "sweep", add_vcov = TRUE)%>%  summary,
                 times = 10)
  microbenchmark(update  = reg_aux.lm(model_heavy, var_main = "depth") %>%  summary,
                 update_lmf  = reg_aux.lm(model_heavy, var_main = "depth", method = "update_lmfit", add_vcov = TRUE)%>%  summary,
                 sweep  = reg_aux.lm(model_heavy, var_main = "depth", method = "sweep", add_vcov = TRUE)%>%  summary,
                 times = 10)
  
}



dec_covar <- function(object, var_main, format = c("wide", "long"),
                      add_coefs = FALSE, conf.int = FALSE, ...) {
  
  format <- match.arg(format)
  if(add_coefs & format == "wide") stop("'add_coefs' only for format = 'long'")
  tidy_quick_arg <-  !conf.int
  
  ## get var names
  formu_obj <- formula(object)
  response_var <- get_response(object)
  var_all <- attr(terms(object), "term.labels")
  var_other <- var_all[! var_all %in% var_main]
  
  ## beta_main: base regression: on main variable(s) only
  if(format == "wide" | add_coefs) {
    formu_base <- drop_terms(f=formu_obj, dr=var_other)
    reg_base <- update(object, as.formula(formu_base), ...)
    # reg_base <- update(object, as.formula(formu_base))
  }
  
  ## gamma: auxiliary regs: covariates on main regs
  if(inherits(object, c("lm", "felm"))) {
    if(inherits(object, "lm")) {
      string_formula <- sprintf("cbind(%s) ~ %s", toString(var_other), paste(var_main, collapse=" + "))
    } else {
      string_formula <- sprintf("%s ~ %s", paste(var_other, collapse=" + "), paste(var_main, collapse=" + "))  
    }
    # reg_aux <- update(object, as.formula(string_formula), ...)  
    reg_aux <- update(object, as.formula(string_formula))  
    gamma_df <- broom::tidy(x=reg_aux, conf.int = conf.int, quick = !conf.int) %>% 
      select(response, term, estimate, contains("conf")) %>% 
      filter(term %in% var_main) %>% 
      rename(variable = term, 
             covariate = response,
             gamma = estimate) %>% 
      rename_with(~str_replace(., "^conf\\.", "gamma_"), contains("conf")) %>% 
      # rename_at(vars(contains("conf")), funs(str_replace(., "^conf\\.", "gamma_"))) %>% 
      select(variable, covariate, contains("gamma"))
    
    
    # gamma_df <- coef(reg_aux)[var_main,, drop=FALSE] %>%
    #   as_tibble() %>%
    #   mutate(variable=var_main) %>%
    #   gather(covariate, gamma, -variable)
    
  } else {
    string_formula2 <- sprintf("%s ~ %s", var_other, paste(var_main, collapse=" + "))
    reg_aux_all <- lapply(string_formula2, function(x) update(object, as.formula(x))) 
    gamma_df <- map2_dfr(reg_aux_all, var_other, ~tidy(.x, conf.int = conf.int, quick = !conf.int) %>%
                           mutate(covariate = .y)) %>%
      filter(term %in% var_main) %>%
      rename(variable = term,
             gamma= estimate) %>%
      rename_with(~str_replace(., "^conf\\.", "gamma_"), contains("conf")) %>% 
      # rename_at(vars(contains("conf")), funs(str_replace(., "^conf\\.", "gamma_"))) %>% 
      select(covariate, variable, gamma, contains("gamma"))
  }
  
  ## Assemble results
  res_df <- tidy(object, conf.int = conf.int) %>% 
    select(term, estimate, contains("conf")) %>% 
    rename(covariate = term,
           beta_K = estimate) %>% 
    # rename_at(vars(contains("conf")), funs(str_replace(., "^conf\\.", "beta_K_"))) %>% 
    rename_with(~str_replace(., "^conf\\.", "beta_K_"), contains("conf")) %>% 
    filter(covariate %in% var_other) %>%
    left_join(gamma_df, by="covariate") %>%
    mutate(delta = gamma * beta_K)
  
  ## add coefs in case:
  if(add_coefs) {
    betas_main <- tibble(model = c("full", "base"),
                             reg =list(full = object,
                                       base = reg_base)) %>%
      mutate(reg = map(reg, tidy, quick = TRUE)) %>%
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
    coef_diffs <- tidy(reg_base) %>%
      filter(term !="(Intercept)") %>%
      left_join(tidy(object), by="term", suffix = c("_base", "_full")) %>%
      mutate(diff = estimate_base-estimate_full) %>%
      select(term, diff)
    
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

plot_dec <- function(x, scales = "free") {
  
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
  
  
  if(n_var>1) pl <- pl + facet_grid(. ~ variable, scales=scales)
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
    geom_vline(xintercept = 0, lty=2) +
    xlab("beta")
  
  if(!legend_size) res <- res +
    scale_size_continuous(guide = "none")
  
  if(add_CI) res <-  res +
    geom_errorbar(aes(ymin = gamma_low, ymax = gamma_high), lty = 2) +
    geom_errorbarh(aes(xmin = beta_K_low, xmax = beta_K_high)) 
  res
}

plot_gam_bet_del <-  function(x, colour = covariate, conf.int = TRUE, scales = "fixed"){
  
  colour_var <-  enquo(colour)
  
  if(!conf.int & !all(c("gamma", "delta") %in% colnames(x))) stop("Make sure called with (format='long')")
  if(conf.int & !"gamma_low" %in% colnames(x)) stop("Make sure called with (format='long', conf.int=TRUE)")
    
  ## reshape data
  x_w <- x %>% 
    select(-starts_with("beta_var")) %>% 
    setNames(str_replace(colnames(.), "beta_K", "beta")) %>% 
    rename(beta_point = beta, gamma_point = gamma,
           delta_point = delta) %>% 
    gather(var_stat, value, contains("beta"), contains("gamma"), delta_point) %>% 
    separate(var_stat, c("parameter", "stat")) %>% 
    spread(stat, value)
  
  ## plot
  x_plot <- x_w %>% 
    mutate(parameter  = factor(parameter, levels = c("beta", "gamma", "delta"))) %>% 
    ggplot(aes(y = covariate, x = point, colour = !!colour_var)) +
    geom_point() +
    facet_grid(.~parameter, scales = scales) +
    geom_vline(xintercept = 0, lty = 2) 
  if(conf.int)  x_plot <- x_plot +
    geom_errorbarh(aes(xmin = low, xmax = high)) +
    theme(legend.position = "none") 
  x_plot
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

# assignInNamespace("process_lm", process_lm_mine, "broom")

#'###############################
#'## Test
################################

if(FALSE){
  # library(devtools)
  # devtools::install_github("tidymodels/broom")
  # devtools::install_github("MatthieuStigler/broom", ref = "dev")
  
  ## lm
  model_full_1 <- lm(y ~ lag.quarterly.revenue + price.index + income.level + market.potential, data=freeny)
  dec_covar(object = model_full_1, var_main = "lag.quarterly.revenue")
  dec_lm1_l <- dec_covar(object = model_full_1, var_main = "lag.quarterly.revenue", format="long", add_coefs=TRUE,
                         conf.int = TRUE)
  dec_lm1_l
  plot_dec(dec_lm1_l)
  plot_gamma_beta(x=dec_lm1_l)
  plot_gam_bet_del(x=dec_lm1_l)
  
  dec_lm1_k2_l <- dec_covar(object = model_full_1, var_main = c("lag.quarterly.revenue", "price.index"),
                            format="long", add_coefs=TRUE, conf.int = TRUE)
  plot_dec(dec_lm1_k2_l)
  plot_gamma_beta(dec_lm1_k2_l)
  plot_gam_bet_del(x=dec_lm1_k2_l, colour = variable)
  
  
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
  dec_lfe_l <- dec_covar(object=model_felm, var_main = "log(pc)", format="long", add_coefs=TRUE, conf.int = TRUE)
  plot_dec(dec_lfe_l)
  plot_gam_bet_del(dec_lfe_l)

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
  df <- tibble(data = list(freeny)) %>%
    mutate(reg = map(data, ~lm(y ~ lag.quarterly.revenue + price.index + income.level + market.potential, data=.)))#,
           # dec = map(reg, ~update(., y ~ . - lag.quarterly.revenue))) 
  
  df %>%
    mutate(reg2 = map2(reg, data,  function(x,y) {x$call$data <- y ; update(x, y ~  lag.quarterly.revenue)})) %>%
    .$reg2
  
  ## get_response
  get_response_alter <-  function(x) attr(as.Formula(formula(x)), "lhs")[[1]]
  
  mod_resps <-  list(model_iv_over = model_iv_over,
                     model_full_1 = model_full_1,
                     model_felm = model_felm)
  lapply(mod_resps, get_response)
  lapply(mod_resps, get_response_alter)
  
  
}
