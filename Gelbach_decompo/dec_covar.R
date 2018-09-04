#' ---
#' title: "Implement Gelbach dec_covarsition"
#' author: "Matthieu Stigler"
#' date: "2018-09-03"
#' output:
#'    pdf_document:
#'       toc: true
#' ---



get_response <- function(x)  
  UseMethod("get_response")

get_response.default <- function(x) {
  deparse(formula(x)[[2]])
}

get_response.felm <- function(x) {
  colnames(x$response)
}

dec_covar <- function(object, var_main) {
  
  
  ## get var names
  formu_obj <- formula(object)
  response_var <- get_response(object)
  var_all <- attr(terms(object), "term.labels")
  var_other <- var_all[var_all!=var_main]
  
  ## reg base
  formu_base <- paste(response_var, "~", var_main)
  reg_base <- update(object, as.formula(formu_base))
  # reg_base <- update(object, reformulate(var_main, response_var ))
  
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
  
  ## Assemble results
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
  model_felm <- felm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp |state, data = Produc)
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
  
  coef(model_iv_just)["log(rprice)"]
  model_iv_just_base <- ivreg(log(packs) ~ log(rprice)  | tdiff ,
                              data = CigarettesSW, subset = year == "1995")
  update(model_iv_just, log(packs) ~ log(rprice) |  tdiff )
  coef(model_iv_just)["log(rprice)"] -
  coef(model_iv_just_base)["log(rprice)"]
  
  ## get_response
  get_response(model_full_1)
  get_response(model_iv_over)
  get_response(model_felm)
  
}
