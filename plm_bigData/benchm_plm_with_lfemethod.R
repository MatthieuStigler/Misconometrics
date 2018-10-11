library(plm)
library(lfe)




# data
T = 50L; 
N = 200000L
N * T
set.seed(123)
data_sim <- data_frame(A = rnorm(N * T),
                       B1 = sample(c(0,1), N * T, replace = TRUE),
                       B2 = rnorm(N * T),
                       individual = rep(1:N, each = T),
                       year = rep(1:T, N))

data_sim

data_small <- filter(data_sim, individual %in% 1:300 & year %in% 1:3)
data_medium_N <- filter(data_sim, individual %in% 1:10000 & year %in% 1:3)
data_medium_T <- filter(data_sim, individual %in% 1:300 & year %in% 1:20)
data_large_N <- filter(data_sim, individual %in% 1:100000 & year %in% 1:3)
data_large_T <- filter(data_sim, individual %in% 1:10000 & year %in% 1:20 )
data_fat_N <- filter(data_sim, year %in% 1:3)
data_fat_N_med_T <- filter(data_sim, year %in% 1:10)


##
as_fake_plm <- function(x) {
  x_plm <- as.matrix(x[, 1:3])
  attr(x_plm, "index") <- attr(x, "index")
  x_plm
}

pd_small <-  pdata.frame(data_small, index = c("individual", "year"))
pd_medium_N <-  pdata.frame(data_medium_N, index = c("individual", "year"))
pd_medium_T <-  pdata.frame(data_medium_T, index = c("individual", "year"))
pd_large_N <-  pdata.frame(data_large_N, index = c("individual", "year"))


pd_small_M <-  as_fake_plm(pd_small)

methods(Within)
plm:::Within.matrix

W_plm <- Within(pd_small_M, effect = "individual")

Within_plm_orig <- function (x, effect, rm.null = TRUE, ...) 
{
  if (is.null(attr(x, "index"))) {
    result <- Within.default(x, effect, ...)
    othervar <- .colSums(abs(x), m = nrow(x), n = ncol(x), 
                         na.rm = TRUE) > 1e-12
    if (rm.null) {
      result <- result[, othervar, drop = FALSE]
      attr(result, "constant") <- character(0)
    }    else {
      result <- result[, drop = FALSE]
      attr(result, "constant") <- colnames(x)[!othervar]
    }
    result
  }  else {
    if (effect %in% c("individual", "time", "group")) 
      result <- x - Between(x, effect)
    if (effect == "twoways") {
      xindex <- attr(x, "index")
      if (is.pbalanced(xindex)) {
        result <- x - Between(x, "individual") 
        - Between(x, "time") 
        + matrix(.colMeans(x, m = nrow(x), n = ncol(x)), nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
      } else {
        time <- index(xindex, "time")
        id <- index(xindex, "individual")
        Dmu <- model.matrix(~time - 1)
        attr(Dmu, "index") <- xindex
        W1 <- Within(x, "individual", rm.null = FALSE)
        WDmu <- Within(Dmu, "individual")
        W2 <- fitted(lm.fit(WDmu, x))
        result <- W1 - W2
      }
    }
  }
  result
}


pd_small_M

attr(pd_small_M, "index")

# Within_plm_out_new <- function (x, effect, rm.null = TRUE, ...) {
#   require(data.table)
#   x2 <-  copy(x)
#   x2[, individual:= attr(x, "index")[,effect]]
#   setDT(x, key= effect)
#   
#   mean_x <- x[, lapply(.SD,mean, na.rm=TRUE), by=individual][,individual:=NULL][x$individual,]
#   x[,individual:=NULL]
#   res <- x - mean_x 
#   setDF(res)
#   return(res)
# }

Within_plm_new <- function (x, effect, rm.null = TRUE, ...) {
  require(lfe)
  l <-  as.list(attr(x, "index")[,effect, drop=FALSE])
  demeanlist(x, fl = l)
}

new <-  Within_plm_new(x = pd_small_M, effect = "individual")
old <-  Within_plm_orig(x = pd_small_M, effect = "individual")

names(attributes(new))
names(attributes(old))
all.equal(new, old, check.attributes=FALSE)

plm:::Within.matrix

assignInNamespace("Within.matrix", Within_plm_new, "plm")
plm_out_new <- plm(A~B1+B2, data=pd_small)

assignInNamespace("Within.matrix", Within_plm_orig, "plm")
plm_out_old <- plm(A~B1+B2, data=pd_small)

all.equal(plm_out_new, plm_out_old)


#### 
library(microbenchmark)
microbenchmark(plm_new = {assignInNamespace("Within.matrix", Within_plm_new, "plm"); 
  plm(A~B1+B2, data=pd_medium_N)}, 
  plm_old = {assignInNamespace("Within.matrix", Within_plm_orig, "plm"); 
    plm(A~B1+B2, data=pd_medium_N)},
  times =5
  )


microbenchmark(plm_new = {assignInNamespace("Within.matrix", Within_plm_new, "plm"); 
  plm(A~B1+B2, data=pd_medium_T)}, 
  plm_old = {assignInNamespace("Within.matrix", Within_plm_orig, "plm"); 
    plm(A~B1+B2, data=pd_medium_T)},
  times =5
)

data_large_N2 <- filter(data_sim, individual %in% 1:50000 & year %in% 1:3)
pd_large_N2 <-  pdata.frame(data_large_N2, index = c("individual", "year"))

microbenchmark(plm_new = {assignInNamespace("Within.matrix", Within_plm_new, "plm"); 
  plm(A~B1+B2, data=pd_large_N2)}, 
  plm_old = {assignInNamespace("Within.matrix", Within_plm_orig, "plm"); 
    plm(A~B1+B2, data=pd_large_N2)},
  times =2
)

library(profvis)
profvis({assignInNamespace("Within.matrix", Within_plm_new, "plm"); 
  plm(A~B1+B2, data=pd_large_N2)}
)
