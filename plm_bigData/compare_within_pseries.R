library(tidyverse)
library(plm)


## data
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

df_all <- data_frame(data_type = c("data_small",  "data_medium_N",  "data_medium_T", "data_large_N", "data_large_T"),
                     data = list(data_small, data_medium_N, data_medium_T, data_large_N, data_large_T)) %>%
  mutate(data_plm =   map(data, ~pdata.frame(., index = c("individual", "year"))),
         n_row = map_int(data, ~nrow(.) %>% as.integer)) %>%
  select(data_type, n_row, everything())

pd_small <-  pdata.frame(data_small, index = c("individual", "year"))
pd_medium_N <-  pdata.frame(data_medium_N, index = c("individual", "year"))
pd_medium_T <-  pdata.frame(data_medium_T, index = c("individual", "year"))
pd_large_N <-  pdata.frame(data_large_N, index = c("individual", "year"))


# Then extract a series, which becomes additionally a pseries
ps_small <-  pd_small$B1
attributes(ps_small) %>%  names



# compute the between and within transformations
Within_matrix <-  function(x, effect) {
  input <- unclass(x)
  dim(input) <- c(length(x), 1)
  res <- plm:::Within.matrix(input, effect= effect) 
  res <-  as.vector(res)
  attributes(res) <-  attributes(x)
  res
}

W_ps <- Within(df_all$data_plm[[1]]$B1)
W_mat <- Within_matrix(x=df_all$data_plm[[1]]$B1, effect="individual")

all.equal(W_ps, W_mat)



#### BENCH
microbenchmark( pseries =  Within(ps_small),
                matrix = Within_matrix(x=ps_small, effect="individual"),
                times = 2)



## within benchmark
benches_within <- df_all %>%
  mutate(n_rep = c(30, 20, 30, 5, 5)) %>% 
  # head(3) %>% 
  mutate(ben = map2(data_plm, n_rep,
                    ~microbenchmark(pseries =  Within(.x$B1),
                                    matrix = Within_matrix(x=.x$B1, effect="individual"),
                                    times=.y)))

benches_within_df <- benches_within %>% 
  mutate(ben = map(ben, ~summary(., unit="ms") %>%
                     mutate(expr = as.character(expr)))) %>%
  unnest(ben)  %>% 
  group_by(data_type) %>% 
  mutate(mean_rel = mean/min(mean)) %>% 
  ungroup()
  
benches_within_df


benches_within_df %>% 
  filter(mean <10000) %>%
  mutate(data_type = factor(data_type, levels = unique(data_type))) %>% 
  ggplot(aes(x = data_type, y = mean_rel, colour = expr)) +
  geom_point(size =3, alpha = I(0.7)) +
  theme(axis.text.x = element_text(angle = 20))

write_csv(benches_within_df,
          "/home/matifou/git/my_github/Misconometrics/plm_bigData/output_benchmark/with_pseries_bench.csv")
write_rds(benches_within,
          "/home/matifou/git/my_github/Misconometrics/plm_bigData/output_benchmark/with_pseries_bench.rda")



### switch methods
plm:::Within.pseries

Within_M_orig <- function (x, effect, rm.null = TRUE, ...) 
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

Within_pseries_orig <- function (x, effect = c("individual", "time", "group", "twoways"), ...) {
  effect <- match.arg(effect)
  if (effect != "twoways") 
    x - Between(x, effect, ...)
  else {
    if (is.pbalanced(x)) 
      x - Between(x, "individual", ...) - Between(x, "time") + 
      mean(x)
    else {
      time <- index(x)[[2]]
      Dmu <- model.matrix(~time - 1)
      attr(Dmu, "index") <- index(x)
      W1 <- Within(x, "individual")
      WDmu <- Within(Dmu, "individual")
      W2 <- fitted(lm.fit(WDmu, x))
      W1 - W2
    }
  }
}

Within_M_new <- function (x, effect, rm.null = TRUE, ...) {
  require(lfe)
  l <-  as.list(attr(x, "index")[,effect, drop=FALSE])
  demeanlist(x, fl = l)
}

W_M_new <-  function() assignInNamespace("Within.matrix", Within_M_new, "plm")
W_M_old <-  function() assignInNamespace("Within.pseries", Within_M_orig, "plm")
W_ps_new <-  function() assignInNamespace("Within.pseries", Within_matrix, "plm")
W_ps_old <-  function() assignInNamespace("Within.pseries", Within_pseries_orig, "plm")

## test


plm_orig <- plm(A~B1+B2, data=pd_small)


plm_new <- plm(A~B1+B2, data=pd_small)

all.equal(plm_new, plm_orig)


#### 
library(microbenchmark)

microbenchmark(plm_oldM_newPS = {W_M_old(); W_ps_new(); plm(A~B1+B2, data=pd_medium_N)}, 
               plm_oldM_oldPS = {W_M_old(); W_ps_old(); plm(A~B1+B2, data=pd_medium_N)},
               plm_newM_oldPS = {W_M_new(); W_ps_old(); plm(A~B1+B2, data=pd_medium_N)},
               plm_newM_newPS = {W_M_new(); W_ps_new(); plm(A~B1+B2, data=pd_medium_N)},
  times =5)



W_M_old();

benches_plm_psnew <- df_all %>%
  mutate(n_rep = c(30, 20, 30, 5, 5)) %>% 
  mutate(ben = map2(data_plm, n_rep,
                    ~microbenchmark(plm_oldM_newPS = {W_ps_new(); plm(A~B1+B2, data=.x)}, 
                                    plm_oldM_oldPS = {W_ps_old(); plm(A~B1+B2, data=.x)},
                                    times=.y))) %>% 
  mutate(ben = map(ben, ~summary(., unit="ms") %>%
                     mutate(expr = as.character(expr)))) %>%
  unnest(ben)  %>% 
  group_by(data_type) %>% 
  mutate(mean_rel = mean/min(mean)) %>% 
  ungroup()

benches_plm_psnew %>% 
  mutate(mean_sec  = 0.001 * mean)


pl_plm_nwpseries <- benches_plm_psnew %>% 
  filter(mean <10000) %>%
  mutate(data_type = factor(data_type, levels = unique(data_type))) %>% 
  ggplot(aes(x = data_type, y = mean, colour = expr)) +
  geom_point(size =3, alpha = I(0.7)) +
  theme(axis.text.x = element_text(angle = 20)) +
  ggtitle("Running plm with Within.pseries (blue) or Within.matrix") +
  ylab("milliseconds")

pl_plm_nwpseries


ggsave(pl_plm_nwpseries,
       filename = "/home/matifou/git/my_github/Misconometrics/plm_bigData/figures/bench_plm_without_pseries.png", width = 7, height = 4)
