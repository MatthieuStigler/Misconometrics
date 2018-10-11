library(data.table)
library(plm)
library(broom)
#library(doBy)
library(lfe)
library(microbenchmark)
library(tidyverse)

## tidyverse
Within1_dplyr <- function(x, vari) {
  vari <- enquo(vari)
  x %>%
    group_by(!!vari) %>%
    mutate_all(funs(. - mean(.))) %>%
    ungroup() %>%
    select(-!!vari)
}

## data table
Within1_dt_1 <- function(x, fe1, x_names = NULL) {

  dt <- as.data.table(x)
  if(is.null(x_names)) {
    x_names <- colnames(dt)
    x_names <- x_names[x_names != fe1]
  }
  
  keep_cols <- c(x_names, fe1)
  model_dt <- dt[, keep_cols, with = FALSE]
  model_dt <- na.omit(model_dt)
  
  colnames_new <- colnames(model_dt)
  colnames_new[colnames_new ==fe1] <- "fe1"
  
  setnames(model_dt, 1:ncol(model_dt), colnames_new)
  
  # Means by FE1:
  setkey(model_dt, fe1)
  res <- model_dt[, lapply(.SD, function(x) x -mean(x)), by= fe1]
  as.data.frame(res[, -"fe1", with=FALSE])
  
}


Within1_dt_2 <- function(x, fe1, x_names = NULL) {
  dt <- as.data.table(x)
  if(is.null(x_names)) {
    x_names <- colnames(dt)
    x_names <- x_names[x_names != fe1]
  }
  x_names_new <-  paste(x_names, "mean", sep="_")
  
  setkeyv(dt, fe1)
  dt[,  paste(x_names, "mean", sep="_"):= lapply(.SD, mean, na.rm = TRUE), by = fe1]
  res <- dt[, x_names, with=FALSE] - dt[, x_names_new, with=FALSE]
  as.data.frame(res)
}

Within1_dt_3a <- function(x) {
  setDT(x,key="individual")
  x_fullnames <- names(x)
  mean_x <- x[,lapply(.SD,mean, na.rm=TRUE),by=individual]
  setcolorder(mean_x,x_fullnames)
  setkey(mean_x,individual)
  res<-x-mean_x[x,][,1:length(names(mean_x))]
  res[,individual:=NULL]
  setDF(res)
  return(res)
}

Within1_dt_3b <- function(t) {
  setDT(t,key="individual")
  x_fullnames <- names(t)
  x_names <- x_fullnames[which(x_fullnames != "individual")]
  mean_x <- t[,lapply(.SD,mean,na.rm=TRUE),by=individual]
  res<-t[mean_x,.(A=x.A-i.A,B1=x.B1-i.B1,B2=x.B2-i.B2),on="individual"]
  setDF(res)
  return(res)
}

Within1_dt_3c <- function(x) {
  setDT(x,key="individual")
  mean_x <- x[, lapply(.SD,mean, na.rm=TRUE), by=individual][,individual:=NULL][x$individual,]
  x[,individual:=NULL]
  res <- x - mean_x 
  setDF(res)
  return(res)
}


Within_doby <- function(x, fe1) {
  formu <- as.formula(paste(". ~ ", fe1))
  scaleBy(formu, data=x, scale = FALSE)
}


as_fake_plm <- function(x, keep = 1:3) {
  x_plm <- as.matrix(x[, keep, drop=FALSE])
  attr(x_plm, "index") <- attr(x, "index")
  x_plm
}



##################################################
########## Generate data, assemble in df of df
##################################################
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



## data frame of data frames

df_all <- data_frame(data_type = c("data_small",  "data_medium_N",  "data_medium_T", "data_large_N", "data_large_T", "data_fat_N", "data_fat_N_med_T"),
                     data = list(data_small, data_medium_N, data_medium_T, data_large_N, data_large_T, data_fat_N, data_fat_N_med_T)) %>%
  mutate(data_plm =   map(data, ~pdata.frame(., index = c("individual", "year"))),
         data_plm_M = map(data_plm, as_fake_plm), 
         n_row = map_int(data, ~nrow(.) %>% as.integer)) %>%
  select(data_type, n_row, everything())

plm(A~ B1 , data=df_all$data_plm[[1]], model="random")

##################################################
########## Check results same
##################################################

head_val <- 3 ## if interactive testing or large R CMD BATCH
times_bench <- 10

## test regressions
regs_test <- df_all %>%
  head(head_val) %>%
  mutate(reg_plm = map(data_plm, ~plm(A~B1 + B2, data=.)),
         reg_felm = map(data, ~felm(A~B1 + B2|individual, data=.)))

regs_test

## check equality of coefs? 
regs_test %>%
  head(head_val) %>% 
  mutate(coefs_plm = map(reg_plm, ~tidy(., quick = TRUE)),
         coefs_felm = map(reg_felm, ~tidy(., quick = TRUE))) %>%
  unnest(coefs_plm, coefs_felm) %>%
  filter(term !=term1)

broom:::tidy.lm(regs_test[2,]$reg_plm[[1]], quick = TRUE)

## test withins 
withins <- df_all %>%
  head(head_val) %>%
  mutate(within_plm = map(data_plm_M, ~ plm:::Within.matrix(., effect = "individual")),
         within_dplyr = map(data, ~ Within1_dplyr(select(., -year), individual)),
         within_felm = map(data, ~ demeanlist(.[,1:3], list(f1= factor(.$individual)))),
         within_dt_1 = map(data, ~ Within1_dt_1(., fe1="individual", x_names = c("A", "B1", "B2"))),
         within_dt_2 = map(data, ~ Within1_dt_2(., fe1="individual", x_names = c("A", "B1", "B2"))),
         within_dt_3a = map(data, ~ Within1_dt_3a(.[,1:4])),
         within_dt_3b = map(data, ~ Within1_dt_3b(.[,1:4])),
         within_dt_3c = map(data, ~ Within1_dt_3c(.[,1:4]))
         )

## Checks
map2_lgl(withins$within_plm, withins$within_dplyr, 
         ~all.equal(as_tibble(.x), .y, check.attributes=FALSE))


map2_lgl(withins$within_dplyr, withins$within_dt_1, 
         ~all.equal(as_tibble(.x), .y, check.attributes=FALSE))

map2_lgl(withins$within_dplyr, withins$within_dt_2, 
     ~all.equal(as_tibble(.x), .y, check.attributes=FALSE))

map2_lgl(withins$within_dplyr, withins$within_dt_3a, 
         ~all.equal(as_tibble(.x), .y, check.attributes=FALSE))

map2_lgl(withins$within_dplyr, withins$within_dt_3b, 
         ~all.equal(as_tibble(.x), .y, check.attributes=FALSE))

map2_lgl(withins$within_dplyr, withins$within_dt_3c, 
         ~all.equal(as_tibble(.x), .y, check.attributes=FALSE))


map2_lgl(withins$within_felm, withins$within_dplyr, 
         ~all.equal(.x, .y, check.attributes=FALSE))

##################################################
########## Benchmarks
##################################################

## regression benchmark
benches_reg <- df_all %>%
  head(head_val) %>%
  mutate(ben = map2(data, data_plm, 
                    ~microbenchmark(plm = plm(A~B1 + B2, data=.y),
                                    lfe = felm(A~B1 + B2|individual, data=.x),
                                    times=20)))

benches_reg$ben

## within benchmark
benches_within <- df_all %>%
  head(head_val) %>%
  mutate(ben = map2(data, data_plm_M, 
                    ~microbenchmark(plm = plm:::Within.matrix(.y, effect = "individual"),
                                    lfe = demeanlist(.x[,1:3], list(f1= factor(.x$individual))),
                                    dplyr = Within1_dplyr(select(.x, -year), individual),
                                    dt_1 = Within1_dt_1(.x, fe1="individual", x_names = c("A", "B1", "B2")),
                                    dt_2 = Within1_dt_2(.x, fe1="individual", x_names = c("A", "B1", "B2")),
                                    dt_3a = map(data, ~ Within1_dt_3a(.[,1:4])),
                                    dt_3b = map(data, ~ Within1_dt_3b(.[,1:4])),
                                    dt_3c = map(data, ~ Within1_dt_3c(.[,1:4])),
                                    times=50)))
benches_within$ben

## assemble benchmarks
bench_df <- data_frame(type = c("regression", "within"),
                       data = list(benches_reg, benches_within)) %>%
  mutate(data = map(data, ~ mutate(., ben = map(ben, ~summary(., unit="ms") %>%
                                                  mutate(expr = as.character(expr)))))) %>%
  unnest(data) %>%
  select(type, data_type, ben) %>%
  unnest(ben) %>%
  mutate(data_type = factor(data_type, ordered=TRUE, levels = head(df_all$data_type, head_val)))

benches_reg$ben[[1]]
benches_within$ben[[1]]

bench_df


## assemble results
# benches_all <- bench_df %>%
#   select(type, data_summary)  %>%
#   unnest(data_summary)
# benches_all <- rbind(benches_reg,
#                      benches_within) %>%
#   mutate(data_type = factor(data_type, ordered=TRUE, levels = head(df_all$data_type, head_val)))

## export
write_csv(bench_df, "benchmark_output.csv")

bench_df %>%
  ggplot(aes(x= data_type, y = mean, colour = expr))+
  geom_point() +
  # ylim(c(0, 200)) +
  facet_grid(type~., scales= "free") +
  ylab("Time (milliseconds)")


# ## Manually
# data_S <- data_small
# data_S_plm_M <- df_all$data_plm_M[[1]]
# bench_with_S <- microbenchmark(plm = plm:::Within.matrix(data_S_plm_M, effect = "individual"),
#                                lfe = demeanlist(data_S[,1:3], list(f1= factor(data_S$individual))),
#                                dplyr = Within1_dplyr(select(data_S, -year), individual),
#                                dt = Within1_dt_1(data_S, fe1="individual", x_names = c("A", "B1", "B2")),
#                                times=10)
# 
# bench_with_S
# 
# bench_reg_S <- microbenchmark(plm = plm(A~B1 + B2, data=df_all$data_plm[[1]]),
#                                lfe = felm(A~B1 + B2|individual, data=data_S),
#                                times=10)
# 
# bench_reg_S
# 
# 
# data_L1 <- df_all$data[[2]]
# data_L1_plm_M <- df_all$data_plm_M[[2]]
# bench_with_L1 <- microbenchmark(plm = plm:::Within.matrix(data_L1_plm_M, effect = "individual"),
#                                 lfe = demeanlist(data_L1[,1:3], list(f1= factor(data_L1$individual))),
#                                 dplyr = Within1_dplyr(select(data_L1, -year), individual),
#                                 dt = Within1_dt_1(data_L1, fe1="individual", x_names = c("A", "B1", "B2")),
#                                 times=10)
# 
# bench_reg_L1 <- microbenchmark(plm = plm(A~B1 + B2, data=df_all$data_plm[[2]]),
#                               lfe = felm(A~B1 + B2|individual, data=data_L1),
#                               times=10)
# 
# bench_S
# print(bench_with_S, unit = "ms")
# print(bench_reg_S, unit = "ms")
# print(bench_with_L1, unit = "ms")
# print(bench_reg_L1, unit = "ms")
# print(bench_reg_L1, unit = "s")
# 
# 
# ## compare
# ## regression benchmark
# benches_reg <- df_all %>%
#   head(head_val) %>%
#   mutate(ben = map2(data, data_plm, 
#                     ~microbenchmark(plm = plm(A~B1 + B2, data=.y),
#                                     lfe = felm(A~B1 + B2|individual, data=.x),
#                                     times=50)))
# 
# 
# bench_reg_S <- microbenchmark(plm = plm(A~B1 + B2, data=df_all$data_plm[[1]]),
#                               lfe = felm(A~B1 + B2|individual, data=data_S),
#                               times=50)
# 
# print(benches_reg$ben[[1]], unit = "ms")
# bench_reg_S
