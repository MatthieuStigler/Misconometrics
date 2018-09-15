library(data.table)
library(broom)
library(lfe)
library(microbenchmark)
library(tidyverse)

## tidyverse
center_dplyr <- function(x, vari) {
  vari <- enquo(vari)
  x %>%
    group_by(!!vari) %>%
    mutate_all(funs(. - mean(.))) %>%
    ungroup() %>%
    select(-!!vari)
}

## data table
center_dt_1 <- function(x, fe1, x_names = NULL) {

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


center_dt_2 <- function(x, fe1, x_names = NULL) {
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

center_dt_3a <- function(x) {
  setDT(x,key="cell")
  x_fullnames <- names(x)
  mean_x <- x[,lapply(.SD,mean, na.rm=TRUE),by=cell]
  setcolorder(mean_x,x_fullnames)
  setkey(mean_x,cell)
  res<-x-mean_x[x,][,1:length(names(mean_x))]
  res[,cell:=NULL]
  setDF(res)
  return(res)
}

center_dt_3b <- function(t) {
  setDT(t,key="cell")
  x_fullnames <- names(t)
  x_names <- x_fullnames[which(x_fullnames != "cell")]
  mean_x <- t[,lapply(.SD,mean,na.rm=TRUE),by=cell]
  res<-t[mean_x,.(A=x.A-i.A,B1=x.B1-i.B1,B2=x.B2-i.B2),on="cell"]
  setDF(res)
  return(res)
}

center_dt_3c <- function(x) {
  setDT(x,key="cell")
  mean_x <- x[, lapply(.SD,mean, na.rm=TRUE), by=cell][,cell:=NULL][x$cell,]
  x[,cell:=NULL]
  res <- x - mean_x 
  setDF(res)
  return(res)
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
                       cell = rep(1:N, each = T),
                       year = rep(1:T, N))

data_sim

data_small <- filter(data_sim, cell %in% 1:300 & year %in% 1:3)
data_medium_N <- filter(data_sim, cell %in% 1:10000 & year %in% 1:3)
data_medium_T <- filter(data_sim, cell %in% 1:300 & year %in% 1:20)
data_large_N <- filter(data_sim, cell %in% 1:100000 & year %in% 1:3)
data_large_T <- filter(data_sim, cell %in% 1:10000 & year %in% 1:20 )
data_fat_N <- filter(data_sim, year %in% 1:3)
data_fat_N_med_T <- filter(data_sim, year %in% 1:10)



## data frame of data frames

df_all <- data_frame(data_type = c("data_small",  "data_medium_N",  "data_medium_T", "data_large_N", "data_large_T", "data_fat_N", "data_fat_N_med_T"),
                     data = list(data_small, data_medium_N, data_medium_T, data_large_N, data_large_T, data_fat_N, data_fat_N_med_T),
                     do_dplyr = str_detect(data_type, "small|medium")) %>%
  mutate(n_row = map_int(data, ~nrow(.) %>% as.integer)) %>%
  select(data_type, n_row, everything())


##################################################
########## Check results same
##################################################

head_val <- 5 ## if interactive testing or large R CMD BATCH
times_bench <- 10

## test center_out 
center_out <- df_all %>%
  head(head_val) %>%
  mutate(center_dplyr = map(data, ~ center_dplyr(select(., -year), cell)),
         center_felm = map(data, ~ demeanlist(.[,1:3], list(f1= factor(.$cell)))),
         center_dt_1 = map(data, ~ center_dt_1(., fe1="cell", x_names = c("A", "B1", "B2"))),
         center_dt_2 = map(data, ~ center_dt_2(., fe1="cell", x_names = c("A", "B1", "B2"))),
         center_dt_3a = map(data, ~ center_dt_3a(.[,1:4])),
         center_dt_3b = map(data, ~ center_dt_3b(.[,1:4])),
         center_dt_3c = map(data, ~ center_dt_3c(.[,1:4]))
         )

## Checks
check_out <- center_out %>%
  gather(method, center_i, contains("_dt_"), contains("_felm")) %>%
  mutate(check_equal = map2_lgl(center_dplyr, center_i, 
           ~all.equal(as.data.frame(.x), as.data.frame(.y), check.attributes=FALSE)))

check_out %>%
  filter(!check_equal)

##################################################
########## Benchmarks
##################################################

## center benchmark
benchm_center <- df_all %>%
  # head(head_val) %>%
  # head(2) %>%
  mutate(ben = map2(data, do_dplyr, 
                    ~microbenchmark(lfe = demeanlist(.x[,1:3], list(f1= factor(.x$cell))),
                                    dplyr =if(.y) center_dplyr(select(.x, -year), cell) else NA,
                                    dt_1 = center_dt_1(.x, fe1="cell", x_names = c("A", "B1", "B2")),
                                    dt_2 = center_dt_2(.x, fe1="cell", x_names = c("A", "B1", "B2")),
                                    dt_3a = center_dt_3a(.x[,1:4]),
                                    dt_3b = center_dt_3b(.x[,1:4]),
                                    dt_3c = center_dt_3c(.x[,1:4]),
                                    times=100)))
benchm_center$ben

## assemble benchmarks
bench_df <-   benchm_center %>%
  mutate(ben = map(ben, ~summary(., unit="ms") )) %>%
  select(-data) %>%
  unnest(ben)

bench_df



write_csv(bench_df, "benchmark_output_centerOnly.csv")

