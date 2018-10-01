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
