library(tidyverse)

levels <- c("data_small", "data_medium_N", "data_medium_T", "data_large_N", "data_large_T", "data_fat_N", "data_fat_N_med_T")

bench_df <- read_csv("benchmark_output.csv") %>%
  mutate(data_type = factor(data_type, ordered=TRUE, levels = levels))

##


## in percentage
bench_df_perc <- bench_df %>%
  group_by(type, data_type) %>%
  mutate(per_min = (mean-min(mean))/min(mean)) %>%
  ungroup()

pl_all <- bench_df %>%
  ggplot(aes(x= data_type, y = mean, colour = expr))+
  geom_point() +
  # ylim(c(0, 200)) +
  facet_grid(type~., scales= "free") +
  ylab("Time (milliseconds)") +
  theme(axis.text = element_text(angle = 20)) +
  ggtitle("Benchmark: comparing regression and within transformations")

pl_all

pl_within <- bench_df %>%
  filter(type == "within") %>%
  ggplot(aes(x= data_type, y = mean, colour = expr))+
  geom_point(alpha=I(.6), size = 2, fill = I("black")) +
  # geom_jitter(width = 0.1, height = 0) +C
  ylab("Time (milliseconds)") +
  ylim(c(0, 2600)) +
  theme(axis.text = element_text(angle = 20)) +
  ggtitle("Benchmark: within transformations, 5 algo at  data sizes")

pl_within

## percentage plot
max_i <- 8
bench_df_perc %>%
  filter(type == "within") %>%
  mutate(per_min =ifelse(per_min < max_i, per_min, max_i +rnorm(1, sd = 0.5))) %>%
  select(data_type, per_min, expr) %>%
  ggplot(aes(x= data_type, y = per_min, colour = expr))+
  geom_point(alpha=I(.6), size = 2, fill = I("black")) +
  # geom_jitter(width = 0.1, height = 0) +C
  ylab("TPercent above best") +
  ylim(c(0, max_i+1)) +
  theme(axis.text = element_text(angle = 20)) +
  ggtitle("Benchmark: within transformations, 5 algo at  data sizes")



pl_all
pl_within %+%
  filter(bench_df, type == "within" & !expr %in% c("dplyr", "dt_3a", "dt_3b")) +
  ylim(c(0, 400))

pl_within %+%
  filter(bench_df, type == "within" & str_detect(expr, "dt_3")) +
  ylim(c(0, 400))


  ggsave(pl_all, filename = "plot_benchmark reg and within.png", height = 4, width = 7)
ggsave(pl_within, filename = "plot_benchmark within.png", height = 4, width = 7)

## numbers
bench_df %>%
  # filter(type =="within") %>%
  group_by(data_type) %>%
  mutate(comp = 100 * (mean - min(mean))/min(mean),
         is_mean  = mean ==min(mean)) %>%
  ungroup() %>%
  filter(expr == "plm" | is_mean)
