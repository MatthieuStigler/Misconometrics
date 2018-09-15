library(tidyverse)

levels <- c("data_small", "data_medium_N", "data_medium_T", "data_large_N", "data_large_T",
            "data_fat_N", "data_fat_N_med_T")

bench_W_df <- read_csv("benchmark_output_centerOnly.csv") %>%
  mutate(data_type = factor(data_type, ordered=TRUE, levels = levels)) %>%
  select(data_type, n_row, do_dplyr, expr, mean, everything()) %>%
  group_by(data_type) %>%
  mutate(perc_min = (mean-min(mean))/min(mean)) %>%
  ungroup()

bench_W_df

## dt3 comparisons
bench_W_df %>%
  filter(str_detect(expr, "dt_3")) %>%
  group_by(data_type) %>%
  mutate(perc_min = (mean-min(mean))/min(mean)) %>%
  ungroup() %>%
  ggplot(aes(x= data_type, y = perc_min, colour = expr))+
  geom_point(alpha=I(.6), size = 2, fill = I("black"))


pl_compars_centerOnly <- bench_W_df %>%
  filter(!expr %in% c("dplyr", "dt_3a", "dt_3b")) %>%
  group_by(data_type) %>%
  mutate(perc_min = (mean-min(mean))/min(mean)) %>%
  ungroup() %>%
  gather(measure, value, mean, perc_min) %>%
  ggplot(aes(x= data_type, y = value, colour = expr))+
  geom_point(alpha=I(.4), size = 4, fill = I("black")) +
  ylab("Time (percent or absolut [millisec]))") +
  facet_grid(measure~., scales = "free") +
  theme(axis.text = element_text(angle = 10)) +
  ggtitle("Benchmark: center transformations at different data sizes (100 reps)")

pl_compars_centerOnly

ggsave(pl_compars_centerOnly, filename = "figures/benchm_center_dt_lfe.png",
       width = 7, height = 4)

