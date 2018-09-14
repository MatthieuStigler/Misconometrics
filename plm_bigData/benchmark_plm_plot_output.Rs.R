library(tidyverse)

levels <- c("data_small", "data_medium_N", "data_medium_T", "data_large_N", "data_large_T", "data_fat_N", "data_fat_N_fat_T")

bench_df <- read_csv("benchmark_output.csv") %>%
  mutate(data_type = factor(data_type, ordered=TRUE, levels = levels))
  

pl_all <- bench_df %>%
  ggplot(aes(x= data_type, y = mean, colour = expr))+
  geom_point() +
  # ylim(c(0, 200)) +
  facet_grid(type~., scales= "free") +
  ylab("Time (milliseconds)") +
  theme(axis.text = element_text(angle = 20)) +
  ggtitle("Benchmark: comparing regression and within transformations")

pl_within <- bench_df %>%
  filter(type == "within") %>%
  ggplot(aes(x= data_type, y = mean, colour = expr))+
  geom_point() +
  ylab("Time (milliseconds)") +
  ylim(c(0, 2600)) +
  theme(axis.text = element_text(angle = 20)) +
  ggtitle("Benchmark: within transformations, 5 algo at  data sizes")

pl_within

pl_all
pl_within

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
