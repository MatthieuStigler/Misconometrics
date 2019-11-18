library(lme4)
library(plm)
library(tidyverse)
library(broom)

library(glmnet)
library(fastDummies)

################################
#'## Create data
################################

N <- 200
T <- 50

set.seed(123)
x_means <- rnorm(N, mean=3)
x_sigmas <- rlnorm(N, sdlog = 0.5)
x <- rnorm(N*T, mean = rep(x_means, each=T),
           sd= rep(x_sigmas, each=T))
alphas <- rnorm(N)
betas <-  rnorm(N, mean=1)
y <-  rep(alphas, each=T) + rep(betas, each=T) * x + rnorm(N*T)

alpha_sd <- sd(alphas)
betas_sd <- sd(betas)

##
df <- tibble(x=x, y=y, indiv = rep(seq_len(N), each=T), time = rep(seq_len(T), N))
df_plm <- pdata.frame(df, index = c("indiv", "time"))

################################
#'## Estimate
################################

fo_rename <- function(df) {
  df %>% 
    rename(intercept = `(Intercept)`)  
}

tidy_wide <- function(df) {
  df %>% 
    tidy(quick=TRUE) %>% 
    select(term, estimate) %>% 
    spread(term, estimate) %>% 
    fo_rename()
}

term_clean <-  function(df) {
 df %>% 
    mutate(term = case_when(term=="(Intercept)"~"intercept",
                            TRUE ~ term)) %>% 
    arrange(term)
}

## lm 
reg_lm <- lm(y~x, data=df)
coefs_lm <- reg_lm%>% 
  tidy(quick=TRUE ) %>% 
  rename(mean=estimate) %>% 
  term_clean

coefs_lm

### plm
reg_plm <- pvcm(y~x, data=df_plm, model = "random")
reg_plm_ols <- pvcm(y~x, data=df_plm, model = "within")
coef_plm_ols <- as_tibble(coef(reg_plm_ols)) %>% 
  fo_rename() %>% 
  gather(term, estimate, everything()) %>% 
  group_by(term) %>% 
  summarise(mean = mean(estimate),
            sd = sd(estimate)) %>% 
    ungroup()

coef_plm_ols


coef_plm <- enframe(coef(reg_plm), "term", "mean") %>% 
  mutate(sd = sqrt(diag(reg_plm$Delta))) %>% 
  term_clean

coef_plm

## lmer
reg_lme <- lmer(y~x+ (x | indiv), data=df)
coef_lme <- reg_lme %>% 
  # df %>% 
  tidy(quick=TRUE) %>% 
  select(term, estimate) %>% 
  head(4) %>% 
  mutate(stat = if_else(str_detect(term, "^sd"), "sd", "mean"),
         term = str_remove_all(term, "^sd_|\\.indiv")) %>% 
  spread(stat, estimate) %>% 
  term_clean

coef_lme

################################
#'## Ridge
################################

X <- fastDummies::dummy_cols(df$indiv)
X <- model.matrix(y~-1+ind_fac*x, data=df %>%  mutate(ind_fac = factor(indiv)))
dim(X)
str_remove_all(colnames(X), "[0-9]") %>%  unique()
X[1:5, 1:4]

tryit <- lm.fit(X, df$y)
enframe(coef(tryit), value = "estimate") %>% 
  mutate(term = if_else(str_detect(name, "x$"), "x", "intercept"),
         indiv = str_extract(name, "[0-9]{1,3}") %>%  as.integer()) %>% 
  select(-name) %>% 
  group_by(term) %>% 
  summarise(mean = mean(estimate),
            sd = sd(estimate)) %>% 
  ungroup()



ridge.mod <- glmnet(x, y, alpha = 0, lambda = lambda)
predict(ridge.mod, s = 0, exact = T, type = 'coefficients')[1:6,]

################################
#'## Compare all
################################

compare <- bind_rows(coef_lme %>% 
            mutate(estimator = "mixed"),
          coef_plm %>% 
            mutate(estimator = "swamy"),
          coef_plm_ols %>% 
            mutate(estimator = "ols_each"),
          coefs_lm %>% 
            mutate(estimator = "ols")) %>% 
  select(estimator, everything())

compare


################################
#'## Visualize
################################

betas_all_df <- tibble(x=betas, estimator = "true") %>% 
  rbind(as_tibble(coef(reg_plm_ols)) %>% 
          select(x) %>% 
          mutate(estimator = "ols_each"))

betas_all_df

betas_all_df %>% 
  ggplot(aes(x=x, fill=estimator)) +
  geom_density(alpha = I(0.4)) +
  geom_vline(xintercept = mean(betas)) 
