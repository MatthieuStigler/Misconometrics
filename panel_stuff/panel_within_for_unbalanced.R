## objective: Check if within transformation is fine with balanced or unbalanced panel?

library(plm)
library(tidyverse)



data("Produc", package = "plm")
Produc <-  as_tibble(Produc)
Produc_unbal <- Produc[-c(1, 2, 25, 87),]

## balancedness?
Produc %>% 
  count(state) %>% 
  count(n)

Produc_unbal %>% 
  count(state) %>% 
  count(n)


# maual within
Produc_wit_bal <- Produc %>% 
  group_by(state, region) %>% 
  mutate_at(vars(-year, -group_cols()), ~. - mean(.)) %>% 
  ungroup()

Produc_wit_unbal <- Produc_unbal %>% 
  group_by(state, region) %>% 
  mutate_at(vars(-year, -group_cols()), ~. - mean(.)) %>% 
  ungroup()

################################
#'## Regressions
################################

## plm
reg_plm_bal <- plm(gsp ~ pcap, data = Produc, index = c("state","year"))
reg_plm_unbal <- plm(gsp ~ pcap, data = Produc_unbal, index = c("state","year"))

## lm
reg_lm_bal <- lm(gsp ~ pcap-1, data=Produc_wit_bal)
reg_lm_unbal <- lm(gsp ~ pcap-1, data=Produc_wit_unbal)

## compare coefs
coef(reg_plm_bal)==coef(reg_lm_bal)
coef(reg_plm_unbal)==coef(reg_lm_unbal)
