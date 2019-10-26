## objective: Check if within transformation is fine with balanced or unbalanced panel?

library(plm)
library(tidyverse)

data("Produc", package = "plm")
Produc <-  as_tibble(Produc)

## balancedness?
Produc %>% 
  count(state) %>% 
  count(n)

################################
#'## Prepare data
################################

# maual within
Produc_wit_bal <- Produc %>% 
  group_by(state, region) %>% 
  mutate_at(vars(-year, -group_cols()), ~. - mean(.)) %>% 
  ungroup()

################################
#'## Regressions
################################

## global
reg_plm_bal <- plm(gsp ~ pcap, data = Produc, index = c("state","year"))


## individual
regs_indiv <- Produc_wit_bal %>% 
  group_by(state) %>% 
  group_modify(~bind_cols(broom::tidy(lm(gsp ~ pcap-1, data=.x), quick=TRUE),
                          summarise_at(.x, c("gsp", "pcap"), list(var=var)))) %>% 
  ungroup()

regs_indiv

################################
#'## Reconstruct weighting
################################

## variances
with(regs_indiv, weighted.mean(estimate, pcap_var))
coef(reg_plm_bal)
