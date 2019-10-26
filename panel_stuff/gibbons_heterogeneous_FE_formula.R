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
                          summarise_at(.x, c("gsp", "pcap", "pc"), list(var=var)))) %>% 
  ungroup()

regs_indiv

################################
#'## Reconstruct weighting
################################

## variances
with(regs_indiv, weighted.mean(estimate, pcap_var))
coef(reg_plm_bal)


################################
#'## More than one regressor
################################

## global
reg_plm_2reg <- plm(gsp ~ pcap + pc, data = Produc, index = c("state","year"))

## frish waugh
data_fw <- Produc_wit_bal %>% 
  mutate(pcap_res = residuals(lm(pcap~pc, data=.)),
         gsp_res = residuals(lm(gsp~pc, data=.)))

reg_plw_fw_2reg <- lm(gsp_res~ pcap_res-1, data = data_fw)
all.equal(coef(reg_plm_2reg)[1], coef(reg_plw_fw_2reg), check.attributes = FALSE)


##
regs_indiv_2reg <- data_fw %>% 
  group_by(state) %>% 
  group_modify(~bind_cols(broom::tidy(lm(gsp_res~ pcap_res-1, data=.x), quick=TRUE) %>% 
                            spread(term, estimate),
                          summarise_at(.x, c("gsp", "pcap", "pcap_res"), list(var=var)))) %>% 

  ungroup()

coef_hetero_2reg <- with(regs_indiv_2reg, weighted.mean(pcap_res, pcap_res_var))

all.equal(coef(reg_plm_2reg)[1], coef_hetero_2reg, check.attributes = FALSE)

