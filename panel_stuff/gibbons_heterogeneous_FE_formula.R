## objective: Check if within transformation is fine with balanced or unbalanced panel?

library(plm)
library(magrittr)
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
  group_by(state, region) %>% 
  group_modify(~bind_cols(broom::tidy(lm(gsp ~ pcap-1, data=.x), quick=TRUE),
                          summarise_at(.x, c("gsp", "pcap", "pc"), list(var=var)))) %>% 
  ungroup()

regs_indiv

################################
#'## Reconstruct weighting
################################

## is same? yes!!
coef_bal_manu <- with(regs_indiv, weighted.mean(estimate, pcap_var))
coef_bal_plm <- coef(reg_plm_bal)[["pcap"]]

all.equal(coef_bal_manu, coef_bal_plm)

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
all.equal(coef(reg_plm_2reg)[["pcap"]], coef(reg_plw_fw_2reg)[["pcap_res"]])


## Estimate individual coefficients
regs_indiv_2reg <- data_fw %>% 
  group_by(state, region) %>% 
  group_modify(~bind_cols(broom::tidy(lm(gsp_res~ pcap_res-1, data=.x), quick=TRUE) %>% 
                            spread(term, estimate),
                          summarise_at(.x, c("gsp", "pcap", "pcap_res"), list(var=var)))) %>% 

  ungroup()

## average coefs, compare to plm
coef_hetero_2reg <- with(regs_indiv_2reg, weighted.mean(pcap_res, pcap_res_var))

all.equal(coef(reg_plm_2reg)[["pcap"]], coef_hetero_2reg)

################################
#'## Higher level groups?
################################

## regions
Produc_wit_bal %>% 
  count(region)

## Can I re-weight groups?
plm_byReg <- Produc %>% 
  group_by(region) %>% 
  group_modify(~broom::tidy(plm(gsp ~ pcap, index = c("state","year"), 
                                data = .))) %>% 
  ungroup() %>% 
  select(region, estimate) %>% 
  rename(coef_plm_byReg = estimate)

plm_byReg

## weight manual by groups
manu_byReg <- regs_indiv %>% 
  group_by(region) %>% 
  summarise(coef_manu_byReg = weighted.mean(estimate, pcap_var),
            pcap_var =sum(pcap_var)) %>% 
  left_join(plm_byReg, by = "region") %>% 
  mutate(is_same = abs(coef_manu_byReg - coef_plm_byReg)< 0.0000000000001)

manu_byReg  

## weight to get full
coef_manu_fromReg <- manu_byReg  %>% 
  summarise(coef_general = weighted.mean(coef_manu_byReg, pcap_var)) %>% 
  pull(coef_general) %>% .[[1]]

all.equal(coef_manu_fromReg, coef_bal_plm)
