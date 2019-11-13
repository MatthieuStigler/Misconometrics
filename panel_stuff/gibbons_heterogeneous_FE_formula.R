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

# maual within demean
Produc_demeaned_bal <- Produc %>% 
  group_by(state, region) %>% 
  mutate_at(vars(-year, -group_cols()), ~. - mean(.)) %>% 
  ungroup()

## create unbalanced data
Produc_unbal <- Produc[-c(1, 2, 25, 17*3+c(1,6, 8), 17*5+c(1,3, 12, 14)),]
count(Produc_unbal, state) %>%  
  count(n)

Produc_demeaned_unbal <- Produc_unbal %>% 
  group_by(state, region) %>% 
  mutate_at(vars(-year, -group_cols()), ~. - mean(.)) %>% 
  ungroup()


################################
#'## Helper function
################################

coefs_by_group <- function(df, .groups=vars(state, region), formula = gsp ~ pcap) {
  df %>% 
    group_by_at(.groups) %>% 
    group_modify(~bind_cols(summarise(.x, n_vals = n()),
                            broom::tidy(lm(formula, data=.x), quick=TRUE) %>% 
                              filter(str_detect(term, "^pcap"))%>% 
                              mutate(term = paste0("coef_", term)) %>%
                              spread(term, estimate),
                            summarise_at(.x, tidyselect::vars_select(names(.), starts_with("pcap")), list(var=var)))) %>% 
    ungroup()
}

################################
#'## Regressions
################################

## global
reg_plm_bal <- plm(gsp ~ pcap, data = Produc, index = c("state","year"))


## individual
regs_indiv <- Produc_demeaned_bal %>% 
  coefs_by_group

regs_indiv 

################################
#'## Reconstruct weighting
################################

## is same? yes!!
coef_bal_manu <- with(regs_indiv, weighted.mean(coef_pcap, pcap_var))
coef_bal_plm <- coef(reg_plm_bal)[["pcap"]]

all.equal(coef_bal_manu, coef_bal_plm)

################################
#'## More than one regressor
################################

## global
reg_plm_2reg <- plm(gsp ~ pcap + pc, data = Produc, index = c("state","year"))

## frish waugh
data_fw <- Produc_demeaned_bal %>% 
  mutate(pcap_res = residuals(lm(pcap~pc, data=.)),
         gsp_res = residuals(lm(gsp~pc, data=.)))

reg_plw_fw_2reg <- lm(gsp_res~ pcap_res-1, data = data_fw)
all.equal(coef(reg_plm_2reg)[["pcap"]], coef(reg_plw_fw_2reg)[["pcap_res"]])


## Estimate individual coefficients
regs_indiv_2reg <- data_fw %>% 
  coefs_by_group(formula = gsp_res~ pcap_res-1)

regs_indiv_2reg


## average coefs, compare to plm
coef_hetero_2reg <- with(regs_indiv_2reg, weighted.mean(coef_pcap_res, pcap_res_var))

all.equal(coef(reg_plm_2reg)[["pcap"]], coef_hetero_2reg)

################################
#'## Higher level groups?
################################

## regions
Produc_demeaned_bal %>% 
  count(region)

## Let's estimate panel by region
plm_byReg <- Produc_demeaned_bal %>% 
  coefs_by_group(.groups = vars(region)) %>% 
  rename(coef_plm_byReg = coef_pcap) %>% 
  select(-pcap_var)

plm_byReg

## weight manual by groups
manu_byReg <- regs_indiv %>% 
  group_by(region) %>% 
  summarise(coef_manu_byReg = weighted.mean(coef_pcap, pcap_var),
            pcap_var =sum(pcap_var)) %>% 
  left_join(plm_byReg, by = "region") %>% 
  mutate(is_same = abs(coef_manu_byReg - coef_plm_byReg)< 0.0000000000001)

manu_byReg  

## weight to get full
coef_manu_fromReg <- manu_byReg  %>% 
  summarise(coef_general = weighted.mean(coef_manu_byReg, pcap_var)) %>% 
  pull(coef_general) %>% .[[1]]

all.equal(coef_manu_fromReg, coef_bal_plm)



################################
#'## Unalanced data
################################


## unbal: global level
coef_plm_unbal <- coef(plm(gsp ~ pcap, data = Produc_unbal, index = c("state","year")))[["pcap"]]

## unbal: reg by state
regs_indiv_unbal <- Produc_unbal %>% 
  coefs_by_group()

regs_indiv_unbal

## Full OLS as re-weight individuals
coef_manu_unbal <- regs_indiv_unbal %$% 
  weighted.mean(coef_pcap, (n_vals-1) * pcap_var)

coef_plm_unbal

## compare global versus weightted individual
all.equal(coef_plm_unbal, coef_manu_unbal)

#### At higher level ####

### individual region regs 
regs_indiv_unbal_byReg <- Produc_demeaned_unbal %>% 
  coefs_by_group(.groups = vars(region))

## reweight states to regions
manu_byReg_unbal <- regs_indiv_unbal %>% 
  group_by(region) %>% 
  summarise(coef_manu_byReg = weighted.mean(coef_pcap, (n_vals-1) * pcap_var),
            pcap_var_sum =sum((n_vals-1) *pcap_var)) %>% 
  left_join(regs_indiv_unbal_byReg, by = "region") %>% 
  mutate(is_same = abs(coef_manu_byReg - coef_pcap)< 0.00000000001)

manu_byReg_unbal


## weight to get full
coef_manu_fromReg <- manu_byReg_unbal  %>% 
  summarise(coef_general = weighted.mean(coef_manu_byReg, pcap_var_sum)) %>% 
  pull(coef_general) %>% .[[1]]

all.equal(coef_manu_fromReg, coef_plm_unbal)
