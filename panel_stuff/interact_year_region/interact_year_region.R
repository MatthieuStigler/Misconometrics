library(plm)
library(lfe)
library(broom)
library(tidyverse)

data("Produc", package = "plm")
zz <- plm(gsp ~ pcap, data = Produc, index = c("state","year"), effect = "twoway")
summary(zz)

## year region
Produc$year_region <- paste(Produc$year, Produc$region, sep = "_")
zz2 <- plm(gsp ~ pcap, data = Produc, index = c("state","year_region"), effect = "twoway")


zz
zz2

FE_1 <- fixef(zz, effect = "time")
FE_2 <- fixef(zz2, effect = "time")
FE_df_1 <- tibble(fe = names(FE_1), value = FE_1) %>% 
  mutate(year = as.integer(fe))
FE_df_2 <- tibble(fe = names(FE_2), value = FE_2) %>% 
  separate(fe, c("year", "region"), remove = FALSE, convert = TRUE)
FE_df_2

## ave
FE_df_2_ave <- FE_df_2 %>% 
  group_by(year) %>% 
  summarise(value = mean(value)) %>% 
  ungroup()

FE_df_both <- rbind(select(FE_df_1, -fe),
                    FE_df_2_ave) %>% 
  mutate(type = rep(c("FE_1", "FE_2"), each = nrow(FE_df_1)))

FE_df_both

## plot
FE_df_both %>% 
  ggplot(aes(x = year, y = value, colour = type))+
  geom_point()


#####################
### lfe
#####################

## year
felm <- felm(gsp ~ pcap|state+year, data = Produc)
all.equal(tidy(felm) %>% as.data.frame(), 
          tidy(zz)%>% as.data.frame())

## year region
felm2 <- felm(gsp ~ pcap|state+year_region, data = Produc)
all.equal(tidy(felm2) %>% as.data.frame(), 
          tidy(zz2)%>% as.data.frame())

Produc2 <- Produc %>% 
  mutate(region = as.factor(region),
         year = as.factor(year))
felm2b <- felm(gsp ~ pcap|state+year:region, data = Produc2)

all.equal(tidy(felm2) %>% as.data.frame(), 
          tidy(felm2b)%>% as.data.frame())

