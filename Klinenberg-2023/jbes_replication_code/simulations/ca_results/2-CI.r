# libraries and workspace
rm(list = ls())
devtools::install_github('johnson-shuffle/mixtape')
load("ca_comp.rdata")
library(mixtape)
ca_dat <- mixtape::smoking
library(tidyverse)
library(CausalImpact)

ca_data <- ca_dat %>%
  dplyr::select(state, year, cigsale) %>% 
  pivot_wider(names_from = state,values_from = cigsale) %>% 
  dplyr::select(`3`, everything()) %>% 
  dplyr::select(-year)

ca <- ca_dat %>% 
  filter(state==3) %>% 
  dplyr::select(cigsale) %>% 
  type_convert()


T1 <- 19 # 18 training periods
T0 <- 31


# CIM data frame
cim_data = as.data.frame(ca_data)
names(cim_data) = c('y',paste0('x',1:(ncol(ca_data)-1)))
cim_post = cim_data[(T1+1):T0,1]

keep <- CausalImpact::CausalImpact(cim_data,pre.period = c(1,T1),post.period = c(T1+1,T0))$series
cim_pred =  cbind(keep$point.pred , keep$point.pred.lower , keep$point.pred.upper)

keep_tvp <- CausalImpact::CausalImpact(cim_data,
                                       pre.period = c(1,T1-1),
                                       post.period = c(T1,T0),
                                       model.args = list(dynamic.regression=T))$series
cim_pred_tvp =  cbind(keep_tvp$point.pred , keep_tvp$point.pred.lower , keep_tvp$point.pred.upper)

save.image('ca_comp.rdata')

# ts.plot(as.vector(cim_pred$`keep$point.pred`), ylim=c(40,140), col="blue")
# lines(cim_pred$`keep$point.pred.lower`, col="blue")
# lines(cim_pred$`keep$point.pred.upper`, col="blue")
# lines(ca$cigsale, col="red")
# abline(v=18)
# 
