# Library
rm(list = ls())
devtools::install_github('johnson-shuffle/mixtape')
load(here::here("04-sim","samar_CA_prop", "ca_comp.rdata"))
library(mixtape)
ca_dat <- mixtape::smoking
library(tidyverse)
library(Synth)

ca_data <- ca_dat %>%
  dplyr::select(state, year, cigsale) %>% 
  pivot_wider(names_from = state,values_from = cigsale) %>% 
  dplyr::select(`3`, everything()) %>% 
  dplyr::select(-year)

ca <- ca_dat %>% 
  filter(state==3) %>% 
  dplyr::select(cigsale) %>% 
  type_convert()


T1 <- 19
T0 <- 31
n1 <- max(ca_dat$state)-1

# Outcome matrix
z1 = as.matrix(ca_dat$cigsale[ca_dat$state=='3'])
z0 = ca_dat$cigsale[ca_dat$state!='3']; z0 = matrix(z0,T0,n1)
nations = factor(unique( ca_dat$state[ca_dat$state!='3']  ))

#  implement sc 3 times  
Z0 = Z1 = X0 = X1 = NULL
sc_weights = rep(NA,n1)
sc_r       = rep(NA,n1) # for the empirical test


# X-matrix to match based on the outcome in the pre-intervention period
X0 = as.matrix( z0[c(1:T1),] )
X1 = as.matrix( z1[c(1:T1),] )
Z0 = as.matrix( z0[c(1:T1),] )
Z1 = as.matrix( z1[c(1:T1),] )


# Fit the SC algorithm
sc_fit = synth( X1=X1 , X0=X0 , Z1=Z1 , Z0=Z0 )  
sc_weights = sc_fit$solution.w


# Prediction
sc_pred = rep(NA,T0)
sc_pred = z0%*%sc_weights


# save the workspace
save.image('ca_comp.rdata')




