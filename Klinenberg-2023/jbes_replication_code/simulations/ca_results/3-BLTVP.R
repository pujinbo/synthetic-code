#https://scunning.com/cunningham_mixtape.pdf
# 3 is CA
# Treatment 1988 (19th)
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) #set working directory to source file
rm(list=ls()) #clear everything
devtools::install_github('johnson-shuffle/mixtape') #download data thank you cunningham
load("ca_comp.rdata")
source(here::here("01-model", "bitto_model_revised.R")) #upload model
library(mixtape) 
bl_dat <- mixtape::smoking
library(haven)
library(tidyverse)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
bl_dat <- bl_dat %>% 
  dplyr::select(year, state, cigsale) %>% 
  type_convert() %>% 
  pivot_wider(names_from = state, values_from = cigsale) %>% 
  dplyr::select(-year) %>% 
  dplyr::select("3", everything())

t_0 <- 19
t<-31

y<- as.vector(bl_dat$`3`)
x<- as.matrix(bl_dat[,-1])

y.bl   <- bltvp_bitto(y,
                      as.data.frame(x), 
                      niter = 10000,
                      nburn = 5000,
                      t_0 = t_0,
                      t=t) 
y.pred <- y.bl$y

save.image('ca_comp.rdata')

