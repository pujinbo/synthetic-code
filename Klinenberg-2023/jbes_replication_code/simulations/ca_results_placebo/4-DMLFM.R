rm(list=ls())
# install.packages('devtools', repos = 'http://cran.us.r-project.org') # if not already installed
# devtools::install_github('liulch/bpCausal')
if(!("pacman" %in% installed.packages())){
  install.packages("pacman")
}
library(pacman)
pacman::p_load("here",
               "tidyverse",
               "mixtape",
               "janitor",
               "beepr",
               "ArCo",
               "glmnet", #for ArCo
               "matrixStats"
)
library(bpCausal)
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) #set working directory to source
load("ca_comp.rdata")

# recreating data
bl_dat <- tidysynth::smoking %>% 
  dplyr::select(year, state, cigsale) %>% 
  # filter(state!="California") %>% 
  type_convert()

bl_dat<-bl_dat %>% 
  janitor::clean_names() %>% 
  mutate(treat=ifelse(state=="California" & year>=1985,1,0)) %>% 
  as.data.frame()

set.seed(5)

bpc<-bpCausal::bpCausal(data=bl_dat,
                   index=c("state","year"),
                   Yname="cigsale",
                   Dname="treat",
                   Xname = NULL,
                   Zname = NULL,
                   Aname = NULL,
                   r=10,
                   re = "both",
                   flasso = 1
                   )

eout1 <- effSummary(bpc,   ## summary treatment effects
                    usr.id = NULL, ## treatment effect for individual treated units, if input NULL, calculate average TT
                    cumu = FALSE,  ## whether to calculate culmulative treatment effects
                    rela.period = TRUE) ## whether to use time relative to the occurence of treatment (1 is the first post-treatment period) or real period (like year 1998, 1999, ...)

dmlfm<-eout1$est.eff %>% 
  as_tibble() %>% 
  mutate(n=sort(unique(bl_dat$year))) 

save.image('ca_comp.rdata')

