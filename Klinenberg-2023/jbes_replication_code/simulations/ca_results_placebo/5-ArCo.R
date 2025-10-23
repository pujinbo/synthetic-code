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
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) #set working directory to source
load("ca_comp.rdata")

# recreating data
arco_dat <- tidysynth::smoking %>% 
  dplyr::select(year, state, cigsale) %>% 
  type_convert() %>% 
  tidyr::pivot_wider(names_from = state, values_from = cigsale) %>% 
  dplyr::select(-year) %>% 
  janitor::clean_names()
# cleaning data for arco
arco_dat<-list("data"=as.matrix(arco_dat))

#lasso for arco

arco_fit<-ArCo::fitArCo(data = arco_dat,# getting the data
                        treated.unit = which(colnames(arco_dat[[1]])=="california"), # getting right row of the data
                        fn=cv.glmnet, #using LASSO
                        p.fn = predict, # recommended prediction
                        t0=15, # when treatment begins
                        VCOV.type = "nw", # preset used in help file
                        boot.cf = TRUE, #bootstrapping
                        R=100
)
a<-as.data.frame(arco_fit$boot.cf)
lb<-matrixStats::rowQuantiles(as.matrix(a),probs=.025)
ub<-matrixStats::rowQuantiles(as.matrix(a),probs=.975)

mean(arco_fit$cf-arco_dat$data[15:31,18])
save.image('ca_comp.rdata')
