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
               "matrixStats",
               "bayesreg",
               "rstan"
)
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) #set working directory to source
load("ca_comp.rdata")

# recreating data
bayesreg_dat <- tidysynth::smoking %>% 
  dplyr::select(year, state, cigsale) %>% 
  type_convert() %>% 
  tidyr::pivot_wider(names_from = state, values_from = cigsale) %>% 
  dplyr::select(-year) %>% 
  janitor::clean_names()


#bayesreg

bayes_reg_setup<-list("N_train"=t_0,
                      "N_test"=t-t_0,
                      "p"=38,
                      "y_train"=bayesreg_dat$california[1:t_0],
                      "X_train"=bayesreg_dat[1:t_0,] %>% 
                        dplyr::select(-california),
                      "X_test"=bayesreg_dat[(t_0+1):t,] %>% 
                        dplyr::select(-california)
)

test1<-stan(
  file = here::here("model", "bscm_horseshoe.stan"),
  data = bayes_reg_setup,
  chains = 4,
  warmup = 1000,
  iter = 2000,
  cores = 1,
  refresh=0
)

bayesreg_pred_dat<-summary(test1)$summary[which(str_detect(rownames(summary(test1)$summary),"y_")),c(1,4,8)] %>% 
  as_tibble() %>% 
  clean_names() %>% 
  rename("pred"="mean",
         "lower_95"="x2_5_percent",
         "upper_95"="x97_5_percent")

# running with bayesreg which *may* work but i'm not sure
# bayesreg_fit<-bayesreg::bayesreg(california~.,
#                                  data = bayesreg_dat[1:t_0,])
# 
# sim<-as.matrix(bayesreg_dat[1:t_0,which(colnames(bayesreg_dat)!="california")])%*%as.matrix(bayesreg_fit$beta)+matrix(1,nrow=t_0,ncol=1)%*%bayesreg_fit$beta0
# 
# bayesreg_pred<-predict(bayesreg_fit,
#                        newdata = bayesreg_dat[t_0+1:t,],
#                        bayes.avg = TRUE,
#                        CI=95)

# bayesreg_pred_dat<-tibble("pred"=c(rowQuantiles(sim,probs=.5),bayesreg_pred[,1]),
#                           "lower_95"=c(rowQuantiles(sim,probs=.025),bayesreg_pred[,3]),
#                           "upper_95"=c(rowQuantiles(sim,probs=.975),bayesreg_pred[,4])
#                           )

save.image('ca_comp.rdata')
