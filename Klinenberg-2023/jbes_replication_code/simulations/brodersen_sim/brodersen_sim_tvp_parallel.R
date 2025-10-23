# Simulation using brodersen setup -----------------------------------

rm(list=ls())
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
               "bpCausal",
               "Synth",
               "rstan",
               # parallelize code
               "parallel",
               "foreach",
               "doParallel"
)
source(here::here("model","bitto_model_revised.R")) #upload model

# Generate data -----------------------------------------------------------

# Using the same length as Abadie 2010
set.seed(1234)
t_0 <- 19 # Treatment occurs
t<-31 # total length
len<-1:t # total length for making x1 and x2
num_sims<-1000 #takes about 4 minutes for 10 simulations aka 24 seconds per simulation.

#bltvp
max.steps=20000
burnin=1000

registerDoParallel(cores=6)  # use multicore, set to the number of our cores
s<-foreach (j=1:num_sims, .combine=rbind) %dopar% {
  dat<-tibble(
    "model"=NA,
    "period"=len,
    "sim"=j,
    "x1"=sin(pi*len/12),
    "beta1"=1,#.1,
    "sqrt_theta1"=2,#.2,
    "betatilde1"=NA,
    "x2"=cos(pi*len/6),
    "beta2"=1,#.1,
    "sqrt_theta2"=0,
    "betatilde2"=NA,
    "x3"=3+sin(len),
    "beta3"=0,
    "sqrt_theta3"=3,#.3,
    "betatilde3"=NA,
    "x4"=cos(len/2),
    "beta4"=0,
    "sqrt_theta4"=0,
    "betatilde4"=NA,
    "local_level"=NA,
    "error"=rnorm(t,mean=0,sd=.1^2),
    "y"=NA,
    "prediction"=NA,
    "lower_95"=NA,
    "upper_95"=NA
  )
  # from brodersen et al.
  dat$betatilde1[1]<-0
  dat$betatilde2[1]<-0
  dat$betatilde3[1]<-0
  dat$betatilde4[1]<-0
  dat$local_level[1]<-0
  
  
  # making the simulated y
  for(i in 2:t){
    dat$betatilde1[i]<-rnorm(1,mean=dat$betatilde1[i-1],1)
    dat$betatilde2[i]<-rnorm(1,mean=dat$betatilde2[i-1],1)
    dat$betatilde3[i]<-rnorm(1,mean=dat$betatilde3[i-1],1)
    dat$betatilde4[i]<-rnorm(1,mean=dat$betatilde4[i-1],1)
    dat$local_level[i]<-rnorm(1,mean=dat$local_level[i-1],.1^2)
  }
  dat$y<-dat$x1*(dat$beta1+dat$betatilde1*dat$sqrt_theta1)+dat$x2*(dat$beta2+dat$betatilde2*dat$sqrt_theta2)+dat$x3*(dat$beta3+dat$betatilde3*dat$sqrt_theta3)+dat$x4*(dat$beta4+dat$betatilde4*dat$sqrt_theta4)+dat$local_level+dat$error
  
  # creating the counterfactuals
  
  # CI
  ci<-CausalImpact::CausalImpact(dat %>% 
                                   dplyr::select(y,x1,x2,x3,x4),
                                 pre.period = c(1,t_0),
                                 post.period = c((t_0+1),t),
                                 model.args = list(niter=max.steps)
  )
  dat_ci<-dat %>% 
    mutate(prediction=as.numeric(ci$series$point.pred),
           lower_95=as.numeric(ci$series$point.pred.lower),
           upper_95=as.numeric(ci$series$point.pred.upper),
           model="CI"
           )
  #ci-tvp
  ci_tvp<-CausalImpact::CausalImpact(dat %>%
                                   dplyr::select(y,x1,x2,x3,x4),
                                 pre.period = c(1,t_0),
                                 post.period = c((t_0+1),t),
                                 model.args = list(niter=max.steps,
                                                   dynamic.regression=T)
  )
  big_dat<-rbind(dat_ci,
        dat %>%
          mutate(prediction=as.numeric(ci_tvp$series$point.pred),
                 lower_95=as.numeric(ci_tvp$series$point.pred.lower),
                 upper_95=as.numeric(ci_tvp$series$point.pred.upper),
                 model="CI-TVP"
          )
        )
  # arco
  ar_data<-dat %>%
    dplyr::select(y,x1,x2,x3,x4) %>%
    as.matrix(.)
  arco_data<-list("data"=ar_data)
  arco_fit<-ArCo::fitArCo(arco_data,
                treated.unit = 1,
                fn=cv.glmnet, #using LASSO
                p.fn = predict, # recommended prediction
                t0=19, # when treatment begins
                VCOV.type = "nw", # preset used in help file
                boot.cf = TRUE, #bootstrapping
                R=100
                )
  a<-as.data.frame(arco_fit$boot.cf)
  big_dat<-rbind(big_dat,dat %>%
    mutate(model="ArCo",
           "prediction"=c(arco_fit$fitted.values,arco_fit$cf),
           "lower_95"=c(rep(NA,t_0-1),matrixStats::rowQuantiles(as.matrix(a),probs=.025)),
           "upper_95"=c(rep(NA,t_0-1),matrixStats::rowQuantiles(as.matrix(a),probs=.975)))
  )

  #DM-LFM
  dm_dat<-dat %>%
    dplyr::select(period,x1,x2,x3,x4,y) %>%
    pivot_longer(-period, names_to="unit", values_to="values") %>%
    mutate(treated=ifelse(unit=="y" & period>=19,1,0)) %>%
    as.data.frame(.)


  bpc<-bpCausal::bpCausal(data=dm_dat,
                          index=c("unit","period"),
                          Yname="values",
                          Dname="treated",
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

  big_dat<-rbind(big_dat,
                 dat  %>%
                   mutate(model="DM-LFM",
                          prediction=eout1$est.eff$estimated_counterfactual,
                          "lower_95"=eout1$est.eff$counterfactual_ci_u, #looks like the naming is a bit off. upper is lower for this package
                          "upper_95"=eout1$est.eff$counterfactual_ci_l
                   )
  )
  # #BL-TVP

  bl_tvp<-bltvp_bitto(dat$y,
                  x=as.data.frame(cbind(dat$x1,dat$x2,dat$x3,dat$x4)),
                  t=t,
                  t_0=t_0,
                  niter=20000,
                  nburn=10000,
                  nthin=1)
  big_dat<-rbind(big_dat,
                 dat  %>%
                   mutate(model="BL-TVP",
                          prediction=rowMedians(bl_tvp$y),
                          "lower_95"=rowQuantiles(bl_tvp$y,probs=.025),
                          "upper_95"=rowQuantiles(bl_tvp$y,probs=.975)
                   )
  )

  cat("BL-TVP done...")

  # SC
  z1 = as.matrix(dat$y)
  z0 = as.matrix(dat %>%
                   dplyr::select(starts_with("x")))

  X0 = as.matrix( z0[c(1:t_0),] )
  X1 = as.matrix( z1[c(1:t_0),] )
  Z0 = as.matrix( z0[c(1:t_0),] )
  Z1 = as.matrix( z1[c(1:t_0),] )

  sc_fit = synth( X1=X1 , X0=X0 , Z1=Z1 , Z0=Z0 )
  sc_weights = sc_fit$solution.w
  sc_pred = z0%*%sc_weights

  big_dat<-rbind(big_dat,
                 dat %>%
                   mutate(model="SC",
                          prediction=sc_pred)
  )
  cat("SC done...")

# horseshoe bayesian thing ------------------------------------------------


  bayes_reg_setup<-list("N_train"=t_0,
                        "N_test"=t-t_0,
                        "p"=4,
                        "y_train"=dat$y[1:t_0],
                        "X_train"=dat[1:t_0,] %>%
                          dplyr::select(starts_with("x")),
                        "X_test"=dat[(t_0+1):t,] %>%
                          dplyr::select(starts_with("x"))
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

  big_dat<-rbind(big_dat,
                 dat %>%
                   mutate(model="bayesreg",
                          prediction=bayesreg_pred_dat$pred,
                          lower_95=bayesreg_pred_dat$lower_95,
                          upper_95=bayesreg_pred_dat$upper_95)
  )
}


# save output ------------------------------------------------------------------

s %>%
  write_rds(here::here("04-sim","brodersen_sim",paste0("brodersen_sim_tvp_",num_sims,".rds")))
