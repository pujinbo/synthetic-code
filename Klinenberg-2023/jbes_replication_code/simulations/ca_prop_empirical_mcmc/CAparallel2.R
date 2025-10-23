#https://scunning.com/cunningham_mixtape.pdf
# 3 is CA
# Treatment 1988 (19th)
rm(list=ls()) #clear everything

# loading data ------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) #set working directory to source file
devtools::install_github('johnson-shuffle/mixtape') #download data thank you cunningham
if(!("pacman" %in% installed.packages())){
  install.packages("pacman")
}
pacman::p_load(
  "tidysynth",
  "mixtape",
  "janitor",
  "haven",
  "tidyverse",
  "CausalImpact",
  "bayesreg",
  "beepr",
  "Synth",
  "shrinkTVP",
  "glmnet", #for ArCo
  "rstan",
  # parallelize code
  "parallel",
  "foreach",
  "doParallel",
  "mlr3",
  "doSNOW"
)
bl_dat <- mixtape::smoking

# cleaning data -----------------------------------------------------------
setwd("../../01-model") #Moving to get my model
source("bitto_model_revised.R") #upload model
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
bl_dat <- tidysynth::smoking %>% 
  dplyr::select(year, state, cigsale) %>% 
  filter(state!="California") %>% 
  type_convert() %>% 
  tidyr::pivot_wider(names_from = state, values_from = cigsale) %>% 
  dplyr::select(-year) 

bl_dat<-bl_dat %>% 
  janitor::clean_names()



# loading time ------------------------------------------------------------------
t_0 <- 19
t<-31

#bltvp
max.steps=2000
burnin=1000

sims_per_obs<-1000
bl_dat_sim<-bl_dat


# Creating dataframe ------------------------------------------------------

placebo_dat<-read_csv(here::here("04-sim","ca_prop_empirical_mcmc", paste0("notvp_placebo_",sims_per_obs,".csv"))) %>% 
  try()

if(any(class(placebo_dat)=="try-error")){
  placebo_dat<-tibble(
    "model"=NA,
    "period"=NA,
    "sim"=NA,
    "unit"=NA,
    "dgp_tvp"=NA,
    "sim_y"=NA,
    "prediction"=NA,
    "lower_95"=NA,
    "upper_95"=NA
  )
}


# running the dataset -----------------------------------------------------
done_sims<-placebo_dat %>% 
  distinct(model,unit) %>% 
  mutate(done=1)

sims_to_do<-expand.grid(c("BLTVP",
                          "CI",
                          "CI-TVP",
                          "SC",
                          "ArCo",
                          "DM-LFM",
                          "BSCM-Horseshoe"),
                        sort(colnames(bl_dat))
                        ) %>% 
  as_tibble() %>% 
  rename("model"=Var1,
         "unit"=Var2) %>% 
  left_join(done_sims) %>% 
  filter(is.na(done)) %>% 
  dplyr::select(-done) %>% 
  arrange(model)

rm(placebo_dat)
for(i in 1:nrow(sims_to_do)){
  set.seed(which(colnames(bl_dat)==sims_to_do$unit[i]))
  form<-as.formula(paste0(sims_to_do$unit[i],"~."))
  m<-bayesreg::bayesreg(formula = form,
                        bl_dat,
                        model = "normal",
                        prior = "lasso",
                        burnin = 1000,
                       n.samples = sims_per_obs
  )
  sim<-as.matrix(bl_dat[,-which(colnames(bl_dat)==sims_to_do$unit[i])])%*%as.matrix(m$beta)+matrix(1,nrow=31,ncol=1)%*%m$beta0
  
  
  
  
  # for(k in 1:sims_per_obs){
  cl <- makeCluster(6)
  registerDoSNOW(cl)  # use multicore, set to the number of our cores
  cat(sims_to_do$unit[i], sims_to_do$model[i],paste0("(",i,"/",nrow(sims_to_do),")"), "started at",lubridate::date(),"...")
  

# CausalImpact ------------------------------------------------------------

  if(sims_to_do$model[i]=="CI"){
    cat("Performing CI...")
    s<-foreach (k=1:sims_per_obs, 
                .combine=rbind,
                
                .packages = c(
                  "tidyverse",
                  "shrinkTVP",
                  "mlr3",
                  "bayesreg",
                  "CausalImpact",
                  "ArCo",
                  "bpCausal",
                  "glmnet",
                  "matrixStats",
                  "Synth"
                )) %dopar% {
                  simk<-sim[,k]
                  bl_dat_sim[,which(colnames(bl_dat)==sims_to_do$unit[i])]<-sim[,k]
                  # creating CI estimate
                  ci<-CausalImpact::CausalImpact(bl_dat_sim %>% 
                                                   relocate(which(colnames(bl_dat)==sims_to_do$unit[i]),.before = 1),
                                                 pre.period = c(1,t_0),
                                                 post.period = c((t_0+1),t),
                                                 model.args = list(niter=max.steps)
                  )
                  cii<-tibble("model"="CI",
                              "period"=1:t,
                              "sim"=k,
                              "unit"=colnames(bl_dat)[which(colnames(bl_dat)==sims_to_do$unit[i])],
                              "dgp_tvp"="no",
                              "sim_y"=simk,
                              "prediction"=as.vector(ci$series$point.pred),
                              "lower_95"=as.vector(ci$series$point.pred.lower),
                              "upper_95"=as.vector(ci$series$point.pred.upper)
                  )
                  
                }

# CI-TVP ------------------------------------------------------------------

    
    }else if(sims_to_do$model[i]=="CI-TVP"){
      cat("Performing CI-TVP...")
                  s<-foreach (k=1:sims_per_obs, 
                              .combine=rbind,
                              
                              .packages = c(
                                "tidyverse",
                                "shrinkTVP",
                                "mlr3",
                                "bayesreg",
                                "CausalImpact",
                                "ArCo",
                                "bpCausal",
                                "glmnet",
                                "matrixStats",
                                "Synth"
                              )) %dopar% {
                                simk<-sim[,k]
                                bl_dat_sim[,which(colnames(bl_dat)==sims_to_do$unit[i])]<-sim[,k]
                                
                                ci_tvp<-CausalImpact::CausalImpact(bl_dat_sim %>% 
                                                                     relocate(which(colnames(bl_dat)==sims_to_do$unit[i]),.before = 1),
                                                                   pre.period = c(1,t_0),
                                                                   post.period = c((t_0+1),t),
                                                                   model.args = list(dynamic.regression=T,
                                                                                     niter=max.steps)
                                )
                                ciitvp<-tibble("model"="CI-TVP",
                                               "period"=1:t,
                                               "sim"=k,
                                               "unit"=colnames(bl_dat)[which(colnames(bl_dat)==sims_to_do$unit[i])],
                                               "dgp_tvp"="no",
                                               "sim_y"=simk,
                                               "prediction"=as.vector(ci_tvp$series$point.pred),
                                               "lower_95"=as.vector(ci_tvp$series$point.pred.lower),
                                               "upper_95"=as.vector(ci_tvp$series$point.pred.upper)
                                )
                              }

# BL-TVP ------------------------------------------------------------------

    
    }else if(sims_to_do$model[i]=="BLTVP"){
      cat("Performing BLTVP...")
      s<-foreach (k=1:sims_per_obs, 
                  .combine=rbind,
                  
                  .packages = c(
                    "tidyverse",
                    "shrinkTVP",
                    "mlr3",
                    "bayesreg",
                    "CausalImpact",
                    "ArCo",
                    "bpCausal",
                    "glmnet",
                    "matrixStats",
                    "Synth"
                  )) %dopar% {
                    simk<-sim[,k]
                    bl_dat_sim[,which(colnames(bl_dat)==sims_to_do$unit[i])]<-sim[,k]
                    y.bl   <- bltvp_bitto(as.matrix(bl_dat_sim[,which(colnames(bl_dat)==sims_to_do$unit[i])]),
                                          as.data.frame(as.matrix(bl_dat_sim[,-which(colnames(bl_dat)==sims_to_do$unit[i])])), 
                                          niter = max.steps,
                                          nburn = burnin,
                                          t_0 = t_0,
                                          t=t) 
                    
                    bltvp<-tibble("model"="BLTVP",
                                  "period"=1:t,
                                  "sim"=k,
                                  "unit"=colnames(bl_dat)[which(colnames(bl_dat)==sims_to_do$unit[i])],
                                  "dgp_tvp"="no",
                                  "sim_y"=simk,
                                  "prediction"=rowMedians(y.bl$y),
                                  "lower_95"=rowQuantiles(y.bl$y,probs = .025),
                                  "upper_95"=rowQuantiles(y.bl$y,probs = .975))
                    
                    
                  }

# SC ----------------------------------------------------------------------

      
    }else if(sims_to_do$model[i]=="SC"){
      cat("Performing SC...")
      s<-foreach (k=1:sims_per_obs, 
                  .combine=rbind,
                  
                  .packages = c(
                    "tidyverse",
                    "shrinkTVP",
                    "mlr3",
                    "bayesreg",
                    "CausalImpact",
                    "ArCo",
                    "bpCausal",
                    "glmnet",
                    "matrixStats",
                    "Synth"
                  )) %dopar% {
                    simk<-sim[,k, drop=FALSE]
                    bl_dat_sim[,which(colnames(bl_dat)==sims_to_do$unit[i])]<-sim[,k, drop=FALSE]
                    # Creating Synthetic Control
                    sc_y<-simk
                    sc_x<-as.matrix(bl_dat_sim[,-which(colnames(bl_dat)==sims_to_do$unit[i]), drop=FALSE])
                    Z0 = Z1 = X0 = X1 = NULL # clearing stuff up
                    
                    X0 = as.matrix( sc_x[c(1:t_0),])
                    X1 = as.matrix( sc_y[c(1:t_0)] )
                    Z0 = as.matrix( sc_x[c(1:t_0),] )
                    Z1 = as.matrix( sc_y[c(1:t_0)])
                    
                    sc_fit = synth( X1=X1 , 
                                    X0=X0 , 
                                    Z1=Z1 , 
                                    Z0=Z0#,
                                    # optimxmethod='All'
                                    ) %>% 
                      try()
                    # sc_fit = try(synth( X1=X1 , X0=X0 , Z1=Z1 , Z0=Z0 ),TRUE)
                    if(class(sc_fit)=="try-error"){
                      sc_pred<-rep(NA,31)
                    }else{
                    sc_weights = sc_fit$solution.w
                    sc_pred = sc_x%*%sc_weights
                    }
                    scc<-tibble("model"="SC",
                                "period"=1:t,
                                "sim"=k,
                                "unit"=colnames(bl_dat)[which(colnames(bl_dat)==sims_to_do$unit[i])],
                                "dgp_tvp"="no",
                                "sim_y"=simk,
                                "prediction"=sc_pred,
                                "lower_95"=sc_pred,
                                "upper_95"=sc_pred
                    )
                  }

# ArCo --------------------------------------------------------------------

    }else if(sims_to_do$model[i]=="ArCo"){
      cat("Performing ArCo...")
      s<-foreach (k=1:sims_per_obs, 
                  .combine=rbind,
                  
                  .packages = c(
                    "tidyverse",
                    "shrinkTVP",
                    "mlr3",
                    "bayesreg",
                    "CausalImpact",
                    "ArCo",
                    "bpCausal",
                    "glmnet",
                    "matrixStats",
                    "Synth"
                  )) %dopar% {
                    simk<-sim[,k]
                    bl_dat_sim[,which(colnames(bl_dat)==sims_to_do$unit[i])]<-simk
                    # Arco
                    
                    arco_dat<-list("data"=as.matrix(bl_dat_sim))
                    arco_fit<-ArCo::fitArCo(data = arco_dat,# getting the data
                                            treated.unit = which(colnames(arco_dat[[1]])==sims_to_do$unit[i]), # getting right row of the data
                                            fn=cv.glmnet, #using LASSO
                                            p.fn = predict, # recommended prediction
                                            t0=20, # when treatment begins
                                            VCOV.type = "nw", # preset used in help file
                                            boot.cf = TRUE, #bootstrapping
                                            R=100
                    )
                    
                    arco_a<-as.data.frame(arco_fit$boot.cf)
                    arco_lb<-matrixStats::rowQuantiles(as.matrix(arco_a),probs=.025)
                    arco_ub<-matrixStats::rowQuantiles(as.matrix(arco_a),probs=.975)
                    
                    arco_data<-tibble(
                      "model"="ArCo",
                      "period"=1:nrow(arco_dat$data),
                      "sim"=k,
                      "unit"=colnames(bl_dat)[which(colnames(bl_dat)==sims_to_do$unit[i])],
                      "dgp_tvp"="no",
                      "sim_y"=simk,
                      "prediction"=c(arco_fit$fitted.values,arco_fit$cf),
                      "lower_95"=c(rep(NA,length(arco_fit$fitted.values)),arco_lb),
                      "upper_95"=c(rep(NA,length(arco_fit$fitted.values)),arco_ub)
                    )
                    
                  }

# DM-LFM ---------------------------------------------------------------

    }else if(sims_to_do$model[i]=="DM-LFM"){
      cat("Performing DMLFM...")
      s<-foreach (k=1:sims_per_obs, 
                  .combine=rbind,
                  
                  .packages = c(
                    "tidyverse",
                    "shrinkTVP",
                    "mlr3",
                    "bayesreg",
                    "CausalImpact",
                    "ArCo",
                    "bpCausal",
                    "glmnet",
                    "matrixStats",
                    "Synth"
                  # )) %do% {
                )) %dopar% {
                    set.seed(k) # set seed so I can recreate maybe?
                    simk<-sim[,k, drop=FALSE]
                    bl_dat_sim[,which(colnames(bl_dat)==sims_to_do$unit[i])]<-sim[,k, drop=FALSE]
                    arco_dat<-list("data"=as.matrix(bl_dat_sim, drop=FALSE))
                    # BpCausal
                    bp_dat<-arco_dat$data %>% 
                      as_tibble() %>% 
                      mutate(t=1:nrow(.)) %>% 
                      pivot_longer(cols=-t, names_to="state",values_to="value") %>% 
                      mutate(treated=ifelse(state==colnames(bl_dat)[which(colnames(bl_dat)==sims_to_do$unit[i])] & t>=19,1,0)) %>% 
                      as.data.frame(., drop=FALSE)
                    
                    bpc<-bpCausal::bpCausal(data=bp_dat,
                                            index=c("state","t"),
                                            Yname="value",
                                            Dname="treated",
                                            Xname = NULL,
                                            Zname = NULL,
                                            Aname = NULL,
                                            r=10,
                                            re = "both",
                                            flasso = 1
                    ) %>% 
                      suppressMessages()
                    eout1 <- bpCausal::effSummary(bpc,   ## summary treatment effects
                                                  usr.id = NULL, ## treatment effect for individual treated units, if input NULL, calculate average TT
                                                  cumu = FALSE,  ## whether to calculate culmulative treatment effects
                                                  rela.period = TRUE) ## whether to use time relative to the occurence of treatment (1 is the first post-treatment period) or real period (like year 1998, 1999, ...)
                    bp_dat<-tibble(
                      "model"="DM-LFM",
                      "period"=1:31,
                      "sim"=k,
                      "unit"=colnames(bl_dat)[which(colnames(bl_dat)==sims_to_do$unit[i])],
                      "dgp_tvp"="no",
                      "sim_y"=simk,
                      "prediction"=eout1$est.eff$estimated_counterfactual,
                      "lower_95"=eout1$est.eff$counterfactual_ci_u, #looks like the naming is a bit off. upper is lower for this package
                      "upper_95"=eout1$est.eff$counterfactual_ci_l
                    )
                  }

# BSCM-Horseshoe ----------------------------------------------------------

    }else if(sims_to_do$model[i]=="BSCM-Horseshoe"){
      s<-foreach (k=1:sims_per_obs, 
                  .combine=rbind,
                  
                  .packages = c(
                    "tidyverse",
                    "shrinkTVP",
                    "mlr3",
                    "bayesreg",
                    "CausalImpact",
                    "ArCo",
                    "bpCausal",
                    "glmnet",
                    "matrixStats",
                    "Synth"
                  )) %dopar% {
                    simk<-sim[,k]
                    bl_dat_sim[,which(colnames(bl_dat)==sims_to_do$unit[i])]<-sim[,k]
                    # bayesreg
                    br<-bayesreg::bayesreg(formula(paste0(sims_to_do$unit[i],"~.")),
                                           data = bl_dat_sim %>% 
                                             slice(1:t_0),
                                           prior = "hs")
                    br_pred<-predict(br,
                                     bl_dat_sim %>% 
                                       dplyr::select(-sims_to_do$unit[i]),
                                     bayes.avg = T,
                                     CI=95
                    )
                    
                    bscm_dat<-tibble(
                      "model"="BSCM-Horseshoe",
                      "period"=1:31,
                      "sim"=k,
                      "unit"=colnames(bl_dat)[which(colnames(bl_dat)==sims_to_do$unit[i])],
                      "dgp_tvp"="no",
                      "sim_y"=simk,
                      "prediction"=br_pred[,1],
                      "lower_95"=br_pred[,3], #looks like the naming is a bit off. upper is lower for this package
                      "upper_95"=br_pred[,4]
                    )
                  }
  }
  

  cat("done at",lubridate::date(),"\n")
  s %>% 
    as.matrix() %>% 
    as_tibble() %>% 
    type_convert() %>%
    write_csv(paste0("notvp_placebo_",sims_per_obs,".csv"), append = T)
}
stopCluster(cl)
