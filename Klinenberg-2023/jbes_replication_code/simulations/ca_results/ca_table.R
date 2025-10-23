# Calculating the ATE for the Germany example
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(CausalImpact)
library(dplyr)
library(tidyverse)
library(tidyr)
library(matrixStats)
library(ggplot2)
load("ca_comp.rdata")

# BLTVP
bltvp <- c("BL-TVP",
           round(mean(rowMeans((y.pred-y))[t_0:T0]),1),
           round(mean(rowQuantiles((y.pred-y), probs = (.025))[t_0:T0]),1),
           round(mean(rowQuantiles((y.pred-y), probs = (.975))[t_0:T0]),1)
)

ci <- c("CI",
        round(mean((cim_pred$`keep$point.pred`-y)[t_0:T0]),1),
        round(mean((cim_pred$`keep$point.pred.lower`-y)[t_0:T0]),1),
        round(mean((cim_pred$`keep$point.pred.upper`-y)[t_0:T0]),1)
)

ci_tvp <- c("CI-TVP",
        round(mean((cim_pred_tvp$`keep_tvp$point.pred`-y)[t_0:T0]),1),
        round(mean((cim_pred_tvp$`keep_tvp$point.pred.lower`-y)[t_0:T0]),1),
        round(mean((cim_pred_tvp$`keep_tvp$point.pred.upper`-y)[t_0:T0]),1)
)

sc <- c("SC",
        round(mean((sc_pred-y)[t_0:T0]),1),
        "NA",
        "NA"
)

dmlfm <- c("DM-LFM",
        round(mean(-1*dmlfm$estimated_ATT[t_0:T0]),1),
        round(mean(-1*dmlfm$estimated_ATT_ci_l[t_0:T0]),1),
        round(mean(-1*dmlfm$estimated_ATT_ci_u[t_0:T0]),1)
)

arco<-c("ArCo",round(-1*c(arco_fit$delta),1))

bayesreg_tab <- c("BSCM-Horseshoe",
           round(mean(-1*(y-bayesreg_pred_dat$pred)[t_0:T0]),1),
           round(mean(-1*(y-bayesreg_pred_dat$lower_95)[t_0:T0]),1),
           round(mean(-1*(y-bayesreg_pred_dat$upper_95)[t_0:T0]),1)
)


res <- tibble(
  "Model" = c(bltvp[1],ci[1],ci_tvp[1],sc[1],dmlfm[1],arco[1],bayesreg_tab[1]),
  "Average GDP Decrease" = c(bltvp[2],ci[2],ci_tvp[2],sc[2],dmlfm[2],arco[3],bayesreg_tab[2]),
  "2.5th Percentile GDP Decrease" = c(bltvp[3],ci[3],ci_tvp[3],sc[3],dmlfm[4],arco[4],bayesreg_tab[3]),
  "97.5th Percentile GDP Decrease"= c(bltvp[4],ci[4],ci_tvp[4],sc[4],dmlfm[3],arco[2],bayesreg_tab[4])
)

# Table 5
write.csv(res, "../pictures_tables/ca_ate.csv")
