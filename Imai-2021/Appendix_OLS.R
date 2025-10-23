rm(list=ls())

## please set the working directory
setwd("~/gitrepos/tscs/code/AJPS_tscs_replication/")

## check/create the folder to put outputs
if (file.exists(file.path(getwd(), "output"))){
    OUT_DIR <- file.path(getwd(), "output")
} else {
    dir.create(file.path(file.path(getwd(), "output")))
    OUT_DIR <- file.path(getwd(), "output")
}


## This script takes the environment saved in Getting_PanelEstimates.R
load("./PanelMatch_temp2.RData")

## loading packages
library(foreign)
library(DataCombine)
library(multiwayvcov)
library(lmtest)
library(stargazer)


## list of variables for creating lagged values
varnames <- c("y", "dem", "logpop", "tradewb", "Populationages014oftotal",
              "Populationages1564oftota", "nfagdp", "unrest")

## creating lagged values (1-4 time periods) for each vars
for(l in 1:4){
  lag <- l
  for(varname in varnames){
    newname <- paste(varname, "_l", lag, sep="")
    cmd <- paste("d2 <- slide(d2, Var = \"", varname, "\", NewVar = \"", newname, "\", slideBy = -",
                 as.character(lag), ", GroupVar = \"wbcode2\"", ", TimeVar = \"year\"", 
                 ")", sep="")
    
    eval(parse(text=cmd))
  }
}

two_way_fe_Ace_con <- lm(y ~ 
             democ + rever + 
               y_l1 + y_l2 + y_l3 + y_l4 + 
     Populationages014oftotal_l1 + Populationages014oftotal_l2 +
     Populationages014oftotal_l3 + Populationages014oftotal_l4 +
     Populationages1564oftota_l1 + Populationages1564oftota_l2 +
     Populationages1564oftota_l3 + Populationages1564oftota_l4 +
     unrest_l1 + unrest_l2 + 
     unrest_l3 + unrest_l4 +
     tradewb_l1 + tradewb_l2 +
     tradewb_l3 + tradewb_l4 + 
     nfagdp_l1 + nfagdp_l2 + nfagdp_l3 + nfagdp_l4 +
     as.factor(wbcode2) + as.factor(year),
   data = d2)

m_vcov <- cluster.vcov(two_way_fe_Ace_con, d2$wbcode2) # extracting the clustered variance covoaraince matrix

## outputs
out <- coeftest(two_way_fe_Ace_con, m_vcov) # coeftest

#
two_way_fe_Ace <- lm(y ~ 
                           democ + rever + 
                           y_l1 + y_l2 + y_l3 + y_l4 +
                           as.factor(wbcode2) + as.factor(year),
                         data = d2)

m_vcov <- cluster.vcov(two_way_fe_Ace, d2$wbcode2) # extracting the clustered variance covoaraince matrix

## outputs
out2 <- coeftest(two_way_fe_Ace, m_vcov) # coeftest


filename <- paste(OUT_DIR, "/Section_S6_Table_S6_1_Columns_1_3.tex", sep="")
sink(file.path(filename))
stargazer(two_way_fe_Ace, 
          #two_way_fe_Ace, 
          two_way_fe_Ace_con,
          #two_way_fe_Ace_con, 
          #font.size = "large", 
          digits = 4, 
          type = "latex",
         # dep.var.labels = "Growth",
          
          # covariate.labels = c("Land-specific arrests", "Land-specific arrests at t-1"),
          dep.var.labels.include = F,
          style = "qje", no.space = FALSE, 
          omit = c("Constant", "wbcode", "year",
                   "Populationages014oftotal", "Populationages1564oftota", "unrest", "tradewb", "nfagdp"),
          omit.yes.no = c("No", "Yes"),
          add.lines = list(c("country FE", "Yes", "Yes"),
                             #"Yes", "Yes"), 
                           c("time FE", "Yes","Yes"),
                             #"Yes","Yes"),
                           c("covariates", "No", "Yes"), 
                             #"Yes", "Yes"),
                           c("estimation", "OLS", "OLS"
                             )),
          # omit.labels = c("Prefecture_FE", "Time_FE"),
          # keep = c("tour"
          # ),
          covariate.labels = 
            c("ATT", "ART", 
              "$\\hat\\rho_1$",
              "$\\hat\\rho_2$","$\\hat\\rho_3$","$\\hat\\rho_4$"
              ),
          star.cutoffs = c(.10, .05, .01), 
          se = list(out2[,2], out[,2]),
          keep.stat = "N",
          header = F, notes = "robust standard errors clustered by prefecture in parentheses", 
          float = F)
sink()





