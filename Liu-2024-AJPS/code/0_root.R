## Replication file: Liu, Wang, Xu (2022)
## install packages and set path

rm(list=ls())
path <- "~/Dropbox/ProjectZ/TSCS_placebo/submission_AJPS/replication"
setwd(path)

########################################
## install packages
########################################

# install "fastplm", "fect" and "panelView" from files (in order)

# install.packages("readstata13")
# install.packages("gridExtra")
# install.packages("ggplot2")
# install.packages("extraDistr")
# install.packages("foreach")
# install.packages("ggplot2")
# install.packages("future")
# install.packages("doParallel")
# install.packages("abind")

########################################
## Three Examples
########################################

## Simulated Example
# takes about 90 minutes.
source("code/1_ex_sim0.R")

## HH2015
# takes about 7 minutes.
source("code/2_ex_HH2015.R")

## FM2015
# takes about 25 minutes.
source("code/3_ex_FM2015.R")

########################################
## Monte Carlo Exercise
########################################

# takes about 1.2 hours.
source("code/4_sim_tests_n100.R")

# takes about 3 hours.
source("code/5_sim_tests_n300.R")


# takes about 16 hours.
source("code/6_sim_ife_mc.R")

# take about 30 minutes.
source("code/7_sim_inference.R")

#### Plot Monte Carlo Results ####

# almost instantly
source("code/8_plot_sim.R")


