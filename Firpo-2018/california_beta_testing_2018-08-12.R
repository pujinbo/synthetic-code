###############################################################################
# Article: Synthetic Control Estimator: A Generalzied Inference Procedure and
#          Confidence Sets
# Authors: Sergio Firpo and Vitor Possebom
# Code by: Vitor Possebom
# Goal: We use the California exercise to beta test the function SCM.CS.
###############################################################################
# Clean and organize the work environment
rm(list = ls())
setwd('C:/Users/VitorPossebom/Dropbox/artigos_proprios/published_articles/SCE_GIP_CS/function_confidence_sets/R_files')
library("Synth")
library("doParallel")
# library("data.table")
# library("statar")
# library("dplyr")
source('function_SCM-CS_v07.R')
# Setup parallel backend to use all but one of cores.
n.cores <- detectCores()-1
cl <- makeCluster(n.cores)
registerDoParallel(cl)
###############################################################################
# Load and organize the data
smoking <- read.csv(file = "smoking_dataset.csv")
smoking <- smoking[order(smoking$year,smoking$state),]
stateid <- as.numeric(rep(1:39,31))
smoking <- cbind(smoking,stateid)
smoking$state <- rep(levels(smoking$state),31)
californiaid <- 3
#######################################
# Estimate the SC weights for each state
results <- foreach(j=1:39, .combine = rbind, .packages = "Synth") %dopar% {
  # Define the comparison regions.
  controlunits <- setdiff(1:39, j)
  # Prepare the data in order to use the synthetic control estimator
  dataprep.out <- dataprep(
    foo = smoking,
    predictors = c("lnincome","beer","age15to24","retprice"),
    predictors.op = "mean",
    time.predictors.prior = seq(from = 1980, to = 1988, by = 1),
    special.predictors = list(
      list("cigsale", seq(from = 1975, to = 1975, by = 1), "mean"),
      list("cigsale", seq(from = 1980, to = 1980, by = 1), "mean"),
      list("cigsale", seq(from = 1988, to = 1988, by = 1), "mean")),
    dependent = "cigsale",
    unit.variable = "stateid",
    unit.names.variable = "state",
    time.variable = "year",
    treatment.identifier = j,
    controls.identifier = controlunits,
    time.optimize.ssr = seq(from = 1970, to = 1988, by = 1),
    time.plot = seq(from = 1970, to = 2000, by = 1))
  # Estimating the synthetic control method
  synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS")
  Y <- dataprep.out$Y1plot
  weights <- synth.out$solution.w
  return(c(Y, NA, weights))
}
# Organize our results
Ymat <- t(results[, 1:31])
weightsmat <- t(results[, 33:70])
#######################################
# Stop parallel backend
stopCluster(cl)
##############################################################################
# Define the options of the function_SCM-CS.
treated <- californiaid
T0 <- 19
precision <- 30
type <- "linear"
phi <- 0
v <- matrix(0, 1, 39)
significance <- 4/39
plot <- TRUE
###############################################################################
# Compute the confidence sets using function_SCM-CS
bounds <- SCM.CS(Ymat, weightsmat, treated, T0, phi, v, precision, type,
                 significance, plot)
