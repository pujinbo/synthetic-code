#####MASTER FILE#######
basepath<-"~/Dropbox/Synthetic Controls/9. replication code/version5"
codepath<-paste0(basepath,"/auxfiles")
inputpath<- basepath
outputpath<- paste0(basepath,"/figures")

#######SIMULATION OUTPUT####:
###NOTE: THESE FILES REQUIRE GUROBI TO RUN
##YOU MAY COMMENT THESE OUT AND REPRODUCE FIGURES USING THE INCLUDED OUTPUT

RUNSIM<-0

setwd(basepath)
if(RUNSIM==1) source("SC_application.R") #requires gurobi, but set up to run locally
if(RUNSIM==2) source("SC_EmpiricalARadjust_SNOW_byK.R") #requires gurobi, set up to run on linux cluster





######PLOTTING RESULTS#####
setwd(basepath)
source("plots_illustrative.R")
setwd(basepath)
source("plots_empiricalMC_descriptive.R")
setwd(basepath)
source("plots_empiricalMC_results.R")
setwd(basepath)
source("plots_application.R")
