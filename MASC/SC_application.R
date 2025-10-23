library(data.table)
library(Synth) #includes the basque data application
library(gurobi)
library(snow)
library(doSNOW)
library(plyr)
library(parallel)
library(optimx)


# California covariates are: 1) log-GDP per cap, 2) retail price, 3) % age 15-24,  4) beer per capita, 5) cig sales per cap, 1988, 6) Cig sales per cap, 1980, 7) Cig sales per cap, 1975

setwd(codepath)
source("Estimator_Code.R")
source("Cross-validation_Code_byK.R")

#cl <- makeCluster(mpi.universe.size()-1, type="MPI",outfile="") #COMMENT OUT if running locally
cl <- makeCluster(4, type="SOCK")
registerDoSNOW(cl)


clusterEvalQ(cl,{
  library(data.table)
  library(Synth)
  library(gurobi)
  library(snow)
  library(doSNOW)
  #  library(Rmpi)
  library(plyr)
  library(RColorBrewer)
  library(kernlab)
  library(optimx)
})

##data and placebo setup:
data(basque)
basque<-as.data.table(basque)
basque<-basque[regionno!=1,]
basque[,regionname:= gsub(" (.*)","",regionname)]
#Adjusting Basque covariates as done by A&G:
# dataprep: prepare data for synth
dataprep.out <-
  dataprep(
    foo = basque
    ,predictors= c("school.illit",
                   "school.prim",
                   "school.med",
                   "school.high",
                   "school.post.high"
                   ,"invest"
    )
    ,predictors.op = c("mean")
    ,dependent     = c("gdpcap")
    ,unit.variable = c("regionno")
    ,time.variable = c("year")
    ,special.predictors = list(
      list("gdpcap",1960:1969,c("mean")),                            
      list("sec.agriculture",seq(1961,1969,2),c("mean")),
      list("sec.energy",seq(1961,1969,2),c("mean")),
      list("sec.industry",seq(1961,1969,2),c("mean")),
      list("sec.construction",seq(1961,1969,2),c("mean")),
      list("sec.services.venta",seq(1961,1969,2),c("mean")),
      list("sec.services.nonventa",seq(1961,1969,2),c("mean")),
      list("popdens",1969,c("mean")))
    ,treatment.identifier  = 17
    ,controls.identifier   = c(2:16,18)
    ,time.predictors.prior = c(1964:1969)
    ,time.optimize.ssr     = c(1960:1969)
    ,unit.names.variable   = c("regionname")
    ,time.plot            = c(1955:1997) 
  )

# 1. combine highest and second highest 
# schooling category and eliminate highest category
dataprep.out$X1["school.high",] <- 
  dataprep.out$X1["school.high",] + 
  dataprep.out$X1["school.post.high",]
dataprep.out$X1                 <- 
  as.matrix(dataprep.out$X1[
    -which(rownames(dataprep.out$X1)=="school.post.high"),])
dataprep.out$X0["school.high",] <- 
  dataprep.out$X0["school.high",] + 
  dataprep.out$X0["school.post.high",]
dataprep.out$X0                 <- 
  dataprep.out$X0[
    -which(rownames(dataprep.out$X0)=="school.post.high"),]

# 2. make total and compute shares for the schooling catgeories
lowest  <- which(rownames(dataprep.out$X0)=="school.illit")
highest <- which(rownames(dataprep.out$X0)=="school.high")

dataprep.out$X1[lowest:highest,] <- 
  (100 * dataprep.out$X1[lowest:highest,]) /
  sum(dataprep.out$X1[lowest:highest,])
dataprep.out$X0[lowest:highest,] <-  
  100 * scale(dataprep.out$X0[lowest:highest,],
              center=FALSE,
              scale=colSums(dataprep.out$X0[lowest:highest,])
  )

basque[,school.high:=school.high+school.post.high]
basque[,school.post.high:=NULL]

setwd(inputpath)

data<-list()
graphnames<-c("Spain")
names <-list()
names[[1]]<- c(unique(basque[regionno==17,regionname]),unique(basque[regionno!=17,regionname]))
data[[1]]<-list()
data[[1]]$outcome<-cbind(basque[regionno==17,gdpcap],t(reshape(basque[regionno!=17,.(regionno,year,gdpcap)], idvar='regionno', timevar='year',direction='wide')[,-"regionno",with=FALSE]))
data[[1]]$covariates<-copy(basque)
data[[1]]$covariates[,time:=year-min(year)+1]
data[[1]]$covariates[,year:=NULL]
data[[1]]$covariates[,regionno:=NULL]
data[[1]]$covariates[,popdens:=mean(popdens,na.rm=TRUE),by=regionname]
for(u in 1:length(names[[1]])){
data[[1]]$covariates[regionname==names[[1]][u],unit:=u]
}
data[[1]]$covariates[,regionname:=NULL]

data[[1]]$covariates<-data[[1]]$covariates[,c(8:11,13,1:7,12,14,15),with=FALSE]
data[[1]]$name<-names[[1]][1]
treatment_times<-list()
ADH_fittimes<-list()
treatment_times[[1]]<-16
ADH_fittimes[[1]]<-6:15




Kvals<-c(1:10)
philength<-100
omits<-list()

nameappend <- c("Sp_")
for(run in c(1)){
  omits[[run]]<-2:length(names[[run]])
  
  for(omit in omits[[run]]){
    control<-length(data)+1
    names[[control]]<-c(names[[run]][omit],names[[run]][-c(1,omit)])
    treatment_times[[control]]<-treatment_times[[run]]
    ADH_fittimes[[control]]<-ADH_fittimes[[run]]
    graphnames<-c(graphnames,paste0(nameappend[run],names[[run]][omit]))

    data[[control]]<-list()
    data[[control]]$outcome<-cbind(data[[run]]$outcome[,omit],data[[run]]$outcome[,-c(1,omit)])
    data[[control]]$covariates<-data[[run]]$covariates[unit!=1,]
    data[[control]]$covariates[unit==omit,unit:=1]
    data[[control]]$name<-names[[run]][omit]
  }
}


basedata<-data
exnames<-graphnames
##############################################
setwd(outputpath)


for(r in 1:length(data)){
  data[[r]]$donors<-basedata[[r]]$outcome[,-1]
  data[[r]]$treated<-as.matrix(basedata[[r]]$outcome[,1])
  data[[r]]$treatment<-treatment_times[[which(graphnames==exnames[r])]]
  data[[r]]$ADH_fittimes<-ADH_fittimes[[which(graphnames==exnames[r])]]
}

#Constructing matrix of outcomes for calculating variance weights
maincovariates<-copy(data[[1]]$covariates)
maincovariates<-maincovariates[time%in%data[[1]]$ADH_fittimes,]
maincovariates[,time:=time-min(time)+1]
allcovariates<-copy(maincovariates)


regularizers.dist=TRUE
regularizers.weight=TRUE

ex<- Filter(function(x) is.function(get(x,.GlobalEnv)),ls(.GlobalEnv))
clusterExport(cl,ex)

params<-foreach(run = 1:length(data))%dopar%{
  print(run)
  output<-list()
  
  dump.info<-TRUE

    #########ESTIMATORS, RESTRICTING TO ADH SAMPLE WINDOW
    ADHdata<-copy(data[[run]])
    ADHdata$outcome<-ADHdata$outcome[-(1:(min(ADHdata$ADH_fittimes)-1)),]
    ADHdata$donors<-ADHdata$donors[-(1:(min(ADHdata$ADH_fittimes)-1)),]
    ADHdata$treated<-as.matrix(ADHdata$treated[-(1:(min(ADHdata$ADH_fittimes)-1)),])
    ADHdata$treatment<-length(ADHdata$ADH_fittimes)+1
    ADHdata$covariates<-ADHdata$covariates[time%in%ADHdata$ADH_fittimes,]
    ADHdata$covariates[,time:=time-min(time)+1]
    ADHdata$firstforecast<-5
    maxphi<-1e5
    maxphi.abadie<-finddist(estimator=solve.locreg,
                            data=ADHdata,
                            tune.pars.list=list(phi=maxphi,K=1,
                                                set.k=NA,
                                                min.preperiods=ADHdata$firstforecast),
                            model='regularize',
                            only.synth=FALSE,only.match=FALSE,dump.info=FALSE)
    
    maxphi.abadiecov<-finddist(estimator=solve.covreg,
                            data=ADHdata,
                            tune.pars.joint=lapply(maxphi, function(x) list(
                              K=1, phi=x, set.k=NA, min.preperiods=ADHdata$firstforecast,
                              Vfun=Cov.Vars,
                              maincovariates=allcovariates,
                              type="all"
                            )),
                            tune.pars.list=NULL,
                            model='regularize',
                            est.options=NULL,
                            only.synth=FALSE,only.match=FALSE,dump.info=dump.info)
    

    
    output$ADH<-solve.synth(data=ADHdata,
                     est.options=list(sigf.ipop=5,
                                    Margin.ipop=5e-04,
                                    Wbar=Wbar),
                     tune.pars=list(phi=0,K=1))

    
    maxphi.abadiecovscweight<-finddist(estimator=solve.covreg,
                               data=ADHdata,
                               tune.pars.joint=lapply(maxphi, function(x) list(
                                 K=1, phi=x, set.k=NA, min.preperiods=ADHdata$firstforecast,
                                 scV=  unlist(output$ADH$tune.pars$V),
                                 Vfun=Cov.Vars,
                                 maincovariates=allcovariates,
                                 type="all"
                               )),
                               tune.pars.list=NULL,
                               model='regularize',
                               est.options=list(sigf.ipop=5,
                                                Margin.ipop=5e-04,
                                                Wbar=Wbar),
                               only.synth=FALSE,only.match=FALSE,dump.info=dump.info)
    
    
    
    output$matchADHAll2 <- cv.solver(estimator=solve.locreg,
                                    data=ADHdata,
                                    model='analytic',
                                    tune.pars.joint=lapply(Kvals,function(x) list(
                                      K=x, phi=1, set.k=NA, min.preperiods=ADHdata$firstforecast,
                                      Vfun=Cov.Vars,
                                      maincovariates=allcovariates,
                                      type="all"
                                    )),
                                    tune.pars.list=NULL,
                                    est.options=list(Wbar=solve.covmatch),
                                    only.synth=FALSE,only.match=TRUE,dump.info=dump.info)
    
    
    
    output$SCAM.CovAll2 <- cv.solver(estimator=solve.synth,
                                    data=ADHdata,
                                    model='analytic',
                                    tune.pars.joint=lapply(Kvals,function(x) list(
                                      K=x, phi=0, set.k=NA, min.preperiods=ADHdata$firstforecast,
                                      Vfun=Cov.Vars,
                                      maincovariates=allcovariates,
                                      type="all"
                                    )),
                                    tune.pars.list=NULL,
                                    est.options=list(Wbar=solve.covmatch,
                                                     sigf.ipop=5,
                                                     Margin.ipop=5e-04),
                                    only.synth=FALSE,only.match=FALSE,dump.info=dump.info)
    
    phis<-seq(0,1,length.out=philength)
    output$SCAM.byphiCovAll2<-cv.estimator.average(estimator=solve.synth,
                                        data=ADHdata,
                                        model='average',
                                        tune.pars=list(
                                          Kvals=Kvals, phis=phis, set.k=NA, min.preperiods=ADHdata$firstforecast,
                                          Vfun=Cov.Vars,
                                          maincovariates=allcovariates,
                                          type="all"
                                        ),
                                        est.options=list(Wbar=solve.covmatch,
                                                         sigf.ipop=5,
                                                         Margin.ipop=5e-04),
                                        only.synth=FALSE,only.match=FALSE,dump.info=dump.info,
                                        forecast.minlength=1,forecast.maxlength=1)

    
    #NOTE: this replicates the SC estimator with phi = 0 and a diagonal V matrix (b/c the SC code normalizes covariates)
    phis<-  sapply(1:length(maxphi.abadiecov[[1]]),function(x) seq(0,maxphi.abadiecov[x],length.out=philength))
    output$SCAM.abadieCov<- cv.solver(estimator=solve.covreg,
                                     data=ADHdata,
                                     model='regularize',
                                     tune.pars.joint=lapply(phis,function(x) list(
                                       K=1, phi=x, set.k=NA, min.preperiods=ADHdata$firstforecast,
                                       Vfun=Cov.Vars,
                                       maincovariates=allcovariates,
                                       type="all"
                                     )),
                                     tune.pars.list=NULL,
                                     est.options=NULL,
                                     only.synth=FALSE,only.match=FALSE,dump.info=dump.info)
    
    phis<-  sapply(1:length(maxphi.abadiecovscweight[[1]]),function(x) seq(0,maxphi.abadiecovscweight[x],length.out=philength))
    output$SCAM.abadieCovSCweight<- cv.solver(estimator=solve.covreg,
                                      data=ADHdata,
                                      model='regularize',
                                      tune.pars.joint=lapply(phis,function(x) list(
                                        K=1, phi=x, set.k=NA, min.preperiods=ADHdata$firstforecast,
                                        scV=  unlist(output$ADH$tune.pars$V),
                                        Vfun=Cov.Vars,
                                        maincovariates=allcovariates,
                                        type="all"
                                      )),
                                      tune.pars.list=NULL,
                                      # est.options=list(sigf.ipop=5,
                                      #                  Margin.ipop=5e-04,
                                      #                  Wbar=Wbar
                                      #                  solver="ipop"
                                      #                  )
                                      #,
                                      only.synth=FALSE,only.match=FALSE,dump.info=dump.info)
    
    

        
    output$SCAM.setCV5<-cv.solver(estimator=solve.synth,
                                     data=ADHdata,
                                     model="analytic",
                                     tune.pars.joint=lapply(Kvals,function(x) list(
                                       K=x, phi=0, set.k=5, min.preperiods=NA,
                                       Vfun=Cov.Vars,
                                       maincovariates=allcovariates,
                                       type="all"
                                     )),
                                     tune.pars.list=NULL,
                                     est.options=list(Wbar=solve.covmatch,
                                                      sigf.ipop=5,
                                                      Margin.ipop=5e-04),
                                     only.synth=FALSE,only.match=FALSE,
                                     dump.info=dump.info,forecast.maxlength=5, forecast.minlength=1)
    
    output$SCAM.setCV3<-cv.solver(estimator=solve.synth,
                                  data=ADHdata,
                                  model="analytic",
                                  tune.pars.joint=lapply(Kvals,function(x) list(
                                    K=x, phi=0, set.k=5:7, min.preperiods=NA,
                                    Vfun=Cov.Vars,
                                    maincovariates=allcovariates,
                                    type="all"
                                  )),
                                  tune.pars.list=NULL,
                                  est.options=list(Wbar=solve.covmatch,
                                                   sigf.ipop=5,
                                                   Margin.ipop=5e-04),
                                  only.synth=FALSE,only.match=FALSE,
                                  dump.info=dump.info,forecast.maxlength=3, forecast.minlength=1)
    


    output$SCAM.multiCV5<-cv.solver(estimator=solve.synth,
                                    data=ADHdata,
                                    model="analytic",
                                    tune.pars.joint=lapply(Kvals,function(x) list(
                                      K=x, phi=0, set.k=NA, min.preperiods=ADHdata$firstforecast,
                                      Vfun=Cov.Vars,
                                      maincovariates=allcovariates,
                                      type="all"
                                    )),
                                    tune.pars.list=NULL,
                                    est.options=list(Wbar=solve.covmatch,
                                                     sigf.ipop=5,
                                                     Margin.ipop=5e-04),
                                    only.synth=FALSE,only.match=FALSE,
                                    dump.info=dump.info,forecast.maxlength=5)
    
    
    output$SCAM.multiCV3<-cv.solver(estimator=solve.synth,
                                    data=ADHdata,
                                    model="analytic",
                                    tune.pars.joint=lapply(Kvals,function(x) list(
                                      K=x, phi=0, set.k=NA, min.preperiods=ADHdata$firstforecast,
                                      Vfun=Cov.Vars,
                                      maincovariates=allcovariates,
                                      type="all"
                                    )),
                                    tune.pars.list=NULL,
                                    est.options=list(Wbar=solve.covmatch,
                                                     sigf.ipop=5,
                                                     Margin.ipop=5e-04),
                                    only.synth=FALSE,only.match=FALSE,
                                    dump.info=dump.info,forecast.maxlength=3)
    
    

  output$name<-data$name
  
  return(output)
}

setwd(basepath)
saveRDS(params,"data_parameters.rds")

names(params[[1]]$ADH)

stopCluster(cl)
mpi.quit()

