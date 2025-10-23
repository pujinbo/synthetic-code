
####PRELIMINARIES#######
library(reshape2)
library(data.table)
library(Synth)
library(gurobi)
library(plyr)
library(forecast)
library(snow)
library(MASS)
library(doSNOW)
library(doRNG)
library(Rmpi)#COMMENT OUT if running locally
library(optimx)

basepath<-"/home/mdkellogg/Synthetic_Controls/replication_code/version4"
codepath<-paste0(basepath,"/auxfiles")
inputpath<- basepath
outputpath<- paste0(basepath,"/figures")

setwd(codepath)
source("Estimator_Code.R")
source("EmpiricalARoutcome_Code.R")
source("Cross-validation_Code_byK.R")

cl <- makeCluster(mpi.universe.size()-1, type="MPI",outfile="") #COMMENT OUT if running locally
# cl <- makeCluster(4, type="SOCK")  #local windows parallelization
registerDoSNOW(cl)
clusterEvalQ(cl,{
  library(reshape2)
  library(data.table)
  library(Synth)
  library(xtable)
  library(gurobi)
  library(plyr)
  library(RColorBrewer)
  library(scales)
  library(snow)
  library(tseries)
  library(forecast)
  library(MASS)
  library(Matrix)
  library(doSNOW)
  library(geometry)
  library(kernlab)
  library(optimx)
})



#####################
#####SETTING UP DATA##

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

firstforecast<-c(8)

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
treatment_times<-list()
ADH_fittimes<-list()
treatment_times[[1]]<-16
ADH_fittimes[[1]]<-6:15
lines <- list(15.5)



omits<-list()
nameappend <- c("Sp_")

basedata<-data
basedata[[1]]$outcome<-basedata[[1]]$outcome[,-1]
basedata[[1]]$covariates<-basedata[[1]]$covariates[unit!=1,]
exnames<-paste0("Sp_",names[[1]][-1])
##############################################


Kvals<-1:10
emp.results<-list()
phivals<-(0:100)/100

estdata<-list()
for(r in 1:length(basedata)){
  estdata[[r]] <- basedata[[r]]
  estdata[[r]]$names<-exnames
  estdata[[r]]$treatment=treatment_times[[r]]
  estdata[[r]]$firstforecast=firstforecast[r]
  estdata[[r]]$ADH_fittimes<-ADH_fittimes[[r]]
}

burn_in=100

###DRAW AND SAVE MULTIVARIATE ERRORS:
#this may take some time...
errors<-NULL
lambdas<-c(3.0,  3.3)
# #lambda=3.0: interactive FE model, no noise
# #lambda=3.3: interactive FE model, white noise process 



r<-1
if(r==1)  file.remove("temp_empiricalMC.csv")
for(lambda in lambdas){
  set.seed(23921)
  errorlambda<-round((lambda%%1)*10,1)
  if(errorlambda==0) draws<-1
  if(errorlambda!=0) draws<-1000
  runsims<-r:(r+draws-1)
  print(paste(r,r+draws-1))
  
  simdata<-estdata[[1]]
  simdata$lambda<-lambda
  
  ex<- Filter(function(x) is.function(get(x,.GlobalEnv)),ls(.GlobalEnv))
  clusterExport(cl,ex)
  if(length(runsims)>1){
    emp.results<-foreach(s = runsims,.combine='rbind') %dorng%{
      setwd(basepath)
      print(paste0("Run ",s," of ",max(runsims)))
      gc()
        result<-try(empirical_MC(s,simdata, philength=100,
                           dump.info=TRUE, Kvals=Kvals,
                           firstforecast=simdata$firstforecast,fold.errors=TRUE))
        while(inherits(result,"try-error")) {
          #change seed value, run a draw of something:
          result<-try(empirical_MC(s,simdata, philength=100,
                                   dump.info=TRUE, Kvals=Kvals,
                                   firstforecast=simdata$firstforecast,fold.errors=TRUE))
        }
        write.table(result, file="temp_empiricalMC.csv",append=TRUE,col.names=FALSE,sep=",")
      return(result)
      
    }
  } else{
    s<-runsims
    setwd(basepath)
    print(paste0("Run ",s," of ",length(runsims)))
    gc()
    emp.results<-empirical_MC(s,simdata, philength=100,
                              dump.info=TRUE, Kvals=Kvals,
                              firstforecast=simdata$firstforecast,fold.errors=TRUE)
  }
  r<-r+draws
  setwd(basepath)
  #Intermediate file:
  saveRDS(emp.results,paste0('empiricalMCadjust_distdatafull_lambda',lambda*100,'.RDS'))
  emp.results<-as.data.table(emp.results)
  #Transform MSPE and MSE into RMSPE and RMSE:
  emp.results[,MSPE4:=MSPE4^0.5]
  emp.results[,MSPEAll:=MSPEAll^0.5]
  emp.results[,fit:=fit^0.5]
  
  emp.results[,meanMSPE4:=mean(MSPE4),by=.(Model,Placebo,lambda,Phi,K)]
  emp.results[Model%in%c("MASC.byphi","Penalized SC.byphi","MASC.byphi.outcome","Penalized SC.byphi.outcome","Penalized SC (SC weights).byphi"),minmeanMSPE4:=min(meanMSPE4),by=.(Model,Placebo,lambda)]
  emp.results[Model%in%c("MASC.byphi","Penalized SC.byphi","MASC.byphi.outcome","Penalized SC.byphi.outcome","Penalized SC (SC weights).byphi"),minPhi4:=min((meanMSPE4==minmeanMSPE4)*Phi+9999*(meanMSPE4!=minmeanMSPE4)),by=.(Model,Placebo,lambda)]
  emp.results[Model%in%c("MASC.byphi","Penalized SC.byphi","MASC.byphi.outcome","Penalized SC.byphi.outcome","Penalized SC (SC weights).byphi"),minK4:=min((meanMSPE4==minmeanMSPE4)*K+9999*(meanMSPE4!=minmeanMSPE4)),by=.(Model,Placebo,lambda)]
  emp.results[Model%in%c("MASC.byphi","MASC"),minMSPE4:=min(MSPE4*(Phi==minPhi4&K==minK4) + 9999*(Phi!=minPhi4|K!=minK4),na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("MASC.byphi.outcome","MASC.outcome"),minMSPE4:=min(MSPE4*(Phi==minPhi4&K==minK4) + 9999*(Phi!=minPhi4|K!=minK4),na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("Penalized SC.byphi","Penalized SC"),minMSPE4:=min(MSPE4*(Phi==minPhi4&K==minK4) + 9999*(Phi!=minPhi4|K!=minK4),na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("Penalized SC.byphi.outcome","Penalized SC.outcome"),minMSPE4:=min(MSPE4*(Phi==minPhi4&K==minK4) + 9999*(Phi!=minPhi4|K!=minK4),na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("Penalized SC (SC weights).byphi","Penalized SC (SC weights)"),minMSPE4:=min(MSPE4*(Phi==minPhi4&K==minK4) + 9999*(Phi!=minPhi4|K!=minK4),na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("MASC.byphi","MASC"),minPhi2:=min(minPhi4,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("MASC.byphi","MASC"),minK2:=min(minK4,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("MASC.byphi.outcome","MASC.outcome"),minPhi2:=min(minPhi4,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("MASC.byphi.outcome","MASC.outcome"),minK2:=min(minK4,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("Penalized SC.byphi","Penalized SC"),minPhi2:=min(minPhi4,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("Penalized SC.byphi","Penalized SC"),minK2:=min(minK4,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("Penalized SC.byphi.outcome","Penalized SC.outcome"),minPhi2:=min(minPhi4,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("Penalized SC.byphi.outcome","Penalized SC.outcome"),minK2:=min(minK4,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("Penalized SC (SC weights).byphi","Penalized SC (SC weights)"),minPhi2:=min(minPhi4,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("Penalized SC (SC weights).byphi","Penalized SC (SC weights)"),minK2:=min(minK4,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  
  emp.results[,minK4:=minK2]
  emp.results[,minPhi4:=minPhi2]
  emp.results[,minPhi2:=NULL]
  emp.results[,minK2:=NULL]
  
  emp.results[,meanMSPEAll:=mean(MSPEAll),by=.(Model,Placebo,lambda,Phi,K)]
  emp.results[Model%in%c("MASC.byphi","Penalized SC.byphi","MASC.byphi.outcome","Penalized SC.byphi.outcome","Penalized SC (SC weights).byphi"),minmeanMSPEAll:=min(meanMSPEAll),by=.(Model,Placebo,lambda)]
  emp.results[Model%in%c("MASC.byphi","Penalized SC.byphi","MASC.byphi.outcome","Penalized SC.byphi.outcome","Penalized SC (SC weights).byphi"),minPhi:=min((meanMSPEAll==minmeanMSPEAll)*Phi+9999*(meanMSPEAll!=minmeanMSPEAll)),by=.(Model,Placebo,lambda)]
  emp.results[Model%in%c("MASC.byphi","Penalized SC.byphi","MASC.byphi.outcome","Penalized SC.byphi.outcome","Penalized SC (SC weights).byphi"),minK:=min((meanMSPEAll==minmeanMSPEAll)*K+9999*(meanMSPEAll!=minmeanMSPEAll)),by=.(Model,Placebo,lambda)]
  emp.results[Model%in%c("MASC.byphi","MASC"),minMSPEAll:=min(MSPEAll*(Phi==minPhi&K==minK) + 9999*(Phi!=minPhi|K!=minK),na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("MASC.byphi.outcome","MASC.outcome"),minMSPEAll:=min(MSPEAll*(Phi==minPhi&K==minK) + 9999*(Phi!=minPhi|K!=minK),na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("Penalized SC.byphi","Penalized SC"),minMSPEAll:=min(MSPEAll*(Phi==minPhi&K==minK) + 9999*(Phi!=minPhi|K!=minK),na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("Penalized SC.byphi.outcome","Penalized SC.outcome"),minMSPEAll:=min(MSPEAll*(Phi==minPhi&K==minK) + 9999*(Phi!=minPhi|K!=minK),na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("Penalized SC (SC weights).byphi","Penalized SC (SC weights)"),minMSPEAll:=min(MSPEAll*(Phi==minPhi&K==minK) + 9999*(Phi!=minPhi|K!=minK),na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("MASC.byphi","MASC"),minPhi2:=min(minPhi,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("MASC.byphi","MASC"),minK2:=min(minK,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("MASC.byphi.outcome","MASC.outcome"),minPhi2:=min(minPhi,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("MASC.byphi.outcome","MASC.outcome"),minK2:=min(minK,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("Penalized SC.byphi","Penalized SC"),minPhi2:=min(minPhi,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("Penalized SC.byphi","Penalized SC"),minK2:=min(minK,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("Penalized SC.byphi.outcome","Penalized SC.outcome"),minPhi2:=min(minPhi,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("Penalized SC.byphi.outcome","Penalized SC.outcome"),minK2:=min(minK,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("Penalized SC (SC weights).byphi","Penalized SC (SC weights)"),minPhi2:=min(minPhi,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  emp.results[Model%in%c("Penalized SC (SC weights).byphi","Penalized SC (SC weights)"),minK2:=min(minK,na.rm=TRUE),by=.(Placebo,lambda,sim)]
  
  emp.results[,minK:=minK2]
  emp.results[,minPhi:=minPhi2]
  emp.results[,minPhi2:=NULL]
  emp.results[,minK2:=NULL]
  

  emp.results<-emp.results[!Model%in%c("Penalized SC.byphi","MASC.byphi","MASC.byphi.outcome","Penalized SC.byphi.outcome","Penalized SC (SC weights).byphi"),]
  emp.results[,minmeanMSPE4:=NULL]
  emp.results[,minmeanMSPEAll:=NULL]
  emp.results[,meanMSPE4:=NULL]

  emp.results[Model=="reg",Model:="SC (Covs)"]
  emp.results[Model=="match",Model:="Matching (Covs)"]
  emp.results[Model=="MASC",Model:="MASC (Covs)"]
  emp.results[Model=="reg.outcome",Model:="SC"]
  emp.results[Model=="match.outcome",Model:="Matching"]
  emp.results[Model=="MASC.outcome",Model:="MASC"]
  emp.results[Model=="Penalized SC",Model:="Pen. SC (Covs)"]
  emp.results[Model=="Penalized SC.outcome",Model:="Pen. SC"]
  emp.results[Model=="Penalized SC (SC weights)",Model:="Pen. SC (SC weights)"]
  
  emp.results[,Placebo:=gsub('Sp_','',Placebo)]
  emp.results[Placebo=="Andalucia",Placebo:="Andalusia"]
  emp.results[Placebo=="Principado",Placebo:="Asturias"]
  emp.results[Placebo=="Baleares",Placebo:="Balearic Isl"]
  emp.results[Placebo=="Canarias",Placebo:="Canary Isl"]
  emp.results[Placebo=="Cataluna",Placebo:="Catalonia"]
  emp.results[Placebo=="Castilla",Placebo:="Cast-Leon"]
  emp.results[Placebo=="Castilla-La",Placebo:="Cast-Mancha"]
  emp.results[Placebo=="Comunidad",Placebo:="Valencia"]
  emp.results[Placebo=="Navarra",Placebo:="Navarre"]
  
  
  saveRDS(emp.results,paste0('empiricalMCadjust_distdata_lambda',lambda*100,'.RDS'))
  
}

lambdas<-c(3.0, 3.3)
emp.results<-NULL
for(lambda in lambdas){
  temp<-readRDS(paste0('empiricalMCadjust_distdata_lambda',lambda*100,'.RDS'))
  emp.results<-rbind(emp.results,temp)
}
saveRDS(emp.results,'empiricalMCadjust_distdata.RDS')


stopCluster(cl)
mpi.quit()



