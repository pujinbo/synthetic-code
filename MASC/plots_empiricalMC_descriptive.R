###########PRELIMINARIES##############
library(reshape2)
library(ggplot2)
library(data.table)
library(Synth)
library(plyr)
library(RColorBrewer)
library(scales)
library(forecast)
library(MASS)
library(ggpubr)
library(ggpattern)


setwd(codepath)
source("Estimator_Code.R")
source("EmpiricalARoutcome_Code.R")
source("Cross-validation_Code_byK.R")
source("plot_rules.R")


theme_update(axis.text.y = element_text(color='black',size=22),
            axis.title.x=element_text(color="black",size=22),
            axis.title.y=element_text(color="black",size=22),
             axis.text.x = element_text(color='black',angle = 0, hjust = 0.5,size=22))


set.seed(23921)


rescale <- c(1000)
meangdpy<-0.1591586*rescale
pregdp<-6.081405 *rescale
basque_fit<-0.0755584*rescale

secaxisydiff<-scale_y_continuous(breaks=pretty_breaks(n=5), sec.axis = sec_axis(~.*100/pregdp, name = "Difference from Trend \n (Pct of GDP per Capita, Basque 1969)", breaks=pretty_breaks(n=5)))


##################SETTING UP DATA###########
setwd(inputpath)

MAXT<-43


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
lines <- list(16)



omits<-list()
nameappend <- c("Sp_")

basedata<-data
basedata[[1]]$outcome<-basedata[[1]]$outcome[,-1]
basedata[[1]]$covariates<-basedata[[1]]$covariates[unit!=1,]
exnames<-paste0("Sp_",names[[1]][-1])
##############################################



omits<-list()
nameappend <- c("Sp_")
graphnames<-c("Spain")
diffupper <- c(5,  75,    10,   20,  30, 50, 50,  100, 20,  50, 30,  75, 20, 20, 90, 60)
difflower <- c(-5,-75,  -10,  -20, -30, -50,-50,-80, -20,-100,-30,-75, -20,-20,-60,-90)


addtime<-1954

setwd(outputpath)


estdata<-list()
for(r in 1:length(basedata)){
  estdata[[r]] <- basedata[[r]]
  estdata[[r]]$names<-exnames
  estdata[[r]]$treatment=treatment_times[[r]]
  estdata[[r]]$firstforecast=firstforecast[r]
  estdata[[r]]$ADH_fittimes<-ADH_fittimes[[r]]
}

burn_in=100


#############PLOTS##############
pdf("empirical simulation plots.pdf")

setwd(inputpath)

exnames<-graphnames[grepl('Sp_',graphnames)]
lambdas<-c(3.3)
for(lambda in lambdas){
  
  simdata<-copy(estdata[[1]])
  simdata$lambda<-lambda

  #FIT FACTOR MODEL TO DATASET:
  dataparams<-empirical_MC(r=1,data=simdata,includeoutput=TRUE)
  
  
  draw<-draw.ar(1,MAXT=MAXT,fittedval=dataparams$fittedval,arcoef=dataparams$arcoef,varcov.ar=dataparams$varcov,
                simdata$covariates)
  
  draw<-data.table(draw$outcome)
  fittedtable<-as.data.table(draw)
  fittedtable[,t:=.I]
  draw[,t:=.I]
  draw<-melt(draw,id.vars='t')
  draw[,variable:=gsub('V','',variable)]
  setnames(draw,'variable','state')
  draw[,state:=as.factor(as.numeric(state))]
  
  fittedtable<-melt(fittedtable,id.vars='t')
  fittedtable[,variable:=gsub('V','',variable)]
  setnames(fittedtable,'variable','state')
  fittedtable[,state:=as.factor(as.numeric(state))]
  

  plotdata<-data.table(simdata$outcome)
  plotdata[,t:=.I]
  plotdata<-melt(plotdata,id.vars='t')
  plotdata[,variable:=gsub('V','',variable)]
  setnames(plotdata,'variable','state')
  plotdata[,state:=as.factor(as.numeric(state))]
  plotdata<-plotdata[t<=MAXT,]
  plotdata[,error:=value-dataparams$fittedval[,state],by=state]
  
  
  mindraw<-NULL
  maxdraw<-NULL
  minerror<-NULL
  maxerror<-NULL
  draws<-data.table()
  errors<-data.table()
  for(d in 1:10){
    draw<-draw.ar(1,MAXT=MAXT,fittedval=dataparams$fittedval,arcoef=dataparams$arcoef,varcov.ar=dataparams$varcov,
                  simdata$covariates)
    draw<-data.table(draw$outcome)
    error<-draw-dataparams$fittedval
    draw[,t:=.I]
    draw<-melt(draw,id.vars='t')
    draw[,variable:=gsub('V','',variable)]
    setnames(draw,'variable','state')
    draw[,state:=as.factor(as.numeric(state))]
    draw[,draw:=d]
    draws<-rbind(draws,draw)
    mindraw<-min(mindraw,draw$value)
    maxdraw<-max(maxdraw,draw$value)
    
    error[,t:=.I]
    error<-melt(error,id.vars='t')
    error[,variable:=gsub('V','',variable)]
    setnames(error,'variable','state')
    error[,state:=as.factor(as.numeric(state))]
    error[,draw:=d]
    errors<-rbind(errors,error)
    minerror<-min(minerror,error$value)
    maxerror<-max(maxerror,error$value)
    
  }
  
  for(region in 1:ncol(dataparams$fittedval)){

    draws[state==region,diff:=value-plotdata[state==region,value],by=.(draw)]
    
    
    #Plotting differences of draws of cataluna from trend:
    print(ggplot(data=errors[state==region&draw<=5,],aes(x=t+addtime,y=value*rescale,group=draw,color=draw))
          +geom_line(size=0.5,linetype='dashed')
          #        +geom_point()
          +geom_line(data=plotdata[state==region,],aes(x=t+addtime,y=error*rescale,group='real'),color='black',size=1)
          +scale_color_continuous(low='grey20',high='grey80')
#          +secaxisydiff
#          +coord_cartesian(ylim=c(difflower[region], diffupper[region]))
          +labs(y='Difference from Trend \n (GDP per Capita, 1986 USD)',x='Year')
          +geom_hline(yintercept=0,color='black',linetype='dashed')
          +guides(color='none')
    )
  }
  
  
  
}



dev.off()


