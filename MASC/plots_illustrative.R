library(reshape2)
library(data.table)
library(Synth)
library(gurobi)
library(plyr)
library(RColorBrewer)
library(scales)
library(forecast)
library(snow)
library(MASS)
library(ggpubr)

setwd(codepath)
source("Estimator_Code.R")
source("EmpiricalARoutcome_Code.R")
source("Cross-validation_Code_byK.R")
source("plot_rules.R")



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
lines <- list(16)



philength<-100
omits<-list()
firstforecast<-c(8)

nameappend <- c("Sp_")


basedata<-data
basedata[[1]]$outcome<-basedata[[1]]$outcome[,-1]
basedata[[1]]$covariates<-basedata[[1]]$covariates[unit!=1,]
exnames<-paste0("Sp_",names[[1]][-1])

setwd(basepath)


##############################################
setwd(outputpath)


estdata<-list()
for(r in 1:length(basedata)){
  estdata[[r]] <- basedata[[r]]
  estdata[[r]]$names<-exnames
  estdata[[r]]$treatment=treatment_times[[r]]
  estdata[[r]]$firstforecast=firstforecast[r]
  estdata[[r]]$ADH_fittimes<-ADH_fittimes[[r]]
}




regularizers.dist=TRUE
regularizers.weight=TRUE


addtoplot<-list(scale_color_manual(values=c("Controls"="gray","MASC"="gray","Penalized SC"="gray","Both MASC and Pen. SC"="gray","Neither"="gray","Treated"="black")),
                scale_size_manual(values=c("Treated"=1,"Controls"=0.6,"MASC"=0.6,"Penalized SC"=0.6,"Synthetic Control"=0.6,"Both MASC and Pen. SC"=0.6,"Neither"=0.6,"Matching"=0.6)),
                scale_shape_manual(breaks=c("Synthetic Control","MASC","Matching","Penalized SC"),values=c("MASC"=21,"Matching"=22,"Penalized SC"=23,"Synthetic Control"=24,"Controls"=NA,"Treated"=NA)),
                scale_fill_manual(breaks=c("Synthetic Control","MASC","Matching","Penalized SC"),values=c("MASC"="black","Matching"="white","Penalized SC"="white","Synthetic Control"="white","Controls"=NA,"Treated"=NA)),
                labs(color="",x="Time",y="Outcome",linetype="",shape="",size="", fill=""))

seeds<-list(c( 10407L, -2134616765L, -1170759606L, -1017897887L,  1556093580L, -1109073681L,  -325238425L),
            c(10407L,   291662443L, -1982957920L,   366965598L, -2107625304L,   324778219L, -1672914269L),
           c(10407L, -1959481905L,  -946753443L, -1731197748L,   902096643L,  1975717137L,  -704206048L))


statenums<-c(17,16,17)

upperte<-c(700,500,200)
lowerte<-c(-400,-500,-200)

rescale <- c(1000)

pregdp<-data[[1]]$outcome[treatment_times[[1]]-1,1]*rescale
secaxisytreat<-scale_y_continuous(breaks=pretty_breaks(n=5), sec.axis = sec_axis(~.*100/pregdp, name = "Difference from Placebo Unit \n (Pct of GDP per Capita, Basque 1969)", breaks=pretty_breaks(n=8)))
secaxisycfact<-scale_y_continuous(breaks=pretty_breaks(n=5), sec.axis = sec_axis(~.*100/pregdp, name = "GDP per Capita Relative to Placebo Unit\n (Pct of GDP per Capita, Basque 1969)", breaks=pretty_breaks(n=8)))


#Figure 1-3: Plot SC, match, MASC
#Figure 4: Plot penalized SC as well


pdf("illustrative_plots.pdf")
s<-1
lambda<-3.3
burn_in<-100

simdata<-copy(estdata[[1]])
simdata$lambda<-lambda

for(statepos in 1:length(statenums)){
  statenum<-statenums[statepos]
  .Random.seed<-seeds[[statepos]]
  Kvals <- c(1:10)
  drawdata<-empirical_MC(s,simdata, philength=100,
                       dump.info=TRUE, Kvals=Kvals,
                       firstforecast=simdata$firstforecast,includeoutput=TRUE)
  draw<-draw.ar(1,43,fittedval=drawdata$fittedval,arcoef=drawdata$arcoef,varcov.ar=drawdata$varcov,
                simdata$covariates)
  
  iter<-statenum-1
  
  drawdata<-copy(simdata)
  drawdata$treated<-as.matrix(draw$outcome[,iter,with=FALSE])
  drawdata$donors<-as.matrix(draw$outcome[,-iter,with=FALSE])
  drawdata$outcome<-as.matrix(cbind(drawdata$treated,drawdata$donors))
  drawdata$covariates<-copy(draw$covariates)
  drawdata$covariates[,unit:=unit+1]
  drawdata$covariates[unit==iter+1,unit:=1]
  
  allcovariates<-copy(drawdata$covariates)
  allcovariates<-allcovariates[time%in%simdata$ADH_fittimes,]
  allcovariates[,time:=time-min(time)+1]
  
  ADHdata<-copy(drawdata)
  ADHdata$outcome<-ADHdata$outcome[-(1:(min(ADHdata$ADH_fittimes)-1)),]
  ADHdata$donors<-ADHdata$donors[-(1:(min(ADHdata$ADH_fittimes)-1)),]
  ADHdata$treated<-as.matrix(ADHdata$treated[-(1:(min(ADHdata$ADH_fittimes)-1)),])
  ADHdata$treatment<-length(ADHdata$ADH_fittimes)+1
  ADHdata$covariates<-ADHdata$covariates[time%in%ADHdata$ADH_fittimes,]
  ADHdata$covariates[,time:=time-min(time)+1]
  ADHdata$firstforecast<-5
  
  ###SOLVING ESTIMATORS:
  outtable<-data.table()
  output<-list(SCAM.abadie=list(),
               SCAM.average=list(), reg=list(), match=list())
  
  phis<-  lapply(0,function(x) 1)
  output$reg<-cv.solver(estimator=solve.locreg,data=drawdata,model='regularize',
                                tune.pars.list=list(phi=phis,
                                                    K=1,
                                                    set.k=NA,
                                                    min.preperiods=firstforecast), 
                                only.synth=TRUE,
                                only.match=FALSE,dump.info=TRUE)

  outtable<-rbind(outtable,data.table(Placebo=simdata$names[iter],
                                      lambda=simdata$lambda,
                                      sim=r,
                                      Model="reg",
                                      MSPE4=mean(output$reg$pred.error[1:4]^2),
                                      fit=mean(output$reg$fit^2),
                                      effect4=mean(output$reg$pred.error[1:4]),
                                      abseffect4=mean(abs(output$reg$pred.error[1:4])),
                                      cv.error=output$reg$cv.error,
                                      numweights=sum(output$reg$weights$pars>0.01),
                                      Phi=output$reg$tune.pars$phi,
                                      K=output$reg$tune.pars$K,
                                      weights=t(output$reg$weights$pars),
                                      pred.error=t(output$reg$pred.error)
  )
  ,fill=TRUE)
  
  
  phis<-  lapply(1:length(Kvals),function(x) 1)
  output$match <- cv.solver(estimator=solve.locreg,data=drawdata,
                            model='regularize',
                            tune.pars.list=list(phi=phis,K=Kvals,
                                                set.k=NA,
                                                min.preperiods=firstforecast),
                            only.synth=FALSE,only.match=TRUE,dump.info=TRUE)
  
  outtable<-rbind(outtable,data.table(Placebo=simdata$names[iter],
                                      lambda=simdata$lambda,
                                      sim=r,
                                      Model="match",
                                      MSPE4=mean(output$match$pred.error[1:4]^2),
                                      fit=mean(output$match$fit^2),
                                      effect4=mean(output$match$pred.error[1:4]),
                                      abseffect4=mean(abs(output$match$pred.error[1:4])),
                                      cv.error=output$match$cv.error,
                                      numweights=sum(output$match$weights$pars>0.01),
                                      Phi=output$match$tune.pars$phi,
                                      K=output$match$tune.pars$K,
                                      weights=t(output$match$weights$pars),
                                      pred.error=t(output$match$pred.error)
  ),fill=TRUE)
  
  
  output$SCAM.average <- cv.solver(estimator=solve.locreg,data=drawdata,
                                   tune.pars.list=list(phi=phis,K=Kvals,
                                                       set.k=NA,
                                                       min.preperiods=firstforecast),
                                   tune.pars.joint=NULL,
                                   only.synth=FALSE,only.match=FALSE,
                                   dump.info=TRUE,model='analytic',fold.errors=TRUE)
  outtable<-rbind(outtable,data.table(Placebo=simdata$names[iter],
                                      lambda=simdata$lambda,
                                      sim=r,
                                      Model="MASC",
                                      MSPE4=mean(output$SCAM.average$pred.error[1:4]^2),
                                      fit=mean(output$SCAM.average$fit^2),
                                      effect4=mean(output$SCAM.average$pred.error[1:4]),
                                      abseffect4=mean(abs(output$SCAM.average$pred.error[1:4])),
                                      cv.error=output$SCAM.average$cv.error,
                                      numweights=sum(output$SCAM.average$weights$pars>0.01),
                                      Phi=output$SCAM.average$tune.pars$phi,
                                      K=output$SCAM.average$tune.pars$K,
                                      weights=t(output$SCAM.average$weights$pars),
                                      pred.error=t(output$SCAM.average$pred.error)
  ),fill=TRUE)
  
  params <- list(SCAM.average =unlist(outtable[Placebo==paste0("Sp_",names[[1]][statenum]) 
                                                      & Model %in% c("MASC"),grepl("weights.",names(outtable)) ,with=FALSE]),
                 reg =unlist(outtable[Placebo==paste0("Sp_",names[[1]][statenum]) 
                                             & Model %in% c("reg"),grepl("weights.",names(outtable)) ,with=FALSE]),
                 match =unlist(outtable[Placebo==paste0("Sp_",names[[1]][statenum]) 
                                               & Model %in% c("match"),grepl("weights.",names(outtable)) ,with=FALSE]))

  
  plotdata<-data.table(t=1:43 + 1954,drawdata$outcome - drawdata$outcome[,1])
  names(plotdata)<-c("t",names[[1]][statenum],names[[1]][-c(1,statenum)])
  plotdata<-melt(plotdata,id.vars="t")
  plotdata<-plotdata[t<=19+1954,]
  
  plotnames<-names[[1]][-c(1,statenum)]
  
  lines<-list()
  
  if(statepos<=2){
    myplot<-(ggplot()
             +geom_line(data=plotdata[t%%2==1& variable%in%plotnames[round(params$reg,2)>0.15],],
                        aes(x=t,y=value*1000,group=variable, size = "Synthetic Control"),linetype="dashed")
             +geom_line(data=plotdata[t%%2==1& variable%in%plotnames[round(params$match,2)>0.15],],
                        aes(x=t,y=value*1000,group=variable, size = "Matching"),linetype="dashed")
             + geom_line(data=plotdata[t%%2==1& variable%in%plotnames[round(params$SCAM.average,2)>0.15],],
                         aes(x=t,y=value*1000,group=variable, size = "MASC"),linetype="dashed")
             +geom_line(data=plotdata[t%%2==1& variable%in%names[[1]][statenum],],
                        aes(x=t,y=value*1000,group=variable, size="Treated"))
             +lines
             +geom_point(data=plotdata[t%%2==1& variable%in%names[[1]][statenum],],
                         aes(x=t,y=value*1000,group=variable, shape="Treated"), size=5)
             +geom_point(data=plotdata[t%%2==1& variable%in%plotnames[round(params$reg,2)>0.15],],
                         aes(x=t,y=value*1000,group=variable, shape = "Synthetic Control", fill="Synthetic Control"), size=5)
             +geom_point(data=plotdata[t%%2==1& variable%in%plotnames[round(params$match,2)>0.15],],
                         aes(x=t,y=value*1000,group=variable, shape = "Matching", fill="Matching"), size=5)
             +geom_point(data=plotdata[t%%2==1& variable%in%plotnames[round(params$SCAM.average,2)>0.15],],
                         aes(x=t,y=value*1000,group=variable, shape = "MASC", fill="MASC"), size=3)
             +geom_vline(xintercept=1954+16,linetype="dashed")
             +addtoplot
             +secaxisycfact
             +labs(y="GDP per Capita Relative to \n Placebo Unit (1986 USD)", x="Year")
    )  
  }
    else{
  myplot<-(ggplot()
           +geom_line(data=plotdata[t%%2==1& variable%in%plotnames[round(params$reg,2)>0.15],],
                      aes(x=t,y=value*1000,group=variable, size = "Synthetic Control"),linetype="dashed")
           +geom_line(data=plotdata[t%%2==1& variable%in%plotnames[round(params$match,2)>0.15],],
                      aes(x=t,y=value*1000,group=variable, size = "Matching"),linetype="dashed")
           + geom_line(data=plotdata[t%%2==1& variable%in%plotnames[round(params$match,2)|round(params$reg,2)>0.15],],
                       aes(x=t,y=value*1000,group=variable, size = "MASC"),linetype="dashed")
           +geom_line(data=plotdata[t%%2==1& variable%in%names[[1]][statenum],],
                      aes(x=t,y=value*1000,group=variable, size="Treated"))
           +lines
           +geom_point(data=plotdata[t%%2==1& variable%in%names[[1]][statenum],],
                       aes(x=t,y=value*1000,group=variable, shape="Treated"), size=5)
           +geom_point(data=plotdata[t%%2==1& variable%in%plotnames[round(params$reg,2)>0.15],],
                       aes(x=t,y=value*1000,group=variable, shape = "Synthetic Control", fill="Synthetic Control"), size=5)
           +geom_point(data=plotdata[t%%2==1& variable%in%plotnames[round(params$match,2)>0.15],],
                       aes(x=t,y=value*1000,group=variable, shape = "Matching", fill="Matching"), size=5)
           +geom_point(data=plotdata[t%%2==1& variable%in%plotnames[round(params$match,2)|round(params$reg,2)>0.15],],
                       aes(x=t,y=value*1000,group=variable, shape = "MASC", fill="MASC"), size=3)
           +geom_vline(xintercept=1954+16,linetype="dashed")
           +addtoplot
           +secaxisycfact
           +labs(y="GDP per Capita Relative to \n Placebo Unit (1986 USD)", x="Year")
  )
    }
  print(myplot+   theme(legend.position="none") + guides(fill=FALSE))
  
  
  
  plotestimates<-data.table(t=1:43 + 1954, 
                            MASC = drawdata$outcome[,-(1)]%*%params$SCAM.average - drawdata$outcome[,(1)],
                            SC = drawdata$outcome[,-(1)]%*%params$reg - drawdata$outcome[,(1)],
                            Matching = drawdata$outcome[,-(1)]%*%params$match - drawdata$outcome[,(1)],
                            Treated =  drawdata$outcome[,(1)] - drawdata$outcome[,(1)]
  )
  
  names(plotestimates)<-c("t","MASC","SC","Matching","Treated")
  
  plotestimates<-plotestimates[t<=19+1954,]
  
  myplot<-(ggplot()
           + geom_line(data = plotestimates, aes(x=t,y=SC*1000))
           + geom_line(data = plotestimates, aes(x=t,y=Matching*1000))
           +geom_line(data = plotestimates, aes(x=t,y=MASC*1000))
           + geom_line(data = plotestimates, aes(x=t,y=Treated*1000),size=1)
           + geom_point(data = plotestimates[t%%2==1,], aes(x=t,y=SC*1000,shape="Synthetic Control", fill="Synthetic Control"),size=5)
           + geom_point(data = plotestimates[t%%2==1,], aes(x=t,y=Matching*1000,shape="Matching", fill="Matching"),size=5)
           + geom_point(data = plotestimates[t%%2==1,], aes(x=t,y=Treated*1000,shape="Treated"),size=5)
           + geom_point(data = plotestimates[t%%2==1,], aes(x=t,y=MASC*1000,shape="MASC", fill ="MASC"),size=3)
           + geom_vline(xintercept=1954+16,linetype="dashed")
           +addtoplot
           +secaxisytreat
           +labs(y="Difference from Placebo Unit \n (GDP Per Capita, 1986 USD)", x="Year")
  )
  print(myplot+   theme(legend.position="none"))
  print(as_ggplot(get_legend(myplot+theme(legend.direction="horizontal") )))
  
  
  
  
  
  
}
dev.off()
