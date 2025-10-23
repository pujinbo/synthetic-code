library(reshape2)
library(ggplot2)
library(data.table)
library(Synth)
library(plyr)
library(RColorBrewer)
library(scales)
library(ggpubr)
library(ggpattern)
library(forecast)
library(forcats)

setwd(codepath)
source("Estimator_Code.R")
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


treatment_times<-list()
treatment_times[[1]]<-16
names<-list()
names[[1]]<- c(unique(basque[regionno==17,regionname]),unique(basque[regionno!=17,regionname]))

data<-list()
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


names[[1]][2]<-"Andalusia"
names[[1]][5]<-"Balearic Isl"
names[[1]][6]<-"Canary Isl"
names[[1]][9]<-"Cast-Mancha"
names[[1]][8]<-"Cast-Leon"
names[[1]][10]<-"Catalonia"
names[[1]][11]<-"Valencia"
names[[1]][16]<-"Navarre"
names[[1]][4]<-"Asturias"


theme_update(axis.text.y = element_text(color='black',size=14),
             axis.text.x = element_text(color='black',angle = 0, hjust = 0.5,size=14),
             strip.text.x=element_text(size=14))


set.seed(23921)

rescale <- c(1000)

pregdp<-6.081405 *rescale

secaxisy<-scale_y_continuous(breaks=pretty_breaks(n=7), sec.axis = sec_axis(~.*100/pregdp, name = "Expected RMSPE \n (Pct of GDP per Capita, Basque 1969)", breaks=pretty_breaks(n=5)))
secaxisydist<-scale_y_continuous(breaks=pretty_breaks(n=7), sec.axis = sec_axis(~.*100/pregdp, name = "RMSPE \n (Pct of GDP per Capita, Basque 1969)", breaks=pretty_breaks(n=4)))
secaxisyscfit<-scale_y_continuous(breaks=pretty_breaks(n=7), sec.axis = sec_axis(~.*100/pregdp, name = "Exected Error \n (Pct of GDP per Capita, Basque 1969)", breaks=pretty_breaks(n=5)))

drop_placebos<-c("Balearic Isl","Extremadura","Madrid")
high_placebo_list<-list(
  c("Andalusia","Asturias","Canary Isl","Cantabria","Cast-Mancha","Catalonia","Murcia","Navarre","Rioja"),
  c("Andalusia","Asturias","Canary Isl","Cantabria","Cast-Mancha","Catalonia","Murcia","Navarre","Rioja"),
  c("Andalusia","Asturias","Canary Isl","Cantabria","Cast-Mancha","Catalonia","Murcia","Navarre","Rioja")
)


plot_Emp_MC<-function(performance){
  
  placebo.order<-list()
  for(dt in sort(unique(performance$lambda))){
    errorlambda<-round((dt%%1)*10,1)
    print(paste0("Lambda Value: ",dt))
    dtpos <- which(dt==levels(as.factor(performance$lambda)))
    high_placebos<-high_placebo_list[[dtpos]]
    meanperformance<-performance[Model%in%c("MASC",
                                            "Matching","SC","PSC-AG","PSC-S")&lambda==dt,.(MSPE4=mean(MSPE4),
                                                                                    MSPEAll=mean(MSPEAll),
                                                                                    fit=mean(fit),
                                                                                    Phi=mean(Phi), K =mean(K),
                                                                                    numweights=mean(numweights),
                                                                                    minMSPE4=mean(minMSPE4),
                                                                                    minMSPEAll=mean(minMSPEAll)),by=.(Placebo,Model)]
    
    meltperformance<-as.data.table(melt(meanperformance,id.vars=c("Placebo","Model")))
    meltperformance[,value:=as.numeric(value)]

    if(dt==3.3){
    #Optimal vs. average selected MSPE, comparing over regions:
    myplot<-(ggplot(data=meltperformance[Model%in%c("MASC")& !Placebo%in%drop_placebos & variable %in%c("MSPEAll","minMSPEAll"),]) +
               geom_col(aes(y=value,x=Placebo,group=factor(variable,levels=c("MSPEAll","minMSPEAll")),fill=variable),position="dodge") +
               secaxisy+
               scale_fill_manual(breaks=c("MSPEAll","minMSPEAll"),values=gray.colors(n=2,start=0.3,end=0.8),
                                 labels=c("Cross-Validated", "Minimal"), name=element_blank())+
               labs(y=expression(paste("Expected RMSPE (GDP per capita)")),x="")
             +textadjust)
    print(myplot+   theme(legend.position="none"))
    print(as_ggplot(get_legend(myplot+guides(color=guide_legend(order=1))+theme(legend.key.size=unit(2,"lines"),legend.direction="horizontal"))))
    


    #Fit vs. performance of SC:
    myplot<-(ggplot(data=meltperformance[Model%in%c("SC")& !Placebo%in%drop_placebos & variable %in%c("MSPEAll","fit"),]) +
               geom_col(aes(y=value,x=Placebo,group=factor(variable,levels=c("MSPEAll","fit")),fill=variable),position="dodge") +
               secaxisyscfit+
               scale_fill_manual(breaks=c("MSPEAll","fit"),values=gray.colors(n=2,start=0.3,end=0.8),
                                 labels=c("RMSPE", "Pre-period fit (RMSE)"), name=element_blank())+
               labs(y=expression(paste("Expected error (GDP per capita)")),x="")
             +textadjust)
    print(myplot+   theme(legend.position="none"))
    print(as_ggplot(get_legend(myplot+guides(color=guide_legend(order=1))+theme(legend.key.size=unit(2,"lines"),legend.direction="horizontal"))))
    

    }
    performance.plot<-copy(meanperformance)
    performance.plot[Model=="SC",Model:="Synth. Control"]

    
    performance.plot[,Model:=factor(Model,levels=c('MASC','Matching','Synth. Control',"PSC-AG","PSC-S"))]
    
    #Stats on MSPE for MASC, in each placebo:
    yupper <-c(1.1, 1, 0.9)*rescale

    if(dt==4.0) {
      performance.plot<-rbind(performance.plot,
                              data.table(Placebo="",Model=c("Synth. Control","Matching","MASC","Penalized SC"),
                                         MSPEAll=0), fill=TRUE)
      performance.plot[,Placebo:=fct_relevel(Placebo,c("",high_placebos),after=Inf)]
      myplot<-(ggplot(data=performance.plot[Model%in%c("MASC","Synth. Control","Penalized SC", "Matching")
                                            & !Placebo%in%drop_placebos
                                            ,]) +
                 geom_col(aes(x=Placebo,
                              y=as.numeric(MSPEAll),
                              fill=factor(Model,levels=c("MASC","Synth. Control","Penalized SC", "Matching")),
                              group = factor(Model,levels=c("MASC","Synth. Control","Penalized SC", "Matching"))),position="dodge",
                          color="black") +
                 geom_vline(xintercept=5,linetype="dashed")+
                 list(scale_fill_manual(breaks=c("MASC","Synth. Control","Penalized SC","Matching"),
                                        values=c('Synth. Control'=colcols[2], 'Matching'="white",
                                                 'MASC'=colcols[1],
                                                 'Penalized SC'=colcols[3]),
                                        name = element_blank()
                 ),
                 scale_color_manual(breaks=c("MASC","Synth. Control","PSC-AG","PSC-S","Matching"),
                                    values=c('Synth. Control'=colcols[2], 'Matching'="white",'MASC'=colcols[1],
                                             'PSC-S'=colcols[4],'PSC-AG'=colcols[3]),
                                    name = element_blank()
                 )
                 )+
                 coord_cartesian(ylim=c(0,yupper[dtpos]))+
                 labs(x='Placebo',y="Expected RMSPE  (GDP per capita)")+
                 secaxisy+
                 textadjust)
      print(myplot+   theme(legend.position="none"))
      print(as_ggplot(get_legend(myplot+guides(fill=guide_legend(nrow=2))+theme(legend.key.size=unit(2,"lines"),legend.direction="horizontal"))))
    }
    if(dt==3.0){
      myplot<-(ggplot(data=performance.plot[Model%in%c("MASC","Synth. Control","PSC-S","PSC-AG", "Matching")
                                            & !Placebo%in%drop_placebos
                                            ,]) +
                 geom_col(aes(x=Placebo,y=as.numeric(MSPEAll), fill=factor(Model,levels=c("MASC","Synth. Control","PSC-S","PSC-AG", "Matching")),
                              group = factor(Model,levels=c("MASC","Synth. Control","PSC-AG","PSC-S", "Matching"))),position="dodge",
                          color="black") +
                 list(scale_fill_manual(breaks=c("MASC","Synth. Control","PSC-AG","PSC-S","Matching"),
                                        values=c('Synth. Control'=colcols[2], 'Matching'="white",
                                                 'MASC'=colcols[1], 'PSC-S'=colcols[4],'PSC-AG'=colcols[3],
                                                 'Penalized SC'=colcols[3]),
                                        name = element_blank()
                 ),
                 scale_color_manual(breaks=c("MASC","Synth. Control","PSC-AG","PSC-S","Matching"),
                                    values=c('Synth. Control'=colcols[2], 'Matching'="white",'MASC'=colcols[1],
                                             'PSC-S'=colcols[4],'PSC-AG'=colcols[3],
                                             'Penalized SC'=colcols[3]),
                                    name = element_blank()
                 )
                 )+
                 coord_cartesian(ylim=c(0,yupper[dtpos]))+
                 labs(x='Placebo',y="Expected RMSPE  (GDP per capita)")+
                 secaxisy+
                 textadjust)
      
    print(myplot+   theme(legend.position="none"))
    print(as_ggplot(get_legend(myplot+guides(fill=guide_legend(nrow=2))+theme(legend.key.size=unit(2,"lines"),legend.direction="horizontal"))))
    }
    
    


  # ###############################################
  #DISTRIBUTIONAL DIFFERENCES ACROSS DRAWS, WITHIN PLACEBOS, AVERAGING OVER PLACEBOS:
  ###############################################
  if(errorlambda!=0){
  f <- function(x) {
    r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  distdata<-performance[Model%in%c("MASC","PSC-S","PSC-AG","SC","Matching")&lambda==dt&
                          !Placebo%in%drop_placebos,
                        .(MSPE.mean=mean(MSPEAll),
                          MSPE.p5=quantile(MSPEAll,probs=c(0.05)),
                          MSPE.p25=quantile(MSPEAll,probs=c(0.25)),
                          MSPE.p50=quantile(MSPEAll,probs=c(0.50)),
                          MSPE.p75=quantile(MSPEAll,probs=c(0.75)),
                          MSPE.p95=quantile(MSPEAll,probs=c(0.95))),
                          by=.(Model)]

  catdata<- as.data.table(performance[Model%in%c("MASC","PSC-S","PSC-AG","SC","Matching")&lambda==dt&
                Placebo=="Valencia",
              .(MSPE.mean=mean(MSPEAll),
                MSPE.p5=quantile(MSPEAll,probs=c(0.05)),
                MSPE.p25=quantile(MSPEAll,probs=c(0.25)),
                MSPE.p50=quantile(MSPEAll,probs=c(0.50)),
                MSPE.p75=quantile(MSPEAll,probs=c(0.75)),
                MSPE.p95=quantile(MSPEAll,probs=c(0.95))),by=.(Model,Placebo)])

  catdata[,Model:=factor(Model,levels=c("PSC-S","PSC-AG","Matching","MASC","SC"))]
  distdata[,Model:=factor(Model,levels=c("PSC-S","PSC-AG","Matching","MASC","SC"))]
  


  yupper <-c(0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7)*rescale
  ylower <- c(0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02)*rescale
  yuppercat <-c(0.37,0.37)*rescale
  ylowercat <- c(0.14,0.14)*rescale

  print(ggplot(data=distdata)
        +geom_boxplot(aes(x=Model,ymin=MSPE.p25,lower=MSPE.p25,middle=MSPE.p50,upper=MSPE.p75,ymax=MSPE.p75),stat="identity",width=0.08)
        +geom_point(data=distdata[Model%in%c("MASC","PSC-S", "PSC-AG","SC","Matching"),],aes(x=Model,y=MSPE.mean),size=2)
        +coord_flip(ylim=c(ylower[dtpos],yupper[dtpos]))
        +secaxisydist
        +labs(x="",y="RMSPE (GDP per capita)")
        +theme(axis.text.y = element_text(color='black'),
               axis.text.x = element_text(color='black',angle = 0, hjust = 0.5),
               axis.title.x = element_text(color='black'),axis.title.y=element_blank())
  )
  print(ggplot(data=catdata)
        +geom_boxplot(aes(x=Model,ymin=MSPE.p25,lower=MSPE.p25,middle=MSPE.p50,upper=MSPE.p75,ymax=MSPE.p75),stat="identity",width=0.08)
        +geom_point(data=catdata[Model%in%c("MASC","PSC-S", "PSC-AG","SC","Matching"),],aes(x=Model,y=MSPE.mean),size=2)
        +coord_flip(ylim=c(ylowercat[dtpos],yuppercat[dtpos]))
        +secaxisydist
        +labs(x="",y="RMSPE (GDP per capita)")
        +theme(axis.text.y = element_text(color='black'),
               axis.text.x = element_text(color='black',angle = 0, hjust = 0.5),
               axis.title.x = element_text(color='black'),axis.title.y=element_blank())
  )

 }

}

}
setwd(basepath)
performance<-readRDS('empiricalMCadjust_distdata_lambda330.RDS')
performance1<-readRDS('empiricalMCadjust_distdata_lambda300.RDS')
performance1<-performance1[,names(performance1)%in%names(performance),with=FALSE]
performance<-rbind(performance,performance1)
performance<-performance[!(lambda==3.3&Model%in%c("MASC","Pen. SC","SC","Matching")),]
performance[lambda==3&Model%in%c("MASC","Pen. SC","SC","Matching"),lambda:=4.0]
performance[,Model:=gsub(" \\(Covs\\)","",Model)]
performance[lambda==3.3&Model=="Pen. SC",Model:="PSC-S"]
performance[lambda==3.3&Model=="Pen. SC (SC weights)",Model:="PSC-AG"]
performance[lambda==3.0&Model=="Pen. SC",Model:="PSC-S"]
performance[lambda==3.0&Model=="Pen. SC (SC weights)",Model:="PSC-AG"]

for(var in c("MSPE4","MSPEAll","fit","minMSPE4","minMSPEAll")){
  performance[,eval(var):=get(var)*rescale]
}

setwd(outputpath)
pdf("Empirical_Sim_results.pdf")
plot_Emp_MC(performance)
dev.off()

